#!/usr/bin/env bash
set -Eeuo pipefail

usage() {
    cat <<'EOF'
Run an isolated, sequential seqnum-4 pathology A/B comparison.

The baseline executable must be built from commit 4d204bd. The optimized
executable should be built from the current dev-graph-bug worktree.
This script never verifies a hash and never reuses an existing run root.

Usage:
  run_seqnum4_pathology_ab.sh \
    --baseline-bin /path/to/ramax-at-4d204bd \
    --optimized-bin /path/to/ramax-optimized \
    --input /path/to/seqnum-4.seqfile \
    --run-root /mnt/d/ramax-ab/seqnum-4-001 \
    [--threads 32] [--sample-seconds 5] \
    [--format maf] [--format gfa] [--format paf] [--format hal] \
    [-- RAMAX_COMMON_ARGS...]

The baseline runs first and the optimized build runs second. Every workdir,
partial file, routing table, chunk manifest, stderr log, and final output is
preserved. Run only when no other ramax, wfmash, or minipoa job is active.
EOF
}

baseline_bin=""
optimized_bin=""
input_path=""
run_root=""
threads=32
sample_seconds=5
declare -a formats=()
declare -a common_args=()
hal_validate_bin=""
hal_stats_bin=""
needs_hal=false

while (($# > 0)); do
    case "$1" in
        --baseline-bin)
            baseline_bin=${2:?missing value for --baseline-bin}
            shift 2
            ;;
        --optimized-bin)
            optimized_bin=${2:?missing value for --optimized-bin}
            shift 2
            ;;
        --input)
            input_path=${2:?missing value for --input}
            shift 2
            ;;
        --run-root)
            run_root=${2:?missing value for --run-root}
            shift 2
            ;;
        --threads)
            threads=${2:?missing value for --threads}
            shift 2
            ;;
        --sample-seconds)
            sample_seconds=${2:?missing value for --sample-seconds}
            shift 2
            ;;
        --format)
            formats+=("${2:?missing value for --format}")
            shift 2
            ;;
        --)
            shift
            common_args=("$@")
            break
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            printf 'Unknown option: %s\n' "$1" >&2
            usage >&2
            exit 2
            ;;
    esac
done

[[ -n "$baseline_bin" ]] || {
    printf '%s\n' '--baseline-bin is required' >&2
    exit 2
}
[[ -n "$optimized_bin" ]] || {
    printf '%s\n' '--optimized-bin is required' >&2
    exit 2
}
[[ -n "$input_path" ]] || {
    printf '%s\n' '--input is required' >&2
    exit 2
}
[[ -n "$run_root" ]] || {
    printf '%s\n' '--run-root is required' >&2
    exit 2
}
[[ -x "$baseline_bin" ]] || {
    printf 'Not executable: %s\n' "$baseline_bin" >&2
    exit 2
}
[[ -x "$optimized_bin" ]] || {
    printf 'Not executable: %s\n' "$optimized_bin" >&2
    exit 2
}
[[ -f "$input_path" ]] || {
    printf 'Input does not exist: %s\n' "$input_path" >&2
    exit 2
}
[[ "$threads" =~ ^[1-9][0-9]*$ ]] || {
    printf 'Invalid thread count: %s\n' "$threads" >&2
    exit 2
}
[[ "$sample_seconds" =~ ^[1-9][0-9]*$ ]] || {
    printf 'Invalid sample interval: %s\n' "$sample_seconds" >&2
    exit 2
}
[[ ! -e "$run_root" ]] || {
    printf 'Refusing to reuse existing run root: %s\n' "$run_root" >&2
    exit 1
}
[[ -x /usr/bin/time ]] || {
    printf '%s\n' '/usr/bin/time is required' >&2
    exit 2
}
[[ -x /usr/bin/flock ]] || {
    printf '%s\n' '/usr/bin/flock is required' >&2
    exit 2
}

if ((${#formats[@]} == 0)); then
    formats=(maf gfa paf hal)
fi
for format in "${formats[@]}"; do
    case "$format" in
        maf|gfa|paf) ;;
        hal) needs_hal=true ;;
        *)
            printf 'Unsupported output format: %s\n' "$format" >&2
            exit 2
            ;;
    esac
done
if "$needs_hal"; then
    hal_validate_bin=$(command -v halValidate || true)
    hal_stats_bin=$(command -v halStats || true)
    [[ -x "$hal_validate_bin" && -x "$hal_stats_bin" ]] || {
        printf '%s\n' 'HAL output requires halValidate and halStats in PATH.' >&2
        exit 2
    }
fi

for process_name in ramax wfmash minipoa; do
    if env -u LD_LIBRARY_PATH /usr/bin/pgrep -x "$process_name" \
        >/dev/null 2>&1; then
        printf 'A %s process is already running; do not start a competing A/B.\n' \
            "$process_name" >&2
        exit 1
    fi
done

mkdir -p "$(dirname "$run_root")"
mkdir "$run_root"
exec 9>"$run_root/seqnum4-pathology-ab.lock"
if ! /usr/bin/flock -n 9; then
    printf 'Another process holds %s\n' \
        "$run_root/seqnum4-pathology-ab.lock" >&2
    exit 1
fi

printf 'expected_baseline_commit\t4d204bd\n' >"$run_root/provenance.tsv"
printf 'label\texecutable\tinput\tthreads\tstarted_at\texit_code\n' \
    >"$run_root/manifest.tsv"
printf 'format\tcomparison\n' >"$run_root/comparison.tsv"

sample_process_memory() {
    local pid=$1
    local output=$2
    local timestamp rss_kib peak_rss_kib virtual_kib swap_kib

    printf 'timestamp\tpid\trss_kib\tpeak_rss_kib\tvirtual_kib\tswap_kib\n' \
        >"$output"
    while kill -0 "$pid" 2>/dev/null; do
        timestamp=$(date --iso-8601=seconds)
        rss_kib=0
        peak_rss_kib=0
        virtual_kib=0
        swap_kib=0
        if [[ -r "/proc/$pid/status" ]]; then
            while read -r key value _; do
                case "$key" in
                    VmRSS:) rss_kib=$value ;;
                    VmHWM:) peak_rss_kib=$value ;;
                    VmSize:) virtual_kib=$value ;;
                    VmSwap:) swap_kib=$value ;;
                esac
            done <"/proc/$pid/status"
            printf '%s\t%s\t%s\t%s\t%s\t%s\n' \
                "$timestamp" "$pid" "$rss_kib" "$peak_rss_kib" \
                "$virtual_kib" "$swap_kib" >>"$output"
        fi
        sleep "$sample_seconds"
    done
}

write_hal_summary() {
    local hal_path=$1
    local output=$2
    local genomes genome

    genomes=$("$hal_stats_bin" --genomes "$hal_path")
    {
        printf '%s\n' '[basic]'
        "$hal_stats_bin" "$hal_path"
        printf '%s\n' '[tree]'
        "$hal_stats_bin" --tree "$hal_path"
        printf '%s\n%s\n' '[genomes]' "$genomes"
        for genome in $genomes; do
            printf '[sequenceStats %s]\n' "$genome"
            "$hal_stats_bin" --sequenceStats "$genome" "$hal_path"
            printf '[numSegments %s]\n' "$genome"
            "$hal_stats_bin" --numSegments "$genome" "$hal_path"
        done
    } >"$output"
}

run_one() {
    local label=$1
    local executable=$2
    local base="$run_root/$label"
    local console_log="$base/logs/console.log"
    local time_log="$base/logs/time-v.txt"
    local pid_file="$base/logs/ramax.pid"
    local rss_log="$base/logs/rss.tsv"
    local started_at time_pid ramax_pid sampler_pid exit_code format
    declare -a output_args=()

    mkdir -p "$base/work" "$base/output" "$base/logs"
    for format in "${formats[@]}"; do
        output_args+=(--output "$base/output/alignment.$format")
    done
    started_at=$(date --iso-8601=seconds)
    printf 'Starting %s at %s\n' "$label" "$started_at" | tee "$console_log"
    OMP_NUM_THREADS=$threads /usr/bin/time -v -o "$time_log" \
        /bin/bash -c \
        'printf "%s\n" "$$" >"$1"; shift; exec "$@"' \
        bash "$pid_file" "$executable" \
        --input "$input_path" \
        "${output_args[@]}" \
        --workdir "$base/work" \
        --threads "$threads" \
        "${common_args[@]}" >>"$console_log" 2>&1 &
    time_pid=$!

    for _ in {1..300}; do
        [[ -s "$pid_file" ]] && break
        kill -0 "$time_pid" 2>/dev/null || break
        sleep 0.1
    done
    [[ -s "$pid_file" ]] || {
        wait "$time_pid" || true
        printf 'RaMAx PID was not recorded for %s\n' "$label" >&2
        return 1
    }
    read -r ramax_pid <"$pid_file"
    sample_process_memory "$ramax_pid" "$rss_log" &
    sampler_pid=$!

    set +e
    wait "$time_pid"
    exit_code=$?
    set -e
    wait "$sampler_pid" 2>/dev/null || true
    printf '%s\t%s\t%s\t%s\t%s\t%s\n' \
        "$label" "$executable" "$input_path" "$threads" \
        "$started_at" "$exit_code" >>"$run_root/manifest.tsv"

    env -u LD_LIBRARY_PATH /bin/grep -E \
        '\[wfmash-router\]|\[cluster-extension|\[stage-memory\].*(suffix-array|cluster-extension)' \
        "$console_log" >"$base/logs/pathology-stages.log" || true
    find "$base/work" -type f \
        \( -name routing.tsv -o -name alignment_chunks.tsv \
           -o -name '*stderr*.log' -o -name 'mappings.paf' \
           -o -name 'alignment.paf' \) \
        -print | LC_ALL=C sort >"$base/logs/router-artifacts.txt"

    ((exit_code == 0)) || {
        printf '%s failed with status %s; all partial evidence is preserved in %s\n' \
            "$label" "$exit_code" "$base" >&2
        return "$exit_code"
    }
    for format in "${formats[@]}"; do
        [[ -s "$base/output/alignment.$format" ]] || {
            printf 'Missing or empty %s output for %s\n' "$format" "$label" >&2
            return 1
        }
        if [[ "$format" == hal ]]; then
            "$hal_validate_bin" "$base/output/alignment.hal" \
                >"$base/logs/halValidate.txt" 2>&1
            write_hal_summary "$base/output/alignment.hal" \
                "$base/output/alignment.hal.stats"
        fi
    done
    if kill -0 "$ramax_pid" 2>/dev/null; then
        printf 'The RaMAx process remained after %s; inspect before continuing.\n' \
            "$label" >&2
        return 1
    fi
    local process_name
    for process_name in wfmash minipoa; do
        if env -u LD_LIBRARY_PATH /usr/bin/pgrep -x "$process_name" \
                >/dev/null 2>&1; then
            printf 'A %s process remained after %s; inspect before continuing.\n' \
                "$process_name" "$label" >&2
            return 1
        fi
    done
    printf '%s\n' complete >"$base/.complete"
    printf 'Completed %s at %s\n' "$label" "$(date --iso-8601=seconds)" \
        | tee -a "$console_log"
}

run_one baseline "$baseline_bin"
run_one optimized "$optimized_bin"

for format in "${formats[@]}"; do
    baseline_output="$run_root/baseline/output/alignment.$format"
    optimized_output="$run_root/optimized/output/alignment.$format"
    if [[ "$format" == hal ]]; then
        baseline_output+=".stats"
        optimized_output+=".stats"
    fi
    if cmp -s "$baseline_output" "$optimized_output"; then
        printf '%s\tequal\n' "$format" >>"$run_root/comparison.tsv"
    else
        printf '%s\tdifferent\n' "$format" >>"$run_root/comparison.tsv"
    fi
done

printf 'A/B artifacts are preserved under %s\n' "$run_root"
printf '%s\n' \
    'Review time-v.txt, rss.tsv, routing/chunk diagnostics, comparison.tsv, coverage, and HAL validation before accepting the optimization.'
