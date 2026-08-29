#!/usr/bin/env bash
set -Eeuo pipefail

usage() {
    cat <<'EOF'
Run the isolated RaMAx dev-perf correctness and scaling acceptance sequence.

The dataset manifest is tab-separated and must contain exactly one seqfile for
each scale:

  10<TAB>/data/chr1/seqfile-10.txt
  50<TAB>/data/chr1/seqfile-50.txt
  100<TAB>/data/chr1/seqfile-100.txt
  250<TAB>/data/chr1/seqfile-250.txt
  500<TAB>/data/chr1/seqfile-500.txt
  1000<TAB>/data/chr1/seqfile-1000.txt

Usage:
  run_dev_perf_server_validation.sh \
    --baseline-bin /path/to/ramax-at-4d204bd7cd28 \
    --optimized-bin /path/to/ramax-dev-perf \
    --datasets /path/to/chr1-datasets.tsv \
    --run-root /mnt/d/ramax-dev-perf/validation-001 \
    [--threads 64] [--memory-limit 900GiB] \
    [--memory-ceiling-gib 900] [--sample-seconds 5] \
    [--mode one-round|multi-round] \
    [--format maf] [--format gfa] [--format paf] [--format hal] \
    [-- RAMAX_COMMON_ARGS...]

The baseline and optimized binaries run sequentially for 10 and 50 samples.
Only the optimized binary runs for 100, 250, 500, and 1000 samples. The script
refuses to reuse a run root, never removes experiment data, requires swap to
remain zero for optimized runs, and creates .complete only after validation.
Run it only on an otherwise idle server.
EOF
}

baseline_bin=""
optimized_bin=""
datasets_manifest=""
run_root=""
threads=64
memory_limit="900GiB"
memory_ceiling_gib=900
sample_seconds=5
run_mode="one-round"
declare -a formats=()
declare -a common_args=()

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
        --datasets)
            datasets_manifest=${2:?missing value for --datasets}
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
        --memory-limit)
            memory_limit=${2:?missing value for --memory-limit}
            shift 2
            ;;
        --memory-ceiling-gib)
            memory_ceiling_gib=${2:?missing value for --memory-ceiling-gib}
            shift 2
            ;;
        --sample-seconds)
            sample_seconds=${2:?missing value for --sample-seconds}
            shift 2
            ;;
        --mode)
            run_mode=${2:?missing value for --mode}
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

for required in baseline_bin optimized_bin datasets_manifest run_root; do
    [[ -n "${!required}" ]] || {
        printf '%s is required\n' "--${required//_/-}" >&2
        exit 2
    }
done
[[ -x "$baseline_bin" ]] || {
    printf 'Not executable: %s\n' "$baseline_bin" >&2
    exit 2
}
[[ -x "$optimized_bin" ]] || {
    printf 'Not executable: %s\n' "$optimized_bin" >&2
    exit 2
}
[[ -f "$datasets_manifest" ]] || {
    printf 'Dataset manifest does not exist: %s\n' "$datasets_manifest" >&2
    exit 2
}
[[ "$threads" =~ ^[1-9][0-9]*$ ]] || {
    printf 'Invalid thread count: %s\n' "$threads" >&2
    exit 2
}
[[ "$memory_ceiling_gib" =~ ^[1-9][0-9]*$ ]] || {
    printf 'Invalid memory ceiling: %s\n' "$memory_ceiling_gib" >&2
    exit 2
}
[[ "$sample_seconds" =~ ^[1-9][0-9]*$ ]] || {
    printf 'Invalid sample interval: %s\n' "$sample_seconds" >&2
    exit 2
}
case "$run_mode" in
    one-round|multi-round) ;;
    *)
        printf 'Unsupported mode: %s\n' "$run_mode" >&2
        exit 2
        ;;
esac
[[ ! -e "$run_root" ]] || {
    printf 'Refusing to reuse existing run root: %s\n' "$run_root" >&2
    exit 1
}

for required_command in /usr/bin/time /usr/bin/flock /usr/bin/setsid \
        /usr/bin/pgrep /usr/bin/ps /usr/bin/awk /usr/bin/cmp; do
    [[ -x "$required_command" ]] || {
        printf 'Required executable is unavailable: %s\n' \
            "$required_command" >&2
        exit 2
    }
done

if ((${#formats[@]} == 0)); then
    formats=(maf gfa paf hal)
fi
needs_hal=false
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

hal_validate_bin=""
hal_stats_bin=""
hal_to_maf_bin=""
if "$needs_hal"; then
    hal_validate_bin=$(command -v halValidate || true)
    hal_stats_bin=$(command -v halStats || true)
    hal_to_maf_bin=$(command -v hal2maf || true)
    for hal_tool in hal_validate_bin hal_stats_bin hal_to_maf_bin; do
        [[ -x "${!hal_tool}" ]] || {
            printf 'HAL validation requires %s in PATH.\n' \
                "${hal_tool%_bin}" >&2
            exit 2
        }
    done
fi

declare -A seqfiles=()
while IFS=$'\t' read -r sample_count seqfile extra; do
    [[ -z "$sample_count" || "$sample_count" == \#* ]] && continue
    [[ -z "${extra:-}" ]] || {
        printf 'Manifest row has extra tab-separated fields: %s\n' \
            "$sample_count" >&2
        exit 2
    }
    [[ "$sample_count" =~ ^[1-9][0-9]*$ && -n "$seqfile" ]] || {
        printf 'Invalid manifest row for sample count %s\n' \
            "$sample_count" >&2
        exit 2
    }
    [[ -z "${seqfiles[$sample_count]+present}" ]] || {
        printf 'Duplicate sample count in manifest: %s\n' \
            "$sample_count" >&2
        exit 2
    }
    [[ -f "$seqfile" ]] || {
        printf 'Seqfile does not exist for %s samples: %s\n' \
            "$sample_count" "$seqfile" >&2
        exit 2
    }
    seqfiles[$sample_count]=$seqfile
done <"$datasets_manifest"

declare -a scales=(10 50 100 250 500 1000)
for sample_count in "${scales[@]}"; do
    [[ -n "${seqfiles[$sample_count]+present}" ]] || {
        printf 'Dataset manifest is missing sample count %s\n' \
            "$sample_count" >&2
        exit 2
    }
    observed_samples=$(/usr/bin/awk '
        /^[[:space:]]*$/ || /^[[:space:]]*#/ { next }
        /^[[:space:]]*\(/ { next }
        { ++count }
        END { print count + 0 }
    ' "${seqfiles[$sample_count]}")
    ((observed_samples == sample_count)) || {
        printf 'Seqfile for scale %s contains %s genome mappings: %s\n' \
            "$sample_count" "$observed_samples" \
            "${seqfiles[$sample_count]}" >&2
        exit 2
    }
done

for process_name in ramax wfmash minipoa; do
    if env -u LD_LIBRARY_PATH /usr/bin/pgrep -x "$process_name" \
            >/dev/null 2>&1; then
        printf 'A %s process is already running; refuse competing validation.\n' \
            "$process_name" >&2
        exit 1
    fi
done

mkdir -p "$(dirname "$run_root")"
mkdir "$run_root"
exec 9>"$run_root/dev-perf-validation.lock"
if ! /usr/bin/flock -n 9; then
    printf 'Another process holds %s\n' \
        "$run_root/dev-perf-validation.lock" >&2
    exit 1
fi

printf 'sample_count\tlabel\tmode\texecutable\tseqfile\tthreads\tstarted_at\texit_code\tmax_group_rss_kib\tmax_group_swap_kib\ttime_max_rss_kib\n' \
    >"$run_root/manifest.tsv"
{
    printf 'baseline_version\t'
    "$baseline_bin" --version
    printf 'optimized_version\t'
    "$optimized_bin" --version
    printf 'kernel\t%s\n' "$(uname -srmo)"
    printf 'threads\t%s\n' "$threads"
    printf 'memory_limit\t%s\n' "$memory_limit"
    printf 'memory_ceiling_gib\t%s\n' "$memory_ceiling_gib"
    printf 'run_mode\t%s\n' "$run_mode"
    command -v lscpu >/dev/null 2>&1 && lscpu
    command -v numactl >/dev/null 2>&1 && numactl --hardware
} >"$run_root/environment.txt" 2>&1

declare -a numa_prefix=()
if command -v numactl >/dev/null 2>&1; then
    numa_nodes=$(numactl --hardware 2>/dev/null |
        /usr/bin/awk '/^available:/ { print $2; exit }')
    if [[ "${numa_nodes:-0}" =~ ^[0-9]+$ ]] && ((numa_nodes > 1)); then
        numa_prefix=(numactl --interleave=all)
    fi
fi

sample_process_group() {
    local leader_pid=$1
    local process_group=$2
    local main_pid=$3
    local output=$4
    local timestamp pid key value
    local group_rss group_swap main_rss main_peak

    printf 'timestamp\tprocess_group\tmain_pid\tgroup_rss_kib\tgroup_swap_kib\tmain_rss_kib\tmain_peak_rss_kib\n' \
        >"$output"
    while kill -0 "$leader_pid" 2>/dev/null; do
        timestamp=$(date --iso-8601=seconds)
        group_rss=0
        group_swap=0
        main_rss=0
        main_peak=0
        while read -r pid; do
            [[ -r "/proc/$pid/status" ]] || continue
            while read -r key value _; do
                case "$key" in
                    VmRSS:) group_rss=$((group_rss + value)) ;;
                    VmSwap:) group_swap=$((group_swap + value)) ;;
                esac
            done <"/proc/$pid/status"
        done < <(/usr/bin/pgrep -g "$process_group" || true)
        if [[ -r "/proc/$main_pid/status" ]]; then
            while read -r key value _; do
                case "$key" in
                    VmRSS:) main_rss=$value ;;
                    VmHWM:) main_peak=$value ;;
                esac
            done <"/proc/$main_pid/status"
        fi
        printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
            "$timestamp" "$process_group" "$main_pid" "$group_rss" \
            "$group_swap" "$main_rss" "$main_peak" >>"$output"
        sleep "$sample_seconds"
    done
}

write_hal_semantic_summary() {
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
    "$hal_to_maf_bin" "$hal_path" "$output.maf"
}

run_one() {
    local sample_count=$1
    local label=$2
    local executable=$3
    local is_optimized=$4
    local seqfile=${seqfiles[$sample_count]}
    local base="$run_root/$run_mode/$sample_count/$label"
    local console_log="$base/logs/console.log"
    local time_log="$base/logs/time-v.txt"
    local rss_log="$base/logs/process-group-rss.tsv"
    local pid_file="$base/logs/ramax.pid"
    local started_at leader_pid main_pid process_group sampler_pid
    local exit_code max_group_rss max_group_swap time_max_rss format
    local ceiling_kib=$((memory_ceiling_gib * 1024 * 1024))
    declare -a output_args=()
    declare -a mode_args=()
    declare -a resource_args=()

    mkdir -p "$base/work" "$base/output" "$base/logs" "$base/temp"
    for format in "${formats[@]}"; do
        output_args+=(-o "$base/output/alignment.$format")
    done
    [[ "$run_mode" == one-round ]] && mode_args+=(--one-round)
    if [[ "$is_optimized" == true ]]; then
        resource_args+=(--memory-limit "$memory_limit" --temp-dir "$base/temp")
    fi

    started_at=$(date --iso-8601=seconds)
    printf 'Starting %s %s-sample %s at %s\n' \
        "$label" "$sample_count" "$run_mode" "$started_at" \
        >"$console_log"

    OMP_NUM_THREADS=$threads \
    OMP_PLACES=${OMP_PLACES:-cores} \
    OMP_PROC_BIND=${OMP_PROC_BIND:-spread} \
        /usr/bin/setsid /usr/bin/time -v -o "$time_log" \
        /bin/bash -c \
        'printf "%s\n" "$$" >"$1"; shift; exec "$@"' \
        bash "$pid_file" "${numa_prefix[@]}" "$executable" \
        -i "$seqfile" \
        "${output_args[@]}" \
        -w "$base/work" \
        -t "$threads" \
        "${mode_args[@]}" \
        "${resource_args[@]}" \
        "${common_args[@]}" >>"$console_log" 2>&1 &
    leader_pid=$!

    for _ in {1..300}; do
        [[ -s "$pid_file" ]] && break
        kill -0 "$leader_pid" 2>/dev/null || break
        sleep 0.1
    done
    [[ -s "$pid_file" ]] || {
        wait "$leader_pid" || true
        printf 'RaMAx PID was not recorded for %s/%s.\n' \
            "$sample_count" "$label" >&2
        return 1
    }
    read -r main_pid <"$pid_file"
    process_group=$(/usr/bin/ps -o pgid= -p "$leader_pid" |
        /usr/bin/awk '{ print $1 }')
    [[ -n "$process_group" ]] || {
        printf 'Process group was not recorded for %s/%s.\n' \
            "$sample_count" "$label" >&2
        return 1
    }
    sample_process_group \
        "$leader_pid" "$process_group" "$main_pid" "$rss_log" &
    sampler_pid=$!

    set +e
    wait "$leader_pid"
    exit_code=$?
    set -e
    wait "$sampler_pid" 2>/dev/null || true

    max_group_rss=$(/usr/bin/awk -F '\t' \
        'NR > 1 && $4 > maximum { maximum = $4 } END { print maximum + 0 }' \
        "$rss_log")
    max_group_swap=$(/usr/bin/awk -F '\t' \
        'NR > 1 && $5 > maximum { maximum = $5 } END { print maximum + 0 }' \
        "$rss_log")
    time_max_rss=$(/usr/bin/awk -F ':' \
        '/Maximum resident set size/ { gsub(/^[[:space:]]+/, "", $2); print $2 }' \
        "$time_log")
    time_max_rss=${time_max_rss:-0}
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "$sample_count" "$label" "$run_mode" "$executable" "$seqfile" \
        "$threads" "$started_at" "$exit_code" "$max_group_rss" \
        "$max_group_swap" "$time_max_rss" >>"$run_root/manifest.tsv"
    /usr/bin/grep -E \
        '\[resource-budget\]|\[spill\]|\[stage-memory\]|\[graph-merge\]|\[extend-ref-nodes\]|\[mask-journal\]' \
        "$console_log" >"$base/logs/stages.log" || true

    ((exit_code == 0)) || {
        printf '%s/%s exited with status %s; all artifacts are preserved.\n' \
            "$sample_count" "$label" "$exit_code" >&2
        return "$exit_code"
    }
    for format in "${formats[@]}"; do
        [[ -s "$base/output/alignment.$format" ]] || {
            printf 'Missing %s output for %s/%s.\n' \
                "$format" "$sample_count" "$label" >&2
            return 1
        }
        if [[ "$format" == hal ]]; then
            "$hal_validate_bin" "$base/output/alignment.hal"
            write_hal_semantic_summary \
                "$base/output/alignment.hal" \
                "$base/output/alignment.hal.stats"
        fi
    done
    if [[ "$is_optimized" == true ]]; then
        ((max_group_rss < ceiling_kib)) || {
            printf 'Optimized %s-sample group RSS %s KiB reached ceiling %s KiB.\n' \
                "$sample_count" "$max_group_rss" "$ceiling_kib" >&2
            return 1
        }
        ((time_max_rss < ceiling_kib)) || {
            printf 'Optimized %s-sample time-v RSS %s KiB reached ceiling %s KiB.\n' \
                "$sample_count" "$time_max_rss" "$ceiling_kib" >&2
            return 1
        }
        ((max_group_swap == 0)) || {
            printf 'Optimized %s-sample run used %s KiB swap.\n' \
                "$sample_count" "$max_group_swap" >&2
            return 1
        }
    fi
    printf '%s\n' "$(date --iso-8601=seconds)" >"$base/.complete"
}

compare_pair() {
    local sample_count=$1
    local baseline="$run_root/$run_mode/$sample_count/baseline/output"
    local optimized="$run_root/$run_mode/$sample_count/optimized/output"
    local format

    for format in "${formats[@]}"; do
        case "$format" in
            maf|gfa|paf)
                /usr/bin/cmp \
                    "$baseline/alignment.$format" \
                    "$optimized/alignment.$format"
                ;;
            hal)
                /usr/bin/cmp \
                    "$baseline/alignment.hal.stats" \
                    "$optimized/alignment.hal.stats"
                /usr/bin/cmp \
                    "$baseline/alignment.hal.stats.maf" \
                    "$optimized/alignment.hal.stats.maf"
                ;;
        esac
    done
    printf '%s\n' "$(date --iso-8601=seconds)" \
        >"$run_root/$run_mode/$sample_count/.equivalent"
}

for sample_count in 10 50; do
    run_one "$sample_count" baseline "$baseline_bin" false
    run_one "$sample_count" optimized "$optimized_bin" true
    compare_pair "$sample_count"
done

for sample_count in 100 250 500 1000; do
    run_one "$sample_count" optimized "$optimized_bin" true
done

printf '%s\n' "$(date --iso-8601=seconds)" >"$run_root/.complete"
printf 'Validation completed; all artifacts are preserved under %s\n' \
    "$run_root"
