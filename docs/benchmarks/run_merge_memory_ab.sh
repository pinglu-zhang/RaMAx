#!/usr/bin/env bash
set -Eeuo pipefail

usage() {
    cat <<'EOF'
Run an isolated, sequential RaMAx merge-memory A/B comparison.

Usage:
  run_merge_memory_ab.sh \
    --baseline-bin /path/to/ramax-at-448de1e \
    --optimized-bin /path/to/ramax-with-merge-memory-fix \
    --input /path/to/Bos.seqfile \
    --run-root /mnt/d/ramax-ab/bos-merge-memory-001 \
    [--threads 32] [--sample-seconds 60] \
    [--format maf] [--format gfa] [--format hal] \
    [-- RAMAX_COMMON_ARGS...]

The script runs baseline first and optimized second. It refuses to reuse an
existing run root, preserves every work/output/log file, and never removes
experiment data. Run it only when no other RaMAx or alignment job is active.
EOF
}

baseline_bin=""
optimized_bin=""
input_path=""
run_root=""
threads=32
sample_seconds=60
declare -a formats=()
declare -a common_args=()
needs_hal_tools=false
hal_validate_bin=""
hal_stats_bin=""

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

[[ -n "$baseline_bin" ]] || { printf '%s\n' '--baseline-bin is required' >&2; exit 2; }
[[ -n "$optimized_bin" ]] || { printf '%s\n' '--optimized-bin is required' >&2; exit 2; }
[[ -n "$input_path" ]] || { printf '%s\n' '--input is required' >&2; exit 2; }
[[ -n "$run_root" ]] || { printf '%s\n' '--run-root is required' >&2; exit 2; }
[[ -x "$baseline_bin" ]] || { printf 'Not executable: %s\n' "$baseline_bin" >&2; exit 2; }
[[ -x "$optimized_bin" ]] || { printf 'Not executable: %s\n' "$optimized_bin" >&2; exit 2; }
[[ -f "$input_path" ]] || { printf 'Input does not exist: %s\n' "$input_path" >&2; exit 2; }
[[ "$threads" =~ ^[1-9][0-9]*$ ]] || { printf 'Invalid thread count: %s\n' "$threads" >&2; exit 2; }
[[ "$sample_seconds" =~ ^[1-9][0-9]*$ ]] || { printf 'Invalid sample interval: %s\n' "$sample_seconds" >&2; exit 2; }
[[ ! -e "$run_root" ]] || { printf 'Refusing to reuse existing run root: %s\n' "$run_root" >&2; exit 1; }
[[ -x /usr/bin/time ]] || { printf '%s\n' '/usr/bin/time is required' >&2; exit 2; }
[[ -x /usr/bin/flock ]] || { printf '%s\n' '/usr/bin/flock is required' >&2; exit 2; }

if ((${#formats[@]} == 0)); then
    formats=(maf gfa hal)
fi
for format in "${formats[@]}"; do
    case "$format" in
        maf|gfa|hal|paf) ;;
        *) printf 'Unsupported output format: %s\n' "$format" >&2; exit 2 ;;
    esac
    if [[ "$format" == hal ]]; then
        needs_hal_tools=true
    fi
done

if "$needs_hal_tools"; then
    hal_validate_bin=$(command -v halValidate || true)
    hal_stats_bin=$(command -v halStats || true)
    [[ -x "$hal_validate_bin" ]] || {
        printf '%s\n' 'HAL output comparison requires halValidate in PATH.' >&2
        exit 2
    }
    [[ -x "$hal_stats_bin" ]] || {
        printf '%s\n' 'HAL output comparison requires halStats in PATH.' >&2
        exit 2
    }
fi

for process_name in ramax wfmash minipoa; do
    if env -u LD_LIBRARY_PATH /usr/bin/pgrep -x "$process_name" >/dev/null 2>&1; then
        printf 'A %s process is already running; do not start a competing A/B.\n' \
            "$process_name" >&2
        exit 1
    fi
done

mkdir -p "$(dirname "$run_root")"
mkdir "$run_root"
exec 9>"$run_root/merge-memory-ab.lock"
if ! /usr/bin/flock -n 9; then
    printf 'Another process holds %s\n' "$run_root/merge-memory-ab.lock" >&2
    exit 1
fi

printf 'label\texecutable\tinput\tthreads\tstarted_at\texit_code\n' \
    >"$run_root/manifest.tsv"
printf '%s\n' \
    $'label\tmerge_index\tstart_timestamp\tcomplete_timestamp\tmerge_wall_seconds\tstart_rss_kib\tmerge_sampled_peak_rss_kib\tcomplete_rss_kib\tmerge_increment_peak_kib\tblocks_created\tblocks_retired\tblocks_released\tsegments_released\tcigar_operations_released\tpending_retired_ids\ttouched_genome_ends' \
    >"$run_root/merge-summary.tsv"

write_merge_summary() {
    local label=$1
    local console_log=$2
    local output=$3

    printf '%s\n' \
        $'label\tmerge_index\tstart_timestamp\tcomplete_timestamp\tmerge_wall_seconds\tstart_rss_kib\tmerge_sampled_peak_rss_kib\tcomplete_rss_kib\tmerge_increment_peak_kib\tblocks_created\tblocks_retired\tblocks_released\tsegments_released\tcigar_operations_released\tpending_retired_ids\ttouched_genome_ends' \
        >"$output"
    /usr/bin/awk -v label="$label" '
        BEGIN { OFS = "\t" }
        function metric(name, i, prefix) {
            prefix = name "="
            for (i = 1; i <= NF; ++i) {
                if (index($i, prefix) == 1) {
                    return substr($i, length(prefix) + 1)
                }
            }
            return "NA"
        }
        function timestamp(line, closing) {
            closing = index(line, "]")
            return closing > 1 ? substr(line, 2, closing - 2) : "NA"
        }
        /\[stage-memory\] stage=graph-merge event=start/ {
            ++merge_index
            start_timestamp = timestamp($0)
            start_rss = metric("rss_kib")
            merge_wall = peak_rss = blocks_created = blocks_retired = "NA"
            blocks_released = segments_released = cigar_released = "NA"
            pending_ids = touched_ends = "NA"
            next
        }
        /\[merge-memory\]/ {
            merge_wall = metric("merge-wall-seconds")
            peak_rss = metric("merge-sampled-peak-rss-kib")
            blocks_created = metric("blocks-created")
            blocks_retired = metric("blocks-retired")
            blocks_released = metric("blocks-released")
            segments_released = metric("segments-released")
            cigar_released = metric("cigar-operations-released")
            pending_ids = metric("pending-retired-ids")
            touched_ends = metric("touched-genome-ends")
            next
        }
        /\[stage-memory\] stage=graph-merge event=complete/ {
            complete_timestamp = timestamp($0)
            complete_rss = metric("rss_kib")
            increment = "NA"
            if (start_rss ~ /^[0-9]+$/ && peak_rss ~ /^[0-9]+$/) {
                increment = peak_rss - start_rss
            }
            print label, merge_index, start_timestamp, complete_timestamp,
                  merge_wall, start_rss, peak_rss, complete_rss, increment,
                  blocks_created, blocks_retired, blocks_released,
                  segments_released, cigar_released, pending_ids, touched_ends
        }
    ' "$console_log" >>"$output"
}

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

write_hal_semantic_summary() {
    local hal_stats_bin=$1
    local hal_path=$2
    local output=$3
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
    local stage_log="$base/logs/merge-stage-memory.log"
    local summary_log="$base/logs/merge-summary.tsv"
    local started_at
    local time_pid ramax_pid sampler_pid exit_code
    declare -a output_args=()
    local format

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
        '\[stage-memory\] stage=graph-merge|\[merge-memory\]' \
        "$console_log" >"$stage_log" || true
    write_merge_summary "$label" "$console_log" "$summary_log"
    /usr/bin/tail -n +2 "$summary_log" >>"$run_root/merge-summary.tsv"

    ((exit_code == 0)) || {
        printf '%s exited with status %s; artifacts were preserved under %s\n' \
            "$label" "$exit_code" "$base" >&2
        return "$exit_code"
    }
    for format in "${formats[@]}"; do
        [[ -s "$base/output/alignment.$format" ]] || {
            printf 'Missing %s output for %s\n' "$format" "$label" >&2
            return 1
        }
        if [[ "$format" == hal ]]; then
            "$hal_validate_bin" "$base/output/alignment.hal"
            write_hal_semantic_summary \
                "$hal_stats_bin" \
                "$base/output/alignment.hal" \
                "$base/output/alignment.hal.stats"
        fi
    done
    printf 'Completed %s at %s\n' "$label" "$(date --iso-8601=seconds)" \
        | tee -a "$console_log"
}

run_one baseline "$baseline_bin"
run_one optimized "$optimized_bin"

for format in "${formats[@]}"; do
    sha256sum \
        "$run_root/baseline/output/alignment.$format" \
        "$run_root/optimized/output/alignment.$format" \
        >"$run_root/$format.sha256"
    if [[ "$format" == maf || "$format" == gfa || "$format" == paf ]]; then
        cmp \
            "$run_root/baseline/output/alignment.$format" \
            "$run_root/optimized/output/alignment.$format"
    elif [[ "$format" == hal &&
            -s "$run_root/baseline/output/alignment.hal.stats" &&
            -s "$run_root/optimized/output/alignment.hal.stats" ]]; then
        cmp \
            "$run_root/baseline/output/alignment.hal.stats" \
            "$run_root/optimized/output/alignment.hal.stats"
    fi
done

printf 'A/B artifacts are preserved under %s\n' "$run_root"
printf '%s\n' \
    'Compare rss.tsv and merge-stage-memory.log; overall time-v peak may be dominated by suffix-array construction.'
