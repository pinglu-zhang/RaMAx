#!/usr/bin/env bash
set -Eeuo pipefail

usage() {
    cat <<'EOF'
Run an isolated, sequential RaMAx merge-performance A/B comparison.

The baseline executable must be built from commit 3ae2009. The optimized
executable should be built from the current dev-graph-bug worktree.

Usage:
  run_merge_performance_ab.sh \
    --baseline-bin /path/to/ramax-at-3ae2009 \
    --optimized-bin /path/to/ramax-dev-graph-bug \
    --input /path/to/anc12.seqfile \
    --run-root /mnt/d/ramax-ab/anc12-merge-001 \
    [--threads 32] [--sample-seconds 5] \
    [--format maf] [--format gfa] [--format paf] [--format hal] \
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
sample_seconds=5
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
        maf|gfa|paf|hal) ;;
        *)
            printf 'Unsupported output format: %s\n' "$format" >&2
            exit 2
            ;;
    esac
    if [[ "$format" == hal ]]; then
        needs_hal_tools=true
    fi
done

if "$needs_hal_tools"; then
    hal_validate_bin=$(command -v halValidate || true)
    hal_stats_bin=$(command -v halStats || true)
    [[ -x "$hal_validate_bin" ]] || {
        printf '%s\n' 'HAL comparison requires halValidate in PATH.' >&2
        exit 2
    }
    [[ -x "$hal_stats_bin" ]] || {
        printf '%s\n' 'HAL comparison requires halStats in PATH.' >&2
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
exec 9>"$run_root/merge-performance-ab.lock"
if ! /usr/bin/flock -n 9; then
    printf 'Another process holds %s\n' \
        "$run_root/merge-performance-ab.lock" >&2
    exit 1
fi

printf 'label\texecutable\tinput\tthreads\tstarted_at\texit_code\n' \
    >"$run_root/manifest.tsv"
merge_header=$'label\tmerge_index\treference\tstart_timestamp\tcomplete_timestamp\twall_seconds\tinitialization_seconds\ttraversal_seconds\tretirement_seconds\tsampling_rebuild_seconds\tboundaries_scanned\toverlaps_found\tparticipant_conflicts\tmerges_committed\tcigar_units_scanned\tsuffix_forward_scan_steps\tblocks_before\tblocks_after\tblocks_created\tblocks_retired\tblocks_released\tsegments_created\tsegments_released\tmaximum_retirement_batch\ttouched_genome_ends\tstart_rss_kib\tend_rss_kib\tsampled_peak_rss_kib\tprocess_peak_rss_kib'
printf '%s\n' "$merge_header" >"$run_root/merge-summary.tsv"

write_merge_summary() {
    local label=$1
    local console_log=$2
    local output=$3

    printf '%s\n' "$merge_header" >"$output"
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
        function epoch(line, closing, value, whole, millis) {
            closing = index(line, "]")
            if (closing <= 1) return 0
            value = substr(line, 2, closing - 2)
            whole = substr(value, 1, 19)
            gsub(/[-:]/, " ", whole)
            millis = substr(value, 21, 3) + 0
            return mktime(whole) + millis / 1000.0
        }
        function reset_metrics() {
            initialization = traversal = retirement = sampling = "NA"
            boundaries = overlaps = conflicts = committed = "NA"
            cigar_scanned = suffix_steps = blocks_before = blocks_after = "NA"
            blocks_created = blocks_retired = blocks_released = "NA"
            segments_created = segments_released = maximum_batch = "NA"
            touched_ends = start_rss = end_rss = sampled_peak = process_peak = "NA"
            graph_wall = "NA"
        }
        /merge multiple genome graphs for / && $0 !~ / done$/ {
            ++merge_index
            reference = $NF
            start_timestamp = timestamp($0)
            start_epoch = epoch($0)
            reset_metrics()
            next
        }
        /\[stage-memory\] stage=graph-merge event=start/ {
            start_rss = metric("rss_kib")
            next
        }
        /\[graph-merge\]/ && /wall_seconds=/ {
            reference = metric("reference")
            graph_wall = metric("wall_seconds")
            initialization = metric("initialization_seconds")
            traversal = metric("traversal_seconds")
            retirement = metric("retirement_seconds")
            sampling = metric("sampling_rebuild_seconds")
            boundaries = metric("boundaries_scanned")
            overlaps = metric("overlaps_found")
            conflicts = metric("participant_conflicts")
            committed = metric("merges_committed")
            cigar_scanned = metric("cigar_units_scanned")
            suffix_steps = metric("suffix_forward_scan_steps")
            blocks_before = metric("blocks_before")
            blocks_after = metric("blocks_after")
            blocks_created = metric("blocks_created")
            blocks_retired = metric("blocks_retired")
            blocks_released = metric("blocks_released")
            segments_created = metric("segments_created")
            segments_released = metric("segments_released")
            maximum_batch = metric("maximum_retirement_batch")
            touched_ends = metric("touched_genome_ends")
            if (start_rss == "NA") start_rss = metric("start_rss_kib")
            end_rss = metric("end_rss_kib")
            sampled_peak = metric("sampled_peak_rss_kib")
            process_peak = metric("process_peak_rss_kib")
            next
        }
        /\[stage-memory\] stage=graph-merge event=complete/ {
            end_rss = metric("rss_kib")
            process_peak = metric("peak_rss_kib")
            next
        }
        /merge multiple genome graphs for / && / done$/ {
            complete_timestamp = timestamp($0)
            measured_wall = epoch($0) - start_epoch
            wall = graph_wall == "NA" ? sprintf("%.3f", measured_wall) : graph_wall
            print label, merge_index, reference, start_timestamp,
                  complete_timestamp, wall, initialization, traversal,
                  retirement, sampling, boundaries, overlaps, conflicts,
                  committed, cigar_scanned, suffix_steps, blocks_before,
                  blocks_after, blocks_created, blocks_retired,
                  blocks_released, segments_created, segments_released,
                  maximum_batch, touched_ends, start_rss, end_rss,
                  sampled_peak, process_peak
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
    local stage_log="$base/logs/merge-stage.log"
    local summary_log="$base/logs/merge-summary.tsv"
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
        'merge multiple genome graphs for|\[stage-memory\] stage=graph-merge|\[graph-merge\]' \
        "$console_log" >"$stage_log" || true
    write_merge_summary "$label" "$console_log" "$summary_log"
    /usr/bin/tail -n +2 "$summary_log" >>"$run_root/merge-summary.tsv"

    ((exit_code == 0)) || {
        printf '%s exited with status %s; artifacts are preserved under %s\n' \
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
    if [[ "$format" == maf || "$format" == gfa || "$format" == paf ]]; then
        cmp \
            "$run_root/baseline/output/alignment.$format" \
            "$run_root/optimized/output/alignment.$format"
    else
        cmp \
            "$run_root/baseline/output/alignment.hal.stats" \
            "$run_root/optimized/output/alignment.hal.stats"
    fi
done

printf 'A/B artifacts are preserved under %s\n' "$run_root"
printf '%s\n' \
    'Review merge-summary.tsv, rss.tsv, time-v.txt, byte comparisons, and HAL validation before accepting the optimization.'
