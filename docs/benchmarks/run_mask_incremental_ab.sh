#!/usr/bin/env bash
set -Eeuo pipefail

usage() {
    cat <<'EOF'
Run an isolated, sequential RaMAx aligned-mask A/B comparison.

The baseline executable must be built from commit f34923c. The optimized
executable should be built from the current dev-graph-improve worktree.

Usage:
  run_mask_incremental_ab.sh \
    --baseline-bin /path/to/ramax-at-f34923c \
    --optimized-bin /path/to/ramax-dev-graph-improve \
    --input /path/to/dataset.seqfile \
    --run-root /mnt/d/ramax-ab/mask-incremental-001 \
    [--threads 32] [--sample-seconds 5] \
    [--format maf] [--format gfa] [--format paf] [--format hal] \
    [-- RAMAX_COMMON_ARGS...]

The script runs baseline first and optimized second. It refuses to reuse an
existing run root, preserves all work/output/log files, and never removes
experiment data. Run it only when no competing alignment job is active.
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
exec 9>"$run_root/mask-incremental-ab.lock"
if ! /usr/bin/flock -n 9; then
    printf 'Another process holds %s\n' \
        "$run_root/mask-incremental-ab.lock" >&2
    exit 1
fi

printf 'label\texecutable\tinput\tthreads\tstarted_at\texit_code\n' \
    >"$run_root/manifest.tsv"
mask_header=$'label\tmask_index\treference\tstart_timestamp\tcomplete_timestamp\twall_seconds\tscanned_blocks\tscanned_segments\tchanged_candidates\tappended\tunchanged_skipped\tmissing_sequences\ttouched_species\ttouched_chromosomes\tincoming_intervals\tnormalized_delta_intervals\tprevious_intervals\tfinal_intervals\tsort_seconds\tmerge_seconds\tmetadata_seconds\tspecies_wall_seconds\tstart_rss_kib\tend_rss_kib\tprocess_peak_rss_kib'
printf '%s\n' "$mask_header" >"$run_root/mask-summary.tsv"

write_mask_summary() {
    local label=$1
    local console_log=$2
    local output=$3

    printf '%s\n' "$mask_header" >"$output"
    /usr/bin/awk -v label="$label" '
        BEGIN { OFS = "\t"; in_stage = 0 }
        function metric(name, i, prefix) {
            prefix = name "="
            for (i = 1; i <= NF; ++i) {
                if (index($i, prefix) == 1) {
                    return substr($i, length(prefix) + 1)
                }
            }
            return "NA"
        }
        function number_or_zero(value) {
            return value == "NA" ? 0 : value + 0
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
            scanned_blocks = scanned_segments = changed = appended = "NA"
            unchanged = missing = touched_species = touched_chr = "NA"
            incoming = normalized = previous = final = 0
            sort_seconds = merge_seconds = metadata_seconds = species_wall = 0
            finalize_seen = 0
            start_rss = end_rss = process_peak = "NA"
        }
        /Adding aligned regions as mask intervals for / {
            ++mask_index
            reference = $NF
            start_timestamp = timestamp($0)
            start_epoch = epoch($0)
            reset_metrics()
            in_stage = 1
            next
        }
        in_stage && /\[stage-memory\] stage=aligned-mask event=start/ {
            start_rss = metric("rss_kib")
            if (process_peak == "NA") process_peak = metric("peak_rss_kib")
            next
        }
        in_stage && /\[mask-finalize\]/ {
            finalize_seen = 1
            incoming += number_or_zero(metric("incoming_intervals"))
            normalized += number_or_zero(metric("normalized_delta_intervals"))
            previous += number_or_zero(metric("previous_intervals"))
            final += number_or_zero(metric("final_intervals"))
            sort_seconds += number_or_zero(metric("sort_seconds"))
            merge_seconds += number_or_zero(metric("merge_seconds"))
            metadata_seconds += number_or_zero(metric("metadata_seconds"))
            species_wall += number_or_zero(metric("wall_seconds"))
            next
        }
        in_stage && /\[mask-journal\]/ {
            scanned_blocks = metric("scanned_blocks")
            scanned_segments = metric("scanned_segments")
            changed = metric("changed_candidates")
            appended = metric("appended")
            unchanged = metric("unchanged_skipped")
            missing = metric("missing_sequences")
            touched_species = metric("touched_species")
            touched_chr = metric("touched_chromosomes")
            next
        }
        in_stage && /\[stage-memory\] stage=aligned-mask event=complete/ {
            end_rss = metric("rss_kib")
            process_peak = metric("peak_rss_kib")
            next
        }
        in_stage && (/Successfully added mask intervals for round with reference / ||
                     /Mask interval addition completed for all species/) {
            complete_timestamp = timestamp($0)
            wall = sprintf("%.3f", epoch($0) - start_epoch)
            if (!finalize_seen) {
                incoming = normalized = previous = final = "NA"
                sort_seconds = merge_seconds = metadata_seconds = species_wall = "NA"
            }
            print label, mask_index, reference, start_timestamp,
                  complete_timestamp, wall, scanned_blocks, scanned_segments,
                  changed, appended, unchanged, missing, touched_species,
                  touched_chr, incoming, normalized, previous, final,
                  sort_seconds, merge_seconds, metadata_seconds, species_wall,
                  start_rss, end_rss, process_peak
            in_stage = 0
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
    local stage_log="$base/logs/mask-stage.log"
    local summary_log="$base/logs/mask-summary.tsv"
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
        'Adding aligned regions as mask intervals|Successfully added mask intervals|Mask interval addition completed|\[stage-memory\] stage=aligned-mask|\[mask-journal\]|\[mask-finalize\]' \
        "$console_log" >"$stage_log" || true
    write_mask_summary "$label" "$console_log" "$summary_log"
    /usr/bin/tail -n +2 "$summary_log" >>"$run_root/mask-summary.tsv"

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

date --iso-8601=seconds >"$run_root/.complete"
printf 'A/B artifacts are preserved under %s\n' "$run_root"
printf '%s\n' \
    'Review mask-summary.tsv, rss.tsv, time-v.txt, byte comparisons, and HAL validation before accepting the optimization.'
