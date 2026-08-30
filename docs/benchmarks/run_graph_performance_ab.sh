#!/usr/bin/env bash
set -euo pipefail

usage() {
    cat <<'EOF'
Run an isolated, sequential RaMAx graph-performance A/B comparison.

The baseline executable must be built from commit 448de1e. The optimized
executable should be built from the current dev-graph-bug worktree.

Usage:
  run_graph_performance_ab.sh \
    --baseline-bin /path/to/ramax-at-448de1e \
    --optimized-bin /path/to/ramax-dev-graph-bug \
    --input /path/to/anc12.seqfile \
    --run-root /mnt/d/ramax-ab/anc12 \
    [--threads 32] [--format maf] [--format gfa] \
    [-- RAMAX_COMMON_ARGS...]

The script runs baseline first and optimized second. It never removes or
overwrites an existing work directory, output directory, or log. Both runs
receive the same input, thread count, formats, and arguments after "--".
EOF
}

baseline_bin=""
optimized_bin=""
input_path=""
run_root=""
threads=32
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
command -v /usr/bin/time >/dev/null || {
    printf '%s\n' '/usr/bin/time is required' >&2
    exit 2
}
command -v flock >/dev/null || {
    printf '%s\n' 'flock is required' >&2
    exit 2
}

if ((${#formats[@]} == 0)); then
    formats=(maf gfa)
fi
for format in "${formats[@]}"; do
    case "$format" in
        maf|hal|paf|gfa) ;;
        *)
            printf 'Unsupported output format: %s\n' "$format" >&2
            exit 2
            ;;
    esac
done

mkdir -p "$run_root"
exec 9>"$run_root/ramax-ab.lock"
if ! flock -n 9; then
    printf 'Another A/B run holds the lock: %s\n' \
        "$run_root/ramax-ab.lock" >&2
    exit 1
fi

if pgrep -x ramax >/dev/null 2>&1; then
    printf '%s\n' \
        'A ramax process is already running; stop competing jobs before A/B.' >&2
    exit 1
fi

prepare_run() {
    local label=$1
    local base="$run_root/$label"
    if [[ -e "$base/work" || -e "$base/output" || -e "$base/logs" ]]; then
        printf 'Refusing to reuse existing A/B path: %s\n' "$base" >&2
        exit 1
    fi
    mkdir -p "$base/work" "$base/output" "$base/logs"
}

run_one() {
    local label=$1
    local executable=$2
    local base="$run_root/$label"
    local console_log="$base/logs/console.log"
    local time_log="$base/logs/time-v.txt"
    local stage_log="$base/logs/stage-memory.log"
    declare -a output_args=()
    local format
    for format in "${formats[@]}"; do
        output_args+=(--output "$base/output/alignment.$format")
    done

    printf 'Starting %s at %s\n' \
        "$label" "$(date --iso-8601=seconds)" | tee "$console_log"
    OMP_NUM_THREADS=$threads \
        /usr/bin/time -v -o "$time_log" \
        "$executable" \
        --input "$input_path" \
        "${output_args[@]}" \
        --workdir "$base/work" \
        --threads "$threads" \
        "${common_args[@]}" 2>&1 | tee -a "$console_log"
    grep -F '[stage-memory]' "$console_log" >"$stage_log" || true

    for format in "${formats[@]}"; do
        test -s "$base/output/alignment.$format"
        if [[ "$format" == hal ]] && command -v halValidate >/dev/null 2>&1; then
            halValidate "$base/output/alignment.hal"
        fi
    done
    printf 'Completed %s at %s\n' \
        "$label" "$(date --iso-8601=seconds)" | tee -a "$console_log"
}

prepare_run baseline
prepare_run optimized
run_one baseline "$baseline_bin"
run_one optimized "$optimized_bin"

printf 'A/B artifacts are preserved under %s\n' "$run_root"
printf '%s\n' \
    'Compare time-v.txt, stage-memory.log, and validated outputs before accepting the optimization.'
