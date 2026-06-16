#!/usr/bin/env bash
set -uo pipefail

# Batch helper for ICIC + PA correction jobs.
#
# Modes:
#   submit  Stage charge data in q_0/scf, run step 1, and submit step 2 jobs.
#   collect Run step 3 in q_0/scf for corrections whose potential jobs finished.
#   status  Print per-directory status only.
#
# Options:
#   --force Ignore the local correction.submitted / correction.done markers.

MODE="${1:-submit}"
FORCE=0
FILTERS=()

shift_mode=0
if [ "${1:-}" = "submit" ] || [ "${1:-}" = "collect" ] || [ "${1:-}" = "status" ]; then
    shift_mode=1
fi
[ "$shift_mode" -eq 1 ] && shift

while [ "$#" -gt 0 ]; do
    case "$1" in
        --force)
            FORCE=1
            shift
            ;;
        --idx|--defect|--only)
            if [ "$#" -lt 2 ]; then
                echo "ERROR: $1 needs a value" >&2
                exit 1
            fi
            IFS=',' read -ra parts <<< "$2"
            for part in "${parts[@]}"; do
                [ -n "$part" ] && FILTERS+=("$part")
            done
            shift 2
            ;;
        --idx=*|--defect=*|--only=*)
            value="${1#*=}"
            IFS=',' read -ra parts <<< "$value"
            for part in "${parts[@]}"; do
                [ -n "$part" ] && FILTERS+=("$part")
            done
            shift
            ;;
        *)
            IFS=',' read -ra parts <<< "$1"
            for part in "${parts[@]}"; do
                [ -n "$part" ] && FILTERS+=("$part")
            done
            shift
            ;;
    esac
done

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

find_project_root() {
    local d="$SCRIPT_DIR"
    while true; do
        if [ -d "$d/calculate" ] || [ -f "$d/config.yaml" ] || [ -f "$d/input" ]; then
            printf '%s\n' "$d"
            return 0
        fi
        local p
        p="$(dirname "$d")"
        if [ "$p" = "$d" ]; then
            printf '%s\n' "$SCRIPT_DIR"
            return 0
        fi
        d="$p"
    done
}

PROJECT_ROOT="$(find_project_root)"

find_tool_dir() {
    for d in "$SCRIPT_DIR" "$PROJECT_ROOT" "$PROJECT_ROOT/script"; do
        if [ -f "$d/1_get_rho.sh" ] && [ -f "$d/2_coulomb_integral.sh" ] && [ -f "$d/3_get_results.sh" ]; then
            printf '%s\n' "$d"
            return 0
        fi
    done
    return 1
}

TOOL_DIR="$(find_tool_dir || true)"
if [ -z "$TOOL_DIR" ]; then
    echo "ERROR: cannot find 1_get_rho.sh / 2_coulomb_integral.sh / 3_get_results.sh" >&2
    exit 1
fi

if [ ! -d "$PROJECT_ROOT/calculate" ]; then
    echo "ERROR: calculate/ not found under project root: $PROJECT_ROOT" >&2
    exit 1
fi

case "$MODE" in
    submit|collect|status) ;;
    *)
        echo "Usage: bash batch_submit_corrections.sh [submit|collect|status] [--force]" >&2
        exit 1
        ;;
esac

defect_index() {
    local name="$1"
    if [[ "$name" =~ ^([0-9]+)(_|$) ]]; then
        printf '%s\n' "${BASH_REMATCH[1]}"
    else
        printf '\n'
    fi
}

matches_filter() {
    local name="$1"
    local idx
    idx="$(defect_index "$name")"
    if [ "${#FILTERS[@]}" -eq 0 ]; then
        return 0
    fi
    local f
    for f in "${FILTERS[@]}"; do
        if [ "$f" = "$name" ] || [ "$f" = "$idx" ]; then
            return 0
        fi
    done
    return 1
}

charge_folder_from_defect_input() {
    local defect_input="$1"
    local q1 sign qabs
    q1="$(awk '/^charged[[:space:]]+state/{print $3; exit}' "$defect_input")"
    if [ -z "$q1" ]; then
        return 1
    fi
    sign="${q1:0:1}"
    qabs="${q1:1}"
    printf 'image-corr_%s%s\n' "$sign" "$qabs"
}

charge_from_input_name() {
    local name="$1"
    case "$name" in
        defect_+*.input)
            printf '+%s\n' "${name#defect_+}" | sed 's/\.input$//'
            ;;
        defect_-*.input)
            printf -- '-%s\n' "${name#defect_-}" | sed 's/\.input$//'
            ;;
        defect_p*.input)
            printf '+%s\n' "${name#defect_p}" | sed 's/\.input$//'
            ;;
        defect_m*.input)
            printf -- '-%s\n' "${name#defect_m}" | sed 's/\.input$//'
            ;;
        *)
            printf '\n'
            ;;
    esac
}

canonical_defect_input() {
    local input_scf_dir="$1"
    local q_state="$2"
    local qabs="${q_state#+}"
    qabs="${qabs#-}"
    if [[ "$q_state" == -* ]]; then
        printf '%s/defect_-%s.input\n' "$input_scf_dir" "$qabs"
    else
        printf '%s/defect_+%s.input\n' "$input_scf_dir" "$qabs"
    fi
}

ensure_stable_defect_input() {
    local source_defect_input="$1"
    local q_state="$2"
    local input_scf_dir
    local stable_input
    input_scf_dir="$(dirname "$source_defect_input")"
    stable_input="$(canonical_defect_input "$input_scf_dir" "$q_state")"
    if [ "$source_defect_input" != "$stable_input" ]; then
        cp "$source_defect_input" "$stable_input"
    fi
    printf '%s\n' "$stable_input"
}

has_ecoul_report() {
    local report="$1"
    [ -f "$report" ] && grep -q "E_Coul(eV)" "$report"
}

is_finished() {
    local scf_dir="$1"
    local corr_dir="$2"
    has_ecoul_report "$scf_dir/$corr_dir/REPORT" && has_ecoul_report "$scf_dir/$corr_dir/REPORT.0"
}

copy_if_exists() {
    local src="$1"
    local dst="$2"
    if [ -f "$src" ]; then
        cp "$src" "$dst"
    fi
}

prepare_local_helpers() {
    local scf_dir="$1"
    copy_if_exists "$TOOL_DIR/generate_occ.py" "$scf_dir/generate_occ.py"
    copy_if_exists "$TOOL_DIR/generate_defect_rho.sh" "$scf_dir/generate_defect_rho.sh"
    copy_if_exists "$TOOL_DIR/vatom.py" "$scf_dir/vatom.py"
}

prepare_occ_files() {
    local charged_scf_dir="$1"
    local run_dir="$2"
    local source_defect_input="$3"
    local run_defect_input="$run_dir/defect.input"
    local neutral_dir
    neutral_dir="$(awk '/^neutral[[:space:]]/{print $2; exit}' "$source_defect_input")"
    if [ -z "$neutral_dir" ]; then
        neutral_dir="$run_dir"
    fi

    if [ ! -f "$neutral_dir/OUT.OCC" ]; then
        echo "  ERROR: neutral OUT.OCC not found: $neutral_dir/OUT.OCC"
        return 1
    fi
    if [ ! -f "$charged_scf_dir/OUT.OCC" ]; then
        echo "  ERROR: charged OUT.OCC not found: $charged_scf_dir/OUT.OCC"
        return 1
    fi

    cp "$source_defect_input" "$run_defect_input"
    cp "$neutral_dir/OUT.OCC" "$run_dir/OUT.OCC0"
    cp "$charged_scf_dir/OUT.OCC" "$run_dir/OUT.OCC1"
}

submit_one() {
    local charged_scf_dir="$1"
    local run_dir="$2"
    local corr_dir="$3"
    local source_defect_input="$4"
    local marker="$run_dir/correction_${corr_dir}.submitted"
    local log="$run_dir/correction_${corr_dir}_submit.log"

    if [ "$FORCE" -ne 1 ] && [ -f "$marker" ] && ! is_finished "$run_dir" "$corr_dir"; then
        echo "  SKIP: already submitted marker exists"
        return 0
    fi

    prepare_occ_files "$charged_scf_dir" "$run_dir" "$source_defect_input" || return 1
    prepare_local_helpers "$run_dir"

    (
        cd "$run_dir" || exit 1
        echo "[$(date '+%F %T')] step 1: generate defect charge density"
        bash "$TOOL_DIR/1_get_rho.sh"
        echo "[$(date '+%F %T')] step 2: submit Coulomb integral jobs"
        bash "$TOOL_DIR/2_coulomb_integral.sh"
    ) > "$log" 2>&1
    local rc=$?
    if [ "$rc" -eq 0 ]; then
        date '+%F %T' > "$marker"
        echo "  SUBMITTED"
    else
        echo "  ERROR: submit failed, see $log"
    fi
    return "$rc"
}

collect_one() {
    local charged_scf_dir="$1"
    local run_dir="$2"
    local corr_dir="$3"
    local source_defect_input="$4"
    local marker="$run_dir/correction_${corr_dir}.done"
    local log="$run_dir/correction_${corr_dir}_collect.log"

    if ! is_finished "$run_dir" "$corr_dir"; then
        echo "  SKIP: potential jobs not finished"
        return 0
    fi
    if [ "$FORCE" -ne 1 ] && [ -f "$marker" ]; then
        echo "  SKIP: already collected marker exists"
        return 0
    fi

    cp "$source_defect_input" "$run_dir/defect.input"
    prepare_local_helpers "$run_dir"
    (
        cd "$run_dir" || exit 1
        echo "[$(date '+%F %T')] step 3: collect ICIC and PA"
        bash "$TOOL_DIR/3_get_results.sh"
    ) > "$log" 2>&1
    local rc=$?
    if [ "$rc" -eq 0 ]; then
        date '+%F %T' > "$marker"
        echo "  COLLECTED"
    else
        echo "  ERROR: collect failed, see $log"
    fi
    return "$rc"
}

total=0
submitted=0
collected=0
skipped=0
failed=0
PROCESSED_KEYS_TEXT=""

echo "Project root: $PROJECT_ROOT"
echo "Tool dir:     $TOOL_DIR"
echo "Mode:         $MODE"
if [ "${#FILTERS[@]}" -gt 0 ]; then
    echo "Filter:       ${FILTERS[*]}"
fi
echo

while IFS= read -r -d '' defect_input; do
    corr_dir="$(charge_folder_from_defect_input "$defect_input" || true)"
    if [ -z "$corr_dir" ]; then
        echo "[$defect_input] ERROR: cannot parse charged state from defect.input"
        failed=$((failed + 1))
        continue
    fi

    input_scf_dir="$(dirname "$defect_input")"
    q_dir="$(basename "$(dirname "$input_scf_dir")")"
    defect_root="$(dirname "$(dirname "$input_scf_dir")")"
    defect_dir="$(basename "$defect_root")"
    if ! matches_filter "$defect_dir"; then
        continue
    fi
    q_state="$(awk '/^charged[[:space:]]+state/{print $3; exit}' "$defect_input")"
    named_q_state="$(charge_from_input_name "$(basename "$defect_input")")"
    if [ -n "$named_q_state" ] && [ "$named_q_state" != "$q_state" ]; then
        echo "[$defect_input] ERROR: file name charge ($named_q_state) differs from content charge ($q_state)"
        failed=$((failed + 1))
        continue
    fi
    q_num="${q_state#+}"

    if [ "$q_dir" = "q_0" ]; then
        run_dir="$input_scf_dir"
        charged_scf_dir="$defect_root/q_${q_num}/scf"
        defect_input="$(ensure_stable_defect_input "$defect_input" "$q_state")"
    else
        charged_scf_dir="$input_scf_dir"
        run_dir="$(awk '/^neutral[[:space:]]/{print $2; exit}' "$defect_input")"
        if [ -z "$run_dir" ]; then
            run_dir="$defect_root/q_0/scf"
        fi
        defect_input="$(ensure_stable_defect_input "$defect_input" "$q_state")"
    fi

    if [ "$q_num" = "0" ]; then
        continue
    fi

    process_key="$defect_root|$q_state"
    if printf '%s' "$PROCESSED_KEYS_TEXT" | grep -Fqx "$process_key"; then
        continue
    fi
    PROCESSED_KEYS_TEXT="${PROCESSED_KEYS_TEXT}${process_key}
"

    if [ ! -d "$charged_scf_dir" ]; then
        echo "[$defect_dir/q_${q_num}] ERROR: charged scf directory not found: $charged_scf_dir"
        failed=$((failed + 1))
        continue
    fi
    if [ ! -d "$run_dir" ]; then
        echo "[$defect_dir/q_${q_num}] ERROR: neutral run directory not found: $run_dir"
        failed=$((failed + 1))
        continue
    fi

    total=$((total + 1))
    printf '[%s/q_%s] input=%s\n' "$defect_dir" "$q_num" "$defect_input"
    printf '  charged: %s\n' "$charged_scf_dir"
    printf '  run_dir: %s\n' "$run_dir"

    if is_finished "$run_dir" "$corr_dir"; then
        if [ "$MODE" = "submit" ]; then
            echo "  SKIP: finished ($corr_dir/REPORT and REPORT.0 exist)"
            skipped=$((skipped + 1))
            continue
        fi
    fi

    case "$MODE" in
        status)
            if is_finished "$run_dir" "$corr_dir"; then
                echo "  STATUS: finished"
            elif [ -f "$run_dir/correction_${corr_dir}.submitted" ]; then
                echo "  STATUS: submitted or running"
            else
                echo "  STATUS: not submitted"
            fi
            ;;
        submit)
            if submit_one "$charged_scf_dir" "$run_dir" "$corr_dir" "$defect_input"; then
                submitted=$((submitted + 1))
            else
                failed=$((failed + 1))
            fi
            ;;
        collect)
            before_done=0
            [ -f "$run_dir/correction_${corr_dir}.done" ] && before_done=1
            if collect_one "$charged_scf_dir" "$run_dir" "$corr_dir" "$defect_input"; then
                if [ "$before_done" -eq 0 ] && [ -f "$run_dir/correction_${corr_dir}.done" ]; then
                    collected=$((collected + 1))
                else
                    skipped=$((skipped + 1))
                fi
            else
                failed=$((failed + 1))
            fi
            ;;
    esac
done < <(
    {
        find "$PROJECT_ROOT/calculate" -path '*/q_0/scf/defect_*.input' -print0
        find "$PROJECT_ROOT/calculate" -path '*/q_*/scf/defect.input' ! -path '*/q_0/scf/defect.input' -print0
        find "$PROJECT_ROOT/calculate" -path '*/q_0/scf/defect.input' -print0
    } | sort -z
)

echo
echo "Summary:"
echo "  scanned:   $total"
echo "  submitted: $submitted"
echo "  collected: $collected"
echo "  skipped:   $skipped"
echo "  failed:    $failed"

if [ "$MODE" = "submit" ]; then
    echo
    echo "After Slurm jobs finish, run:"
    echo "  bash $SCRIPT_DIR/$(basename "$0") collect"
fi
