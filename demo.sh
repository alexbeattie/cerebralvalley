#!/usr/bin/env bash
#
# demo.sh — guided tour of the causality review (step 2).
#
# Each case asks the clinician's question: does a reported variant actually
# explain THIS patient? It hits the public HPO API, so expect a couple seconds
# per case.
#
# Usage:
#   ./demo.sh            # pause between cases
#   ./demo.sh --no-pause # run straight through

set -u

if [[ -x ".venv/bin/python" ]]; then
  PY=".venv/bin/python"
elif [[ -x "../biohack/.venv/bin/python" ]]; then
  PY="../biohack/.venv/bin/python"   # shared venv from the main worktree
else
  PY="python3"
fi

PAUSE=1
[[ "${1:-}" == "--no-pause" ]] && PAUSE=0

if [[ -t 1 ]]; then
  BOLD=$'\033[1m'; DIM=$'\033[2m'; CYAN=$'\033[36m'; GREEN=$'\033[32m'; RESET=$'\033[0m'
else
  BOLD=""; DIM=""; CYAN=""; GREEN=""; RESET=""
fi

step=0
case_run() {
  local title="$1"; shift
  local point="$1"; shift
  step=$((step + 1))
  echo
  echo "${CYAN}########################################################################${RESET}"
  echo "${CYAN}#${RESET} ${BOLD}DEMO ${step}: ${title}${RESET}"
  echo "${CYAN}#${RESET} ${DIM}Why it matters:${RESET} ${point}"
  echo "${CYAN}########################################################################${RESET}"
  "$PY" -m causality_review.cli "$@"
  if [[ "$PAUSE" -eq 1 ]]; then
    echo; read -r -p "${DIM}[enter] for the next case...${RESET}" _
  fi
}

echo "${GREEN}${BOLD}"
echo "  Causality Review — step 2: does a reported variant explain THIS patient?"
echo "${RESET}${DIM}  The lab hands over 2-5 suspicious variants. The clinician has the"
echo "  patient's full phenotype. This runs the match in reverse and ranks the"
echo "  reported variants by how well each explains this specific patient —"
echo "  separating explained from unexplained features, and flagging when the"
echo "  leftovers point to a second cause or a genome re-analysis.${RESET}"

# 1. Clean best-fit: the epilepsy gene explains the whole epilepsy picture -----
case_run \
  "Clean best fit — the right gene explains everything" \
  "An epileptic-encephalopathy picture against SCN1A (Dravet) vs MYH7 (a cardiomyopathy gene). SCN1A should explain every feature and rank far above the cardiac gene. This is the 'obvious' case done objectively." \
  --variant "SCN1A:c.3637C>T" --variant "MYH7:c.1063G>A" \
  --hpo "Seizure" --hpo "Global developmental delay" --hpo "Febrile seizure" --hpo "Ataxia"

# 2. Partial + orphan feature: the second-cause / re-analysis trigger ----------
case_run \
  "Partial fit — one feature no reported variant explains" \
  "Same case, but the patient also has a cleft palate. SCN1A still ranks first, but now explains only 4/5 — and the tool flags the orphan feature as possible second cause or genome re-analysis, instead of letting the clinician wave it away." \
  --variant "SCN1A:c.3637C>T" --variant "MYH7:c.1063G>A" \
  --hpo "Seizure" --hpo "Global developmental delay" --hpo "Febrile seizure" --hpo "Ataxia" --hpo "Cleft palate"

# 3. Two partial explanations: the overfitting trap made visible --------------
case_run \
  "Two partial matches — resisting the overfit" \
  "A cardiac patient with two plausible reported variants (KCNQ1 long-QT vs MYH7 cardiomyopathy). Neither explains the whole picture; the tool shows each as a partial, explicit fit rather than letting intuition force one to 'fit'." \
  --variant "KCNQ1:c.477+1G>A" --variant "MYH7:c.1063G>A" \
  --hpo "Hypertrophic cardiomyopathy" --hpo "Prolonged QT interval" --hpo "Syncope"

echo
echo "${GREEN}${BOLD}  Demo complete.${RESET}"
echo "${DIM}  Every gene-phenotype link carries an openable HPO source URL. Ranking is"
echo "  by fit to THIS patient, not by the lab's ordering. Unexplained features are"
echo "  surfaced, not buried — re-analysis is a feature, not a failure.${RESET}"
echo
