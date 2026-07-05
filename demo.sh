#!/usr/bin/env bash
#
# demo.sh — a guided tour of the variant-curation assistant.
#
# Runs five real variants end to end (VEP -> gnomAD -> ClinVar), each chosen to
# make one point. Before every run it prints WHY the case matters, so the demo
# tells a story instead of just dumping output. Each variant hits three live
# public APIs, so expect a few seconds per case.
#
# Usage:
#   ./demo.sh            # run the whole reel, pausing between cases
#   ./demo.sh --no-pause # run straight through (good for recording)

set -u

# --- interpreter: prefer the project venv, fall back to python3 -------------
if [[ -x ".venv/bin/python" ]]; then
  PY=".venv/bin/python"
else
  PY="python3"
fi

PAUSE=1
[[ "${1:-}" == "--no-pause" ]] && PAUSE=0

# --- pretty-printing helpers (fall back to plain text if not a TTY) ---------
if [[ -t 1 ]]; then
  BOLD=$'\033[1m'; DIM=$'\033[2m'; CYAN=$'\033[36m'; GREEN=$'\033[32m'; RESET=$'\033[0m'
else
  BOLD=""; DIM=""; CYAN=""; GREEN=""; RESET=""
fi

step=0

# case "<title>" "<what this proves>" <cli args...>
case_run() {
  local title="$1"; shift
  local point="$1"; shift
  step=$((step + 1))

  echo
  echo "${CYAN}########################################################################${RESET}"
  echo "${CYAN}#${RESET} ${BOLD}DEMO ${step}: ${title}${RESET}"
  echo "${CYAN}#${RESET} ${DIM}Why it matters:${RESET} ${point}"
  echo "${CYAN}#${RESET} ${DIM}\$ ${PY} -m variant_curator.cli $*${RESET}"
  echo "${CYAN}########################################################################${RESET}"

  "$PY" -m variant_curator.cli "$@"

  if [[ "$PAUSE" -eq 1 ]]; then
    echo
    read -r -p "${DIM}[enter] for the next case...${RESET}" _
  fi
}

echo "${GREEN}${BOLD}"
echo "  Variant-Curation Assistant — live demo"
echo "${RESET}${DIM}  One gene + one HGVS change in, a fully sourced evidence bundle out."
echo "  The hour a curator spends opening VEP, gnomAD, and ClinVar by hand,"
echo "  done in seconds — with every line linked back to its source, and a"
echo "  hard refusal to assert anything the data does not support.${RESET}"

# 1. Flagship: pathogenic missense, ancestry gap visible ---------------------
case_run \
  "Pathogenic missense — and the ancestry gap it exposes" \
  "28 of 30 observed alleles are Non-Finnish European; ZERO in African, East/South Asian, Latino, Ashkenazi, Middle Eastern. For a non-European patient, 'rare = suspicious' has almost no population behind it — the tool shows that instead of hiding it in a global number." \
  --gene BRCA1 --hgvs "c.181T>G"

# 2. Contrast: common benign polymorphism ------------------------------------
case_run \
  "Common benign polymorphism — the mirror image" \
  "Global allele frequency ~23% (roughly 1 in 4 people carry it), so it cannot be causing disease; ClinVar agrees (Benign). Note it honestly reports the ClinVar conflict flag and the lower East-Asian frequency rather than pretending the data is tidy." \
  --gene SCN5A --hgvs "c.1673A>G"

# 3. Strongest evidence class: a truncating (nonsense) variant ---------------
case_run \
  "Truncating (nonsense) variant — the highest-confidence break" \
  "A stop-gain chops the protein short: in BRCA1 that is the strongest 'this is broken' signal a clinician has. Ultra-rare (6 alleles in ~1.6M), Pathogenic in ClinVar. Proves both supported variant types work, not just missense." \
  --gene BRCA1 --hgvs "c.5251C>T"

# 4. Trust: refuse an unsupported gene ---------------------------------------
case_run \
  "Out-of-scope GENE — refuses instead of guessing" \
  "EGFR is not on the reviewed allowlist, so the tool stops with a clear message. An assistant that says 'I won't touch this' is trustworthy in a way one that confidently guesses on everything is not." \
  --gene EGFR --hgvs "c.2369C>T"

# 5. Trust: refuse an out-of-scope consequence -------------------------------
case_run \
  "Out-of-scope CONSEQUENCE — abstains after resolving" \
  "It resolves the variant, sees it is synonymous (outside the missense/nonsense MVP scope), and abstains. This is the 'abstain where the evidence isn't there' premise, already working." \
  --gene BRCA1 --hgvs "c.2082C>T"

echo
echo "${GREEN}${BOLD}  Demo complete.${RESET}"
echo "${DIM}  Reminder: today the tool GATHERS and SOURCES the evidence; the final"
echo "  ACMG pathogenic/benign call is the next layer. The win shown here is the"
echo "  hour of manual lookup collapsed to seconds, fully cited, with no bluffing.${RESET}"
echo
