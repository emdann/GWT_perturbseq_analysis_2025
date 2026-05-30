#!/usr/bin/env bash
# Download ChIP-seq peak calls (BED) from ChIP Atlas using experiment metadata table
# Genome: hg38 | ChIP Atlas: https://chip-atlas.dbcls.jp
#
# Usage: ./download_chipatlas_tracks.sh [outdir] [manifest_tsv] [threshold]
# threshold: q-value cutoff as 05/10/20/50 (q < 1e-N); default 20

set -euo pipefail

OUTDIR="${1:-/mnt/oak/users/emma/data/ChIP_atlas_data}"
MANIFEST="${2:-$(dirname "$0")/chipatlas_cd4_tf_experiments.tsv}"
GENOME="hg38"
THRESHOLD="${3:-20}"  # q-value threshold: 05, 10, 20, or 50 (q < 1e-N)
BASE_URL="https://chip-atlas.dbcls.jp/data/${GENOME}/eachData/bed${THRESHOLD}"

mkdir -p "${OUTDIR}"

if [[ ! -f "${MANIFEST}" ]]; then
  echo "[ERROR] Manifest not found: ${MANIFEST}"
  exit 1
fi

# Count data rows (skip header)
total=$(tail -n +2 "${MANIFEST}" | wc -l)
echo "Downloading ${total} peak call BED files (q < 1e-${THRESHOLD}) to ${OUTDIR}/"
echo "Manifest: ${MANIFEST}"
echo ""

# ---------------------------------------------------------------------------
# Download: read srx_id (col 1) and tf (col 2) from TSV, skip header
# ---------------------------------------------------------------------------
while IFS=$'\t' read -r srx tf _rest; do
  tfdir="${OUTDIR}/${tf}"
  mkdir -p "${tfdir}"
  outfile="${tfdir}/${srx}.${THRESHOLD}.bed"
  url="${BASE_URL}/${srx}.${THRESHOLD}.bed"

  if [[ -f "${outfile}" ]]; then
    echo "[SKIP] ${srx} (${tf}) — already exists"
    continue
  fi

  echo "[DL] ${srx} (${tf}) → ${outfile}"
  wget -q --show-progress -O "${outfile}" "${url}" || {
    echo "[ERROR] Failed to download ${srx}; removing partial file"
    rm -f "${outfile}"
  }
done < <(tail -n +2 "${MANIFEST}")

echo ""
echo "Done. Peak BED files organized in ${OUTDIR}/ by TF."
echo ""
echo "Summary of downloaded files:"
# List unique TFs from manifest and count downloaded tracks per TF
tail -n +2 "${MANIFEST}" | awk -F'\t' '{print $2}' | sort -u | while read -r tf; do
  count=$(ls "${OUTDIR}/${tf}/"*.bed 2>/dev/null | grep -c "\.${THRESHOLD}\.bed" || true)
  echo "  ${tf}: ${count} tracks"
done
