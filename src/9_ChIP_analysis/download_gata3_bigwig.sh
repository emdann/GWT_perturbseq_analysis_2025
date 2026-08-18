#!/usr/bin/env bash
# Download GATA3 ChIP-seq bigwig from ChIP Atlas (SRX092314, Th1 Cells, hg38)
# ChIP Atlas: https://chip-atlas.dbcls.jp
#
# Usage: ./download_gata3_bigwig.sh [outdir]

set -euo pipefail

OUTDIR="${1:-/mnt/oak/users/emma/data/ChIP_atlas_data/GATA3}"
GENOME="hg38"
SRX_ID="SRX092314"
TF="GATA3"
BASE_URL="https://chip-atlas.dbcls.jp/data/${GENOME}/eachData/bw"

mkdir -p "${OUTDIR}"

OUTFILE="${OUTDIR}/${SRX_ID}.bw"
URL="${BASE_URL}/${SRX_ID}.bw"

if [[ -f "${OUTFILE}" ]]; then
    echo "[SKIP] ${SRX_ID} (${TF}) bigwig — already exists at ${OUTFILE}"
    exit 0
fi

echo "[DL] ${SRX_ID} (${TF}) → ${OUTFILE}"
echo "URL: ${URL}"
wget -q --show-progress -O "${OUTFILE}" "${URL}" || {
    echo "[ERROR] Failed to download ${SRX_ID} bigwig; removing partial file"
    rm -f "${OUTFILE}"
    exit 1
}

echo ""
echo "Done. GATA3 bigwig saved to: ${OUTFILE}"
