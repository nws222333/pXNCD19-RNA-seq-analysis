#!/usr/bin/env bash
set -euo pipefail
ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
DATA=${1:?Usage: bash workflow/run_mapping.sh /path/to/provider_good_FASTQ_directory}
THREADS=${THREADS:-6}
mkdir -p "$ROOT/outputs/mapping"
for host in C600 TH2; do
  index="$ROOT/outputs/mapping/index_${host}"
  mkdir -p "$index"
  STAR --runThreadN "$THREADS" --runMode genomeGenerate --genomeDir "$index" \
    --genomeFastaFiles "$ROOT/inputs/references/${host}_plus_pXNCD19.fa" \
    --sjdbGTFfile "$ROOT/inputs/references/${host}_plus_pXNCD19.gtf" \
    --sjdbOverhang 149 --genomeSAindexNbases 10
done
tail -n +2 "$ROOT/inputs/sample_manifest.tsv" | while IFS=$'\t' read -r sample host comparison group code read1 read2 unique; do
  dest="$ROOT/outputs/mapping/$sample";mkdir -p "$dest"
  STAR --runThreadN "$THREADS" --genomeDir "$ROOT/outputs/mapping/index_${host}" \
    --readFilesIn "$DATA/$read1" "$DATA/$read2" --readFilesCommand zcat \
    --outFileNamePrefix "$dest/" --outTmpDir "/tmp/pXNCD19_STAR_${sample}_$$" \
    --outSAMtype None --quantMode GeneCounts --outFilterMultimapNmax 1 \
    --outFilterMismatchNoverLmax 0.04 --alignIntronMax 1 --genomeLoad NoSharedMemory
done
