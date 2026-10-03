"""FPKM and CPM are descriptive. Neither is an input to DESeq2."""
from pathlib import Path
import csv
ROOT=Path(__file__).resolve().parents[1]
def read(path):
    with path.open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f,delimiter='\t'))
manifest=read(ROOT/'inputs/sample_manifest.tsv');out=ROOT/'outputs/abundance';out.mkdir(parents=True,exist_ok=True)
for host,cmp in [('C600','C600_vs_CD19C'),('TH2','TH2_vs_TH2C')]:
    metadata={r['gene_id']:r for r in read(ROOT/'inputs/references'/f'{host}_plus_pXNCD19_gene_metadata.tsv')}
    samples=[s for s in manifest if s['host']==host]
    rows=[]
    for g in read(ROOT/'inputs/counts'/f'{cmp}_composite_counts.tsv'):
        m=metadata[g['gene_id']];length=int(m['length_bp'])
        for s in samples:
            c=int(g[s['sample']]);unique=int(s['uniquely_mapped_read_pairs'])
            rows.append({'gene_id':g['gene_id'],'sample':s['sample'],'reference_component':m['reference_component'],
                         'length_bp':length,'count':c,'uniquely_mapped_read_pairs':unique,
                         'CPM':c*1e6/unique,'FPKM':c*1e9/(length*unique)})
    with (out/f'{cmp}_descriptive_abundance.tsv').open('w',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
