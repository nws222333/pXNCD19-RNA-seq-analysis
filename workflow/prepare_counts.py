"""Extract STAR reverse-strand GeneCounts and the current unique-read denominators."""
from pathlib import Path
import csv
ROOT=Path(__file__).resolve().parents[1]
def read(path):
    with path.open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f,delimiter='\t'))
def write(path,fields,rows):
    with path.open('w',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter='\t');w.writeheader();w.writerows(rows)
manifest=read(ROOT/'inputs/sample_manifest.tsv');cache={}
for sample in manifest:
    folder=ROOT/'outputs/mapping'/sample['sample'];counts={}
    with (folder/'ReadsPerGene.out.tab').open(encoding='utf-8') as f:
        for line in f:
            p=line.rstrip().split('\t')
            if not p[0].startswith('N_'):counts[p[0]]=int(p[3])
    cache[sample['sample']]=counts
    for line in (folder/'Log.final.out').read_text(encoding='utf-8').splitlines():
        p=line.split('|')
        if len(p)==2 and p[0].strip()=='Uniquely mapped reads number':sample['uniquely_mapped_read_pairs']=int(p[1].strip())
for host,cmp in [('C600','C600_vs_CD19C'),('TH2','TH2_vs_TH2C')]:
    names=[s['sample'] for s in manifest if s['host']==host]
    metadata=read(ROOT/'inputs/references'/f'{host}_plus_pXNCD19_gene_metadata.tsv')
    for kind in ['composite','host']:
        chosen=metadata if kind=='composite' else [r for r in metadata if r['reference_component']==host]
        rows=[{'gene_id':r['gene_id'],**{s:cache[s].get(r['gene_id'],0) for s in names}}for r in chosen]
        write(ROOT/'inputs/counts'/f'{cmp}_{kind}_counts.tsv',['gene_id',*names],rows)
write(ROOT/'inputs/sample_manifest.tsv',list(manifest[0]),manifest)
