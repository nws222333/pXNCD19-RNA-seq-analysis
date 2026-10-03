"""Exact-input RNA ORA. Python standard library only; BH over all represented terms."""
from pathlib import Path
from collections import defaultdict
import csv,math,json,sys

ROOT=Path(__file__).resolve().parents[1]
def read(path):
    with path.open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f,delimiter='\t'))
def write(path,rows):
    path.parent.mkdir(parents=True,exist_ok=True)
    if not rows:raise ValueError('Empty analysis table '+str(path))
    with path.open('w',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
def bh(values):
    order=sorted(range(len(values)),key=values.__getitem__);out=[1.]*len(values);last=1.
    for j in range(len(order)-1,-1,-1):
        i=order[j];last=min(last,values[i]*len(values)/(j+1));out[i]=last
    return out
def hypergeom_upper(k,N,K,n):
    if k<=max(0,n+K-N):return 1.
    hi=min(K,n)
    if k>hi:return 0.
    logchoose=lambda a,b:math.lgamma(a+1)-math.lgamma(b+1)-math.lgamma(a-b+1)
    logp=logchoose(K,k)+logchoose(N-K,n-k)-logchoose(N,n);logs=[logp]
    for i in range(k,hi):
        logp+=math.log(K-i)-math.log(i+1)+math.log(n-i)-math.log(N-K-n+i+1);logs.append(logp)
    peak=max(logs);return min(1.,math.exp(peak)*sum(math.exp(v-peak) for v in logs))
def main():
    annotations=read(ROOT/'inputs/annotations/gene_sets.tsv');summaries=[]
    for cmp in ['C600_vs_CD19C','TH2_vs_TH2C']:
        result=read(ROOT/'outputs/DESeq2'/f'{cmp}_DESeq2_results.tsv');tested=set();up=set();down=set()
        for r in result:
            try:q=float(r['padj']);lfc=float(r['log2FoldChange'])
            except (ValueError,TypeError):continue
            if not math.isfinite(q):continue
            tested.add(r['gene_id'])
            if q<.01 and math.isfinite(lfc) and abs(lfc)>=1:(up if lfc>0 else down).add(r['gene_id'])
        for family in ['BP','MF','CC','KEGG']:
            terms=defaultdict(set);descriptions={};background=set()
            for r in annotations:
                if r['comparison']!=cmp or r['family']!=family or r['gene_id'] not in tested:continue
                terms[r['term_id']].add(r['gene_id']);background.add(r['gene_id']);descriptions[r['term_id']]=r['description']
            foreground=(up|down)&background;N=len(background);n=len(foreground);rows=[]
            for term,members in sorted(terms.items()):
                overlap=foreground&members;k=len(overlap);K=len(members)
                p=hypergeom_upper(k,N,K,n)
                direction='up&down' if overlap&up and overlap&down else 'up' if overlap&up else 'down' if overlap&down else ''
                ef=(k/n)/(K/N) if n and N and K else 0.
                rows.append({'ID':term,'Description':descriptions[term],'GeneRatio':f'{k}/{n}','BgRatio':f'{K}/{N}',
                             'enrich_factor':ef,'pvalue':p,'qvalue':0.,'geneID':';'.join(sorted(overlap)),
                             'gene_number':k,'Contained':direction,'Ontology':{'BP':'Biological Process','MF':'Molecular Function','CC':'Cellular Component','KEGG':'KEGG'}[family],
                             'GeneRatio_num':k/n if n else 0.,'BgRatio_num':K/N if N else 0.,
                             'negLog10p':-math.log10(max(p,1e-300)),'negLog10q':0.,'direction':direction,
                             'background_genes':N,'foreground_genes':n,'background_term_genes':K,'BH_family_terms':len(terms)})
            adjusted=bh([r['pvalue']for r in rows])
            for r,q in zip(rows,adjusted):r['qvalue']=q;r['negLog10q']=-math.log10(max(q,1e-300))
            rows.sort(key=lambda r:(r['qvalue'],r['pvalue'],r['ID']))
            write(ROOT/'outputs/ORA'/f'{cmp}_{family}_all_terms.tsv',rows)
            write(ROOT/'outputs/ORA'/f'{cmp}_{family}_background_genes.tsv',
                  [{'gene_id':g,'is_DEG_foreground':int(g in foreground),'direction':'up' if g in up else 'down' if g in down else ''} for g in sorted(background)])
            # All k=0 and k=1 terms stay in the BH family; no overlap-dependent prefilter.
            summaries.append({'comparison':cmp,'family':family,'background_genes':N,'foreground_genes':n,
                              'BH_family_terms':len(rows),'terms_FDR_lt_0_05':sum(r['qvalue']<.05 for r in rows)})
    (ROOT/'outputs/ORA/family_summary.json').write_text(json.dumps(summaries,indent=2),encoding='utf-8')
    (ROOT/'outputs/ORA/python_version.txt').write_text(sys.version+'\n',encoding='utf-8')
    print(json.dumps(summaries))
if __name__=='__main__':main()
