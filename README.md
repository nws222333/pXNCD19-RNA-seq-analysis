# pXNCD19 RNA-seq analysis

This release reproduces the recipient-reference RNA-seq contrasts for C600/CD19C and TH2/TH2C. It includes frozen composite references, integer counts, annotation memberships, sample labels, and the analysis and Figure 3/4 plotting scripts. Release 2026-10-03 corrects the RNA enrichment BH families. Both contrasts retain their original DEG sets (312 and 181 loci).

## Software

The completed analysis used STAR 2.7.11b, R 4.6.1, DESeq2 1.52.0 and Python 3.12.14. Python scripts use the standard library. `expected_results/R_sessionInfo.txt` records the R/Bioconductor dependencies. Figure rendering additionally uses ggplot2, dplyr, tidyr, patchwork, svglite and ragg; `expected_results/plot_sessionInfo.txt` records their installed versions.

STAR mapping requires Linux or WSL with STAR and zcat on PATH. DESeq2, abundance calculation, enrichment, and plotting also run on Windows. Files use UTF-8. Run commands from the repository root.

## Reproduce from the supplied counts

```sh
Rscript workflow/run_deseq2.R
python workflow/calculate_abundance.py
python workflow/run_enrichment.py
Rscript workflow/plot_figures.R
```

The outputs go to `outputs/`. `expected_results/` supplies the completed statistical tables and software sessions. `SHA256SUMS.tsv` identifies the release files.

## Optional rerun from FASTQ

RNA-seq data are deposited under BioProject PRJNA1524675. This repository excludes FASTQ. The original mapping used the provider's preliminary quality-filtered `*_good_*` paired-end files, not the separately reported rRNA-filtered files. The historical mapping input FASTQ SHA256 manifest is in `inputs/provenance/`. The sample manifest gives the exact original filenames. Place or symlink the twelve pairs into one directory, then run

```sh
bash workflow/run_mapping.sh /path/to/provider_good_FASTQ_directory
python workflow/prepare_counts.py
Rscript workflow/run_deseq2.R
python workflow/calculate_abundance.py
python workflow/run_enrichment.py
Rscript workflow/plot_figures.R
```

`prepare_counts.py` regenerates the count files and mapped-pair denominators in `inputs/`. Save the frozen release before regenerating those files. SRA downloads may differ from the provider-filtered input. Exact reconstruction of the provider's preprocessing requires its unavailable filtering script, SOAP version, and NCBI/SILVA database releases. Do not substitute guessed preprocessing parameters or describe raw-SRA remapping as a byte-identical historical rerun.

## References and counting

C600 uses NCBI assembly GCF_003367885.1 (ASM336788v1); TH2 uses the supplied Prokka-annotated XNTH2 assembly. Each includes donor assembly contig 4, the 95,080-bp pXNCD19 sequence. The frozen FASTA/GTF and metadata are the exact mapping references. Component-prefixed gene IDs distinguish recipient and incoming-plasmid loci. GTF exon features derive from GFF CDS intervals; the metadata length for a gene is the union of its CDS intervals. Local TH2 and donor assembly public accessions and original Prokka version were unavailable.

STAR builds indices with `--sjdbOverhang 149 --genomeSAindexNbases 10`. Mapping uses `--quantMode GeneCounts --outFilterMultimapNmax 1 --outFilterMismatchNoverLmax 0.04 --alignIntronMax 1 --outSAMtype None --genomeLoad NoSharedMemory`. Threads default to six. Counts come from column four (reverse-strand counts) of `ReadsPerGene.out.tab`, with its first four summary rows excluded. No separate featureCounts step occurs.

## Differential expression and descriptive abundance

Each six-library contrast fits `~ group`, with recipient as the reference level. DESeq2 excludes incoming pXNCD19 loci and all-zero recipient-reference loci. The script specifies the historical DESeq defaults explicitly: median-of-ratios size factors (`sfType="ratio"`), parametric dispersion fit, Wald tests, `betaPrior=FALSE`, and an unshrunk transconjugant/recipient log2 fold change. `results()` uses `alpha=0.01`, independent filtering, `pAdjustMethod="BH"`, and `lfcThreshold=0`. The DEG definition subsequently requires padj <0.01 and absolute log2 fold change >=1. The Wald null tests zero log2 fold change, not a fold-change threshold of one. The model contains no batch or host-by-plasmid interaction term.

CPM = count * 10^6 / uniquely mapped paired reads. FPKM = count * 10^9 / (CDS-union length in bp * uniquely mapped paired reads). The denominator includes both components of the matching composite reference. FPKM and CPM provide descriptive abundance and plots; DESeq2 takes integer counts and estimates its own size factors. Neither calculation normalizes RNA to plasmid DNA copy number.

## Annotation and over-representation

`inputs/annotations/*_exact_CDS_crosswalk.tsv` documents the transfer of Biomarker Technologies annotations across exact CDS sequence matches. Ambiguous matches are excluded. `gene_sets.tsv` freezes the actual memberships used, so inference does not depend on a live database. GO omits obsolete labels and the BP/MF/CC root terms. It uses supplied annotation memberships without additional GO ancestor propagation. The provider's GO release is unknown.

KEGG memberships originated from the E. coli `eco` gene list and gene-to-pathway links. The historical mapping matches the entire symbol string before the semicolon in the provider's KEGG annotation against the corresponding entire symbol field in the eco gene list. It does not expand synonym lists, perform KO-based annotation or reconstruct all pathways. This restricted matching limits coverage. Only the derived gene memberships are distributed here, not the downloaded KEGG source database files. The provenance manifest records their retrieval timestamps as filesystem evidence and hashes, not an official database release.

For each contrast and each family (GO BP, MF, CC and KEGG), the background consists of recipient-reference genes with finite DESeq2 padj and valid annotations in that family. The foreground is the union of higher- and lower-abundance DEGs intersected with this background. The test is the one-sided hypergeometric upper tail, computed in Python's standard library. BH includes every term represented in that background, including terms overlapping zero or one DEG. BH runs separately within each of the eight contrast/family combinations. FDR <0.05 defines enrichment. P <0.05 without this FDR threshold is nominal only. Each output contains N, n, K, k, P, q, overlap gene IDs and the BH family size. Background/foreground gene and term lists accompany the outputs.

The previous implementation filtered terms with fewer than two DEG overlaps before BH. Release 2026-10-03 removes this filter. C600 quorum sensing now has q=0.1049346773 and lipid A biosynthetic process has q=1. TH2 has two significant BP terms, two MF terms and one CC term. Neither contrast has a significant KEGG term. These RNA enrichments do not establish pathway activity or per-DNA-copy regulation.

## Figure sources and scope

Figure 3/4 scripts derive RNA abundance and fold changes from the supplied counts and current DESeq2 outputs, and read enrichment from the corrected all-term ORA tables. Candidate-locus selections are frozen in `inputs/figure_sources/`. Figure 4c additionally uses the completed descriptive DNA support table from the TH2/TH2C analysis. This repository covers the RNA analysis and those affected figures.
