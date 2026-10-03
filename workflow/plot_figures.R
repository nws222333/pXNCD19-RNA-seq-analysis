try(suppressWarnings(Sys.setlocale('LC_ALL','en_US.UTF-8')), silent=TRUE)
suppressPackageStartupMessages({library(ggplot2);library(patchwork);library(dplyr);library(tidyr);library(svglite);library(ragg)})
args_all<-commandArgs(trailingOnly=FALSE)
script<-sub('^--file=','',args_all[grepl('^--file=',args_all)][1])
root<-normalizePath(file.path(dirname(script),'..'),winslash='/')
out<-file.path(root,'outputs','figures');dir.create(out,recursive=TRUE,showWarnings=FALSE)
read_tsv<-function(path)read.delim(file.path(root,path),check.names=FALSE,stringsAsFactors=FALSE)
source_table<-function(f)read_tsv(file.path('inputs/figure_sources',f))
theme_set(theme_classic(base_size=7,base_family='Arial')+theme(plot.title=element_text(face='bold',size=8),axis.text=element_text(size=6.5),axis.title=element_text(size=7),legend.text=element_text(size=6.5),legend.title=element_text(size=7),strip.text=element_text(size=7),plot.margin=margin(4,4,4,4)))
save_fig<-function(p,stem,w,h){
  svglite(file.path(out,paste0(stem,'.svg')),width=w/25.4,height=h/25.4);print(p);dev.off()
  cairo_pdf(file.path(out,paste0(stem,'.pdf')),width=w/25.4,height=h/25.4,family='Arial');print(p);dev.off()
  agg_png(file.path(out,paste0(stem,'.png')),width=w,height=h,units='mm',res=300);print(p);dev.off()
  agg_tiff(file.path(out,paste0(stem,'.tiff')),width=w,height=h,units='mm',res=600,compression='lzw');print(p);dev.off()
}
abund600<-read_tsv('outputs/abundance/C600_vs_CD19C_descriptive_abundance.tsv')
abundTH2<-read_tsv('outputs/abundance/TH2_vs_TH2C_descriptive_abundance.tsv')
res600<-read_tsv('outputs/DESeq2/C600_vs_CD19C_DESeq2_results.tsv')
repplot<-function(selection,abundance,groups,cols,title,ncol){
  selection$display<-gsub('_','-',selection$display,fixed=TRUE)
  x<-abundance%>%inner_join(selection%>%select(gene_id,display),by='gene_id')
  x$Group<-sub('-[123]$','',x$sample);x$Group<-factor(x$Group,levels=groups)
  x$display<-factor(x$display,levels=selection$display)
  ggplot(x,aes(Group,log2(FPKM+1),colour=Group))+
    geom_point(position=position_jitter(width=.06,height=0,seed=1),size=1.3)+
    stat_summary(fun=mean,geom='point',shape=18,size=2)+facet_wrap(~display,ncol=ncol)+
    scale_colour_manual(values=cols)+labs(title=title,x=NULL,y='log2(FPKM + 1)')+
    theme(legend.position='none',axis.text.x=element_text(size=6),strip.background=element_blank(),strip.text=element_text(face='italic',size=6.5))
}
cols600<-c(C600='#E64B35',CD19C='#00A087')
arn<-source_table('Fig3a_arn_source.tsv');env<-source_table('Fig3c_envelope_source.tsv');met<-source_table('Fig3d_metabolism_source.tsv')
pa<-repplot(arn,abund600,c('C600','CD19C'),cols600,'a  arn loci',3)
pc<-repplot(env,abund600,c('C600','CD19C'),cols600,'c  Envelope and signalling candidates',3)
bp<-read_tsv('outputs/ORA/C600_vs_CD19C_BP_all_terms.tsv');kg<-read_tsv('outputs/ORA/C600_vs_CD19C_KEGG_all_terms.tsv')
e<-bind_rows(kg[kg$ID=='eco02024',],bp[bp$ID=='GO:0009245',])
e$term<-factor(c('Quorum sensing (KEGG)','Lipid A biosynthesis (GO)'),levels=c('Lipid A biosynthesis (GO)','Quorum sensing (KEGG)'))
e$detail<-sprintf('k = %d; P = %.3g; FDR = %.3g',e$gene_number,e$pvalue,e$qvalue)
pb<-ggplot(e,aes(enrich_factor,term))+
  geom_segment(aes(x=0,xend=enrich_factor,yend=term),colour='grey65',linewidth=.4)+
  geom_point(aes(size=gene_number,colour=-log10(qvalue)),shape=18)+
  geom_text(aes(label=detail),hjust=0,nudge_x=.35,size=2.3)+
  scale_size_continuous(name='DEG count',range=c(2,4),breaks=c(5,11))+
  scale_colour_gradient(low='#4292C6',high='#B2182B',name='-log10(FDR)',limits=c(0,1),breaks=c(0,.5,1))+
  scale_x_continuous(limits=c(0,9),expand=expansion(mult=c(.01,.01)))+
  labs(title='b  Candidate enrichment terms',subtitle='Both terms have FDR >= 0.05',x='Enrichment factor',y=NULL)+
  guides(colour=guide_colourbar(barwidth=unit(25,'mm'),barheight=unit(2,'mm'),title.position='top'))+
  theme(plot.subtitle=element_text(size=6.5),axis.text.y=element_text(size=6),legend.position='bottom',legend.box='vertical',legend.key.size=unit(3,'mm'))
met<-met%>%select(display,gene_id)%>%left_join(res600,by='gene_id')
met$display<-factor(gsub('_','-',met$display,fixed=TRUE),levels=gsub('_','-',met$display,fixed=TRUE))
pd<-ggplot(met,aes(log2FoldChange,display))+
  geom_segment(aes(xend=0,yend=display),colour='grey65',linewidth=.4)+
  geom_point(aes(size=-log10(padj)),colour='#3C5488')+
  geom_vline(xintercept=0,colour='grey60',linewidth=.3)+
  scale_size_continuous(name='-log10(DEG FDR)',range=c(1.5,3.5),breaks=c(3,5,7))+
  labs(title='d  Metabolism-associated loci',x='RNA log2 fold change (CD19C / C600)',y=NULL)+
  theme(axis.text.y=element_text(face='italic'),legend.position='bottom')
fig3<-(pa/pb+plot_layout(heights=c(1.5,1)))|(pc/pd+plot_layout(heights=c(1.6,1.3)))
save_fig(fig3,'Fig3_C600_corrected',183,178)
write.table(e,file.path(out,'Fig3b_enrichment_source.tsv'),sep='\t',quote=FALSE,row.names=FALSE)

cols<-c(TH2='#4DBBD5',TH2C='#3C5488')
module<-source_table('Fig4_module_DNA_RNA.tsv');module$display[module$gene_id=='TH2_HKNPDCAD_05211']<-'dfrA12'
module$display<-gsub('_','-',module$display,fixed=TRUE)
samples<-c('TH2-1','TH2-2','TH2-3','TH2C-1','TH2C-2','TH2C-3')
wide<-abundTH2%>%select(gene_id,sample,FPKM)%>%pivot_wider(names_from=sample,values_from=FPKM)
fpkm<-as.matrix(wide[match(module$gene_id,wide$gene_id),samples]);mode(fpkm)<-'numeric'
z<-t(scale(t(fpkm)));z[!is.finite(z)]<-0
ml<-as.data.frame(z);ml$display<-module$display
ml<-pivot_longer(ml,all_of(samples),names_to='Sample',values_to='z')
ml$display<-factor(ml$display,levels=rev(module$display));ml$Sample<-factor(ml$Sample,levels=samples)
module$display<-factor(module$display,levels=rev(module$display))
state_labels<-c(reference_copy_not_detected='Undetected',ambiguous_repetitive='Repeat',retained_DNA='Retained',partial_or_low_DNA_support='Low support')
state_cols<-c(Undetected='#B24745',Repeat='#D8BD79',Retained='#528B82','Low support'='#989BA0')
module$DNA<-state_labels[module$DNA_state]
prna<-ggplot(ml,aes(Sample,display,fill=z))+geom_tile(colour='white',linewidth=.1)+
  scale_fill_gradient2(low='#2166AC',mid='white',high='#B2182B',midpoint=0,breaks=c(-1,0,1),name='RNA row z score')+
  guides(fill=guide_colourbar(barwidth=unit(27,'mm'),barheight=unit(2.5,'mm')))+
  labs(x=NULL,y=NULL,title='c  RNA abundance and DNA support')+
  theme(axis.line=element_blank(),axis.ticks=element_blank(),axis.text.x=element_text(angle=45,hjust=1),axis.text.y=element_text(size=6.7,face='italic'),legend.position='bottom')
pdna<-ggplot(module,aes('DNA',display,fill=DNA))+geom_tile(colour='white',linewidth=.1)+
  scale_fill_manual(values=state_cols,name='DNA support',drop=FALSE)+labs(x=NULL,y=NULL)+
  theme(axis.line=element_blank(),axis.ticks=element_blank(),axis.text.y=element_blank(),legend.position='bottom',plot.margin=margin(17,2,4,0),axis.text.x=element_text(angle=45,hjust=1))+
  guides(fill=guide_legend(nrow=1,byrow=TRUE))
heat<-prna+pdna+plot_layout(widths=c(6,1))
bp<-read_tsv('outputs/ORA/TH2_vs_TH2C_BP_all_terms.tsv');bp<-bp[bp$qvalue<.05,]
bp$Description<-factor(bp$Description,levels=rev(bp$Description))
pa<-ggplot(bp,aes(enrich_factor,Description,size=gene_number,colour=-log10(qvalue)))+
  geom_point()+scale_colour_gradient(low='#4292C6',high='#B2182B',name='-log10(FDR)',labels=function(x)sprintf('%.2f',x))+
  scale_size_continuous(range=c(1.5,4),breaks=sort(unique(bp$gene_number)),name='DEG count')+
  labs(title='a  GO biological process',x='Enrichment factor',y=NULL)+theme(legend.position='bottom',axis.text.y=element_text(size=6.3))
stress<-source_table('Fig4d_stress_source.tsv')
pb<-repplot(stress,abundTH2,c('TH2','TH2C'),cols,'b  Surface and stress loci',3)
fig4<-(pa/pb+plot_layout(heights=c(.8,1.65)))|heat
fig4<-fig4+plot_layout(widths=c(1.28,1),guides='collect')&theme(legend.position='bottom',legend.box='vertical',legend.spacing.y=unit(0,'mm'),legend.key.size=unit(2,'mm'))
save_fig(fig4,'Fig4_TH2_DNA_corrected',183,165)
write.table(bp,file.path(out,'Fig4a_enrichment_source.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
writeLines(capture.output(sessionInfo()),file.path(root,'expected_results/plot_sessionInfo.txt'))
cat('Saved Figures 3 and 4 in PDF, SVG, PNG and TIFF formats\n')
