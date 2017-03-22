rm(list=ls())
options(java.parameters = "-Xmx8000m")
library('optparse')
library(ggplot2)
library(reshape)
library(plyr)
library(cowplot)
library(RColorBrewer)
library(wesanderson)
library(scales)
library(qqman)
library(HMMcopy)
library(gridExtra)
library(grid)

theme_set(theme_bw(base_size = 12))

out_dir = 'f:/Cornell/experiment/c.glabrata/manuscript/analysis/result/'

# pca
fly.drug = read.csv("f:/Cornell/experiment/c.glabrata/manuscript/analysis/scripts/fly.drug.csv",header=T)
pnas.drug = read.csv("f:/Cornell/experiment/c.glabrata/manuscript/analysis/scripts/pnas.path.csv",header=T)

plink.eigenv = read.table('f:/Cornell/experiment/c.glabrata/manuscript/analysis/scripts/Cgla.pca.eigenval')
pc1 = plink.eigenv[1]/sum(plink.eigenv)
qplot(1:nrow(plink.eigenv[1]),plink.eigenv[1])+labs(x="Eigen Value Ranking",y="Eigen Value")

plink.pca = read.table("f:/Cornell/experiment/c.glabrata/manuscript/analysis/scripts/Cgla.pca.eigenvec",header=F)
plink.pca = plink.pca[,c(1,3,4,5)]
colnames(plink.pca) = c('strain','PC1','PC2','PC3')
plink.pca = merge(plink.pca,fly.drug)
plink.pca = merge(plink.pca,pnas.drug,all.x=T)
plink.pca$fluconazole[is.na(plink.pca$fluconazole)]="Not Available"
plink.pca$fluconazole=factor(plink.pca$fluconazole,c('Not Available','14.4','43.2','129.6','388.8'))

#figure 2 add cluster number later in paint
cpalette = c("lightgrey", "#EABE94", "#0B775E", "#35274A", "#F2300F")
set.seed(20)
p=ggplot(plink.pca,aes(PC1,PC2,color=factor(fluconazole),shape=Clade))+geom_jitter(size=I(4),alpha=0.6,width = 0.02,height = 0.015)+xlab("PC1 (15.0%)")+ylab("PC2 (13.5%)")+
    scale_color_manual(name="Fluconazole MIC",values = cpalette)
                       #values = wes_palette(5,name='Rushmore'))
    
#ggplot(plink.pca,aes(PC1,PC2,color=factor(place2)))+geom_jitter(size=I(4),alpha=0.5,width = 0.01,height = 0.01)+xlab("PC1 (15.0%)")+ylab("PC2 (13.5%)")
ggsave(p,filename=paste0(out_dir,"/","48.pca.tiff"), width=174, height=120, units='mm', dpi=500)

# pca test
#drugtest=matrix(c(0, 6, 3.5, 3.5 ),nrow = 2,byrow = F)
#fisher.test(drugtest,alternative = 'l')

# exom diversity
exom10k = read.table("f:/Cornell/experiment/c.glabrata/manuscript/analysis/scripts/Cgla_exon_theta.10k",header = T)
#summary(lm(theta~exon,data=exom10k))
cor.test(exom10k$theta,exom10k$exon)
p=ggplot(exom10k,aes(exon,theta))+geom_point(alpha=0.2,color='darkblue')+geom_smooth(method = 'lm',se=F)+xlab("Exon Propotion")+ylab(expression(paste("Watterson's ",theta)))+
    annotate("text", x = c(0.15,0.15), y = c(11,12), label = c("rho = -0.17","p value = 2.14e-9"),size=6)
#ggsave(p,filename=paste0(out_dir,"/","exon.theta.tiff"), width=174, height=120, units='mm', dpi=500)

# diversity by chromosome
exom10k$chr=factor(exom10k$chr,c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13"))

p1=ggplot(exom10k,aes(start,pi,color=chr))+geom_line(size=I(0.5))+facet_grid(~chr,scales="free_x",space='free_x')+xlab("")+ylab(expression(paste("Watterson's ",theta)))+
    scale_x_continuous(breaks = seq(0, 1000000, 500000),labels=c('0k','500k','1000k'))+scale_color_discrete(guide=F)+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.text.x = element_text(size = 6),
          strip.background = element_blank(),
          panel.border = element_rect(colour = "gray"),
          axis.text.x = element_text(size = 4),axis.text.y = element_text(size = 8),
          axis.title.y =element_text(size=8),
          panel.spacing.x = unit(0.1, "mm"))
#ggsave(p,filename=paste0(out_dir,"/","theta.by.chr.tiff"), width=174, height=30, units='mm', dpi=500)


fly.drug=merge(fly.drug,pnas.drug[,c('name','fluconazole')],all.x = T)
fly.drug = fly.drug[order(fly.drug$strain),]
fly.drug = fly.drug[-1,]


### figure 3
cnv.cds = read.table("f:/Cornell/experiment/c.glabrata/manuscript/analysis/scripts/Cgla_CNV_05cds.out",header = F)
#df=gwas(Xa=cnv.data2,phe=fly.drug$survive,pc1=plink.pca$PC1,pc2=plink.pca$PC2,pc3=plink.pca$PC3,maf=0.05)
bin.cds.freq=(apply(cnv.cds[,-c(1,2,3)],1,function(x){sum(x)/length(x)}))
bin.cds.df = data.frame(chr=cnv.cds[,1],pos=cnv.cds[,2],pos2=cnv.cds[,3],freq=bin.cds.freq)
bin.cds.df$chr=factor(bin.cds.df$chr,c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13"))
bin.cds.df$chr=mapvalues(bin.cds.df$chr,from=c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13"),
          to=c("chrA","chrB","chrC","chrD","chrE","chrF","chrG","chrH","chrI","chrJ","chrK","chrL","chrM"))

exom10k$chr=mapvalues(exom10k$chr,from=c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13"),
                      to=c("chrA","chrB","chrC","chrD","chrE","chrF","chrG","chrH","chrI","chrJ","chrK","chrL","chrM"))

p1=ggplot(exom10k,aes(start,pi,color=chr))+geom_line(size=I(0.5))+facet_grid(~chr,scales="free_x",space='free_x')+xlab("")+ylab(expression(paste("Genetic Diversity (",pi,")")))+
    scale_x_continuous(breaks = seq(0, 1000000, 500000),labels=c('0k','500k','1000k'))+scale_color_discrete(guide=F)+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.text.x = element_text(size = 7),
          strip.background = element_blank(),
          panel.border = element_rect(colour = "gray"),
          axis.ticks.x = element_blank(),axis.text.x = element_blank(),axis.text.y = element_text(size = 8),
          axis.title.y =element_text(size=8),
          panel.spacing.x = unit(0.1, "mm"))

fmt_dcimals <- function(decimals=0){function(x) format(x,nsmall = decimals,scientific = FALSE)}
p2=ggplot(bin.cds.df,aes(pos,freq,color=chr))+geom_line(size=I(0.5))+facet_grid(~chr,scales="free_x",space='free_x')+xlab("")+ylab("Duplication Frequency")+
    scale_color_discrete(guide=F)+scale_y_continuous(labels = fmt_dcimals(3))+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.text.x = element_text(size = 7),
          strip.background = element_blank(),
          panel.border = element_rect(colour = "gray"),
          axis.ticks.x = element_blank(),axis.text.x = element_blank(),axis.text.y = element_text(size = 8),
          axis.title.y =element_text(size=8),
          panel.spacing.x = unit(0.1, "mm"))

cor.test(exom10k$pi,exom10k$exon)
p3=ggplot(exom10k,aes(exon,pi))+geom_point(alpha=0.2,color='darkblue')+geom_smooth(method = 'lm',se=F)+xlab("Exon Propotion")+ylab(expression(paste("Genetic Diversity (",pi,")")))+#ylab(expression(paste("Watterson's ",theta)))+
    annotate("text", x = c(0.2,0.2), y = c(0.0045,0.005), label = c("rho = -0.15","P value = 6.44e-8"),size=3)+
    theme(axis.title=element_text(size=9))

chrCNV = read.table("f:/Cornell/experiment/c.glabrata/manuscript/analysis/scripts/chromLength_CNV.txt",header=T)
summary(lm(chrCNV$Length~chrCNV$Total))
cor.test(chrCNV$Length,chrCNV$Total)
chrCNV$Lengthk = chrCNV$Length/1000
p4=ggplot(chrCNV,aes(Lengthk,Total))+geom_point()+geom_smooth(method = 'lm',se=F)+xlab("Chromosome Length (kb)")+
    annotate("text", x = 1000, y = 7.5, label = c("rho = -0.3329"),size=3)+ theme(axis.title=element_text(size=9))+
    ylab("Count of Large CNVs") + scale_y_continuous(breaks = seq(0,10,2))


#up_row <- plot_grid(p1, p2, labels = c('A', 'B'), align = 'v')
bottom_row <- plot_grid(p3, p4, labels = c('C', 'D'), align = 'h',scale = 0.95)
plot2=plot_grid(p1,p2, bottom_row, labels = c('A', 'B',''), ncol = 1,rel_heights = c(1,1,1.2),scale = 1)


#plot2=plot_grid(p1,p2,nrow = 2,labels = c('A',"B"),label_size = 6)
ggsave(plot2,filename=paste0(out_dir,"/","fig3.tiff"), width=174, height=160, units='mm', dpi=500)

bin.cds.high = bin.cds.df[bin.cds.df$freq>0.3,]
write.table(bin.cds.high,file=paste0(out_dir,"/cnv.0.3.bed"),row.names = F,col.names = F,quote = F)

# exonic diverstiy test
chisq.test(matrix(c(8045535,125739,12338305,261338),nrow=2))

# adhesin gene test
gff = read.csv("f:/Cornell/experiment/c.glabrata/manuscript/analysis/scripts/gff/C_glabrata_current.gene.csv",header=F)
gff.sub = gff[,c(9,10)]
colnames(gff.sub)=c("Gene","functions")

NSS = read.csv("f:/Cornell/experiment/c.glabrata/manuscript/analysis/scripts/countNS_gene.out.csv",header=T)
NSS = merge(NSS,gff.sub)

NSS = NSS[order(NSS$NS_Kb,decreasing = T),]
length(grep('adhesin',NSS$functions,ignore.case = T))
length(grep('adhesin',NSS$functions[1:100],ignore.case = T))

# total adhesin gene 68
fisher.test(matrix(c(5199,68,100,8),nrow = 2))
write.table(NSS[1:100,],file = "f:/Cornell/experiment/c.glabrata/manuscript/analysis/result/TableS2.xls",quote = F,row.names = F,col.names = F,sep='\t')

# Fig S Watterson's theta
p1=ggplot(exom10k,aes(start,theta,color=chr))+geom_line(size=I(0.5))+facet_grid(~chr,scales="free_x",space='free_x')+xlab("")+ylab("Watterson's theta")+
    scale_x_continuous(breaks = seq(0, 1000000, 500000),labels=c('0k','500k','1000k'))+scale_color_discrete(guide=F)+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.text.x = element_text(size = 7),
          strip.background = element_blank(),
          panel.border = element_rect(colour = "gray"),
          axis.ticks.x = element_blank(),axis.text.x = element_blank(),axis.text.y = element_text(size = 8),
          axis.title.y =element_text(size=8),
          panel.spacing.x = unit(0.1, "mm"))
ggsave(p1,filename=paste0(out_dir,"/sup/","figS_theta.tiff"), width=174, height=80, units='mm', dpi=500)


# gwas figure, fig7
xdata <- read.table("f:/Cornell/experiment/c.glabrata/manuscript/analysis/scripts/Cgla_both_fluco.gwas",header=T)
#xdata <- x[,c(1,2,3,11)]
#colnames(xdata) <- c("CHR", "SNP", "BP", "P");
xdata = xdata[xdata$CHR!=14,]
#write.table(xdata, file= "merged_suvival.gwas", quote=F, row.names=F)

#xdata$CHR=factor(xdata$CHR,c("1","2","3","4","5","6","7","8","9","10","11","12","13"))
#xdata$CHR=mapvalues(xdata$CHR,from=c("1","2","3","4","5","6","7","8","9","10","11","12","13"),
#                         to=c("chrA","chrB","chrC","chrD","chrE","chrF","chrG","chrH","chrI","chrJ","chrK","chrL","chrM"))


#manhattan(xdata, col=c("blue3", "orange4"), main = "Manhattan Plot", cex = 0.5, cex.axis = 0.8)
#manhattan(xdata,ylim=c(0,9),cex.axis = 1, col=brewer.pal(8, "Set2"))

tiff(filename=paste0(out_dir, '/',"manhattan_plot.tiff"),width=174, height=100, units="mm",res=500)
manhattan(xdata,ylim=c(0,8),cex.axis = 1, col=brewer.pal(8, "Set2"), suggestiveline=F, genomewideline=6.59,
          chrlabs = c("A","B","C","D","E","F","G","H","I","J","K","L","M"))
dev.off()

# significant signal
pick.sig = xdata[xdata$P<(0.05/nrow(xdata)),]

xdata[xdata$P<(0.05/nrow(xdata))*100,]
xdata[xdata$CHR==6 & xdata$P<(0.05/nrow(xdata))*1000,]

numofsig = function(chr,bp,l,xdata){
    xtemp = xdata[xdata$CHR==chr,]
    xnear = xtemp[xtemp$BP > (bp-l) & xtemp$BP < (bp+l),]
    nsnp = nrow(xnear)
    nsnp
}

x1=numofsig(6,824799,25000,xdata)
x2=numofsig(8,595532,25000,xdata)
x3=numofsig(11,283,25000,xdata)
x4=numofsig(12,1451204,25000,xdata)
x5=numofsig(13,1241240,25000,xdata)
x6=numofsig(13,10070,25000,xdata)
x1+x2+x3+x4+x5+x6

#qqplot
p.qq = xdata$P
o = -log10(sort(p.qq,decreasing=F))
e = -log10(ppoints(length(p.qq)))
tiff(filename=paste0(out_dir, '/',"qq_plot.tiff"),width=174, height=174, units="mm",res=500)
plot(e,o,xlab="Expected -log10(p)",ylab="Observed -log10(p)")
abline(0,1,col='red')
dev.off()


# 012 plot
gwas.012 = read.csv("f:/Cornell/experiment/c.glabrata/manuscript/analysis/scripts/chrI_490216_012.csv",header=T)
gwas.phe = merge(plink.pca,gwas.012)
ggplot(gwas.phe[gwas.phe$fluconazole!="Not Available",],aes(genotype,as.numeric(as.character(fluconazole))))+geom_jitter(height = 0,width = 0.15,size=2)+geom_smooth(method = 'lm')+
    ylim(0,200)

# CG figure fig4
dt <- read.table("f:/Cornell/experiment/c.glabrata/manuscript/analysis/scripts/MutRate.txt", header=T)
dt$trinuc <- factor(dt$trinuc, levels=c("AAA", "AAC", "AAG", "AAT", "CAA", "CAC", "CAG", "CAT", "GAA", "GAC", "GAG", "GAT", "TAA", "TAC", "TAG", "TAT", "ACA", "ACC", "ACG", "ACT", "CCA", "CCC", "CCG", "CCT", "GCA", "GCC", "GCG", "GCT", "TCA", "TCC", "TCG", "TCT"))
colvec <- ifelse(grepl("CG", dt$trinuc),"red", "black")

#theme_set(theme_bw(base_size=12))
#q<-qplot(data=dt, x=trinuc, y=rate2, geom=c("point"), color=type)
#q + theme(axis.text.x = element_text(size=8, angle=90),legend.position="none")
p = ggplot(data=dt, aes(x=trinuc, y=rate2, color=type)) +geom_point(size=3) + 
    theme(axis.text.x = element_text(size=8, angle=90, vjust = 0.5, colour=colvec),legend.position = 'top',legend.title = element_text(size=10))+
    labs(title="", x="", y="Relative Mutation Rate")+ 
    scale_color_discrete("Base at Second Position\n",labels = c("A/T", "C/G"))+guides(color=guide_legend(title.position="top"))
label.df = data.frame(trinuc = c("ACG","CCG",'GCG','TCG'), rate2 = c(1.75,2.25,1.9,1.9),type = c('c',"c",'c',"c"))
p=p + geom_text(data=label.df,label = "*",color = I('black'))
ggsave(p,filename=paste0(out_dir,"/fig4.tiff"), width=174, height=120, units='mm', dpi=500)


# CNV figure
grid_arrange_shared_legend <- function(..., ncol = length(list(...)), nrow = 1, position = c("bottom", "right")) {
    
    plots <- list(...)
    position <- match.arg(position)
    g <- ggplotGrob(plots[[1]] + theme(legend.position = position))$grobs
    legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
    lheight <- sum(legend$height)
    lwidth <- sum(legend$width)
    gl <- lapply(plots, function(x) x + theme(legend.position="none"))
    gl <- c(gl, ncol = ncol, nrow = nrow)
    
    combined <- switch(position,
                       "bottom" = arrangeGrob(do.call(arrangeGrob, gl),
                                              legend,
                                              ncol = 1,
                                              heights = unit.c(unit(1, "npc") - lheight, lheight)),
                       "right" = arrangeGrob(do.call(arrangeGrob, gl),
                                             legend,
                                             ncol = 2,
                                             widths = unit.c(unit(1, "npc") - lwidth, lwidth)))
    
    grid.newpage()
    grid.draw(combined)
    
    # return gtable invisibly
    invisible(combined)
    
}

# strain Y625, No14 Fig 6
column.types <- c("character", "numeric", "numeric", "numeric")
all.data <- read.table("f:/Cornell/experiment/c.glabrata/manuscript/analysis/scripts/14_cnv.seg", header=FALSE, colClasses=column.types)
colnames(all.data)<-c("Chrom","Start","End","CopyNum")
chr.data <- subset(all.data, Chrom=="chr4")

p1=ggplot(data=chr.data, aes(x= Start, y= CopyNum))+ geom_point(aes(color=CopyNum))+
    scale_x_continuous(labels = scales::unit_format("Kb", 1e-3), breaks=scales::pretty_breaks(n=3))+
    geom_vline(aes(xintercept= 600000), linetype="dashed", size=1.2)+ ylim(c(-1.5,2))+
    scale_colour_gradient2(name='log2(Copy Number)',limits=c(-2,2),low = 'purple', mid = "blue", high = "red")+ggtitle("Strain Y625 ChrD")+
    ylab("log2(Copy Number)")	+ xlab("") +theme(plot.title = element_text(hjust = 0.5))

column.types <- c("character", "numeric", "numeric", "numeric")
all.data <- read.table("f:/Cornell/experiment/c.glabrata/manuscript/analysis/scripts/11_cnv.seg", header=FALSE, colClasses=column.types)
colnames(all.data)<-c("Chrom","Start","End","CopyNum")

chr.data <- subset(all.data, Chrom=="chr4")

p2=ggplot(data=chr.data, aes(x= Start, y= CopyNum))+ geom_point(aes(color=CopyNum))+
    scale_x_continuous(labels = scales::unit_format("Kb", 1e-3), breaks=scales::pretty_breaks(n=3))+
    geom_vline(aes(xintercept= 600000), linetype="dashed", size=1.2)+ ylim(c(-1.5,2))+
    scale_colour_gradient2(name='log2(Copy Number)',limits=c(-2,2),low = "purple",mid="blue", high = "red")+ggtitle("Strain Y622 ChrD")+
    ylab("log2(Copy Number)")+ xlab("") +theme(plot.title = element_text(hjust = 0.5))

plot12=grid_arrange_shared_legend(p1, p2, ncol = 1, nrow = 2)
ggsave(plot12,filename=paste0(out_dir,"/","fig6B.tiff"), width=154, height=120, units='mm', dpi=500)

###HMMCopy adapted from Xuepeng 
##
in_dir = "f:/Cornell/experiment/c.glabrata/manuscript/analysis/scripts/"
Cagl_uncorrected_reads <- wigsToRangedData(paste0(in_dir,"14_readcounts.wig"), paste0(in_dir,"Cagl.gc.wig"), paste0(in_dir,"Cagl.map.wig"))
Cagl_corrected_copy <- correctReadcount(Cagl_uncorrected_reads)
param <- HMMsegment(Cagl_corrected_copy, getparam = TRUE)
param$mu <- log(c(1, 1.4, 2, 2.7, 3, 4.5) / 2, 2)
param$m <- param$mu
segmented_copy <- HMMsegment(Cagl_corrected_copy, param)
rangedDataToSeg(Cagl_corrected_copy, file = "Cagl_corrected_copy.seg")
dat<-read.table("Cagl_corrected_copy.seg",header=T,sep="\t")
chrOrder<-paste("chr",1:13,sep="")
dat$chr <-factor(dat$chr, levels=chrOrder)
newdata <-dat[order(dat$chr),]
newdata = newdata[newdata$chr!="chr14",]

newdata$chr=factor(newdata$chr,c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13"))
chr<-summary(newdata[,2])

newdata$chr=mapvalues(newdata$chr,from=c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13"),
                         to=c("chrA","chrB","chrC","chrD","chrE","chrF","chrG","chrH","chrI","chrJ","chrK","chrL","chrM"))

savedata=newdata[,2:5]

wcnd.df = data.frame(x=seq(1,length(newdata[,5])),value=newdata[,5])
p3=ggplot(wcnd.df,aes(x,value,color=value))+geom_point()+ylim(c(-2,2))+ylab("log2(Copy Number)")+xlab("")+
    scale_colour_gradient2(name='log2(Copy Number)',limits=c(-2,2),low = "purple",mid="blue", high = "red")+theme(legend.position ='none')+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text.x = element_blank(),axis.ticks.x = element_blank())

p3=p3+geom_vline(xintercept = 0,linetype=2)

st=0
p3=p3+annotate("text",x=-250,y=-Inf,label='Chr:',size=4,vjust=-1,hjust=0.8)
end<-st;st<-st+chr[[1]];pos<-(end+st)/2;p3=p3+geom_vline(xintercept = st,linetype=2);p3=p3+annotate("text",x=pos,y=-Inf,label=sub('chr','',as.character(newdata[st,2][1])),size=4,vjust=-1)
end<-st;st<-st+chr[[2]];pos<-(end+st)/2;p3=p3+geom_vline(xintercept = st,linetype=2);p3=p3+annotate("text",x=pos,y=-Inf,label=sub('chr','',as.character(newdata[st,2][1])),size=4,vjust=-1)
for(i in 3:length(chr))
{
    end<-st
    st<-st+chr[[i]]
    pos<-(end+st)/2
    p3=p3+geom_vline(xintercept = st,linetype=2)
    p3=p3+annotate("text",x=pos,y=-Inf,label=sub('chr','',as.character(newdata[st,2][1])),size=4,vjust=-1)
    
}
p3
ggsave(p3,filename=paste0(out_dir,"/","fig6A2.tiff"), width=174, height=80, units='mm', dpi=500)
         
p4=plot_grid(p3,plot12,nrow = 2,rel_heights = c(1,2))
ggsave(p4,filename=paste0(out_dir,"/","fig6.tiff"), width=174, height=180, units='mm', dpi=500)

###gwas for cnv
gwas=function(Xa,phe,pc1,pc2,pc3,maf=0.05){
    MAF = apply(Xa[,-c(1,2)],1,function(x){sum(x)/length(x)})
    Xa = Xa[which(MAF >= maf & MAF <= 1-maf),]
    N=nrow(Xa)
    pvalue=vector(length=N)
    F.vec=vector(length = N)
    for (i in 1:N){
        lm.sum=summary(lm(phe~as.numeric(Xa[i,-c(1,2)])+pc1+pc2+pc3))
        f=lm.sum$fstatistic
        F.vec[i]=f[1]
        pvalue[i]=pf(f[1],f[2],f[3],lower.tail=FALSE)
    }
    df=data.frame(CHR=Xa[,1],BP=Xa[,2],P=pvalue)
    return(df)
    #return(list(pvalue,lm.all))
}

df=gwas(Xa=cnv.data,phe=fly.drug$survive,pc1=plink.pca$PC1,pc2=plink.pca$PC2,pc3=plink.pca$PC3,maf=0.05)

df2=gwas(Xa=cnv.data,phe=fly.drug$fluconazole,pc1=plink.pca$PC1,pc2=plink.pca$PC2,pc3=plink.pca$PC3,maf=0.05)


cnv.data2 = read.table("f:/Cornell/experiment/c.glabrata/manuscript/analysis/scripts/Cgla_CNV_05.out",header = T)
df=gwas(Xa=cnv.data2,phe=fly.drug$survive,pc1=plink.pca$PC1,pc2=plink.pca$PC2,pc3=plink.pca$PC3,maf=0.05)

df2=gwas(Xa=cnv.data2,phe=fly.drug$fluconazole,pc1=plink.pca$PC1,pc2=plink.pca$PC2,pc3=plink.pca$PC3,maf=0.05)


# suplementary 
# Fig S1
draw.cnv = function(sup_dir,out_dir,note,name){
    Cagl_uncorrected_reads <- wigsToRangedData(paste0(sup_dir,note,"_readcounts.wig"), paste0(sup_dir,"Cagl.gc.wig"), paste0(sup_dir,"Cagl.map.wig"))
    Cagl_corrected_copy <- correctReadcount(Cagl_uncorrected_reads)
    param <- HMMsegment(Cagl_corrected_copy, getparam = TRUE)
    param$mu <- log(c(1, 1.4, 2, 2.7, 3, 4.5) / 2, 2)
    param$m <- param$mu
    segmented_copy <- HMMsegment(Cagl_corrected_copy, param)
    rangedDataToSeg(Cagl_corrected_copy, file = "Cagl_corrected_copy.seg")
    dat<-read.table("Cagl_corrected_copy.seg",header=T,sep="\t")
    chrOrder<-paste("chr",1:13,sep="")
    dat$chr <-factor(dat$chr, levels=chrOrder)
    newdata <-dat[order(dat$chr),]
    newdata = newdata[newdata$chr!="chr14",]
    chr<-summary(newdata[,2])
    
    newdata$chr=factor(newdata$chr,c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13"))
    newdata$chr=mapvalues(newdata$chr,from=c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13"),
                      to=c("chrA","chrB","chrC","chrD","chrE","chrF","chrG","chrH","chrI","chrJ","chrK","chrL","chrM"))

    wcnd.df = data.frame(x=seq(1,length(newdata[,5])),value=newdata[,5])
    p3=ggplot(wcnd.df,aes(x,value,color=value))+geom_point()+ylim(c(-2,2))+ylab("log2(Copy Number)")+xlab("")+
        scale_colour_gradient2(name='log2(Copy Number)',limits=c(-2,2),low = "purple",mid="blue", high = "red")+theme(legend.position ='none')+
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text.x = element_blank(),axis.ticks.x = element_blank(),plot.title = element_text(hjust = 0.5))+
        ggtitle(name)

    p3=p3+geom_vline(xintercept = 0,linetype=2)
    #chr<-summary(newdata[,2])
    st=0
    p3=p3+annotate("text",x=-250,y=-Inf,label='Chr:',size=4,vjust=-1,hjust=0.8)
    end<-st;st<-st+chr[[1]];pos<-(end+st)/2;p3=p3+geom_vline(xintercept = st,linetype=2);p3=p3+annotate("text",x=pos,y=-Inf,label=sub('chr','',as.character(newdata[st,2][1])),size=4,vjust=-1)
    end<-st;st<-st+chr[[2]];pos<-(end+st)/2;p3=p3+geom_vline(xintercept = st,linetype=2);p3=p3+annotate("text",x=pos,y=-Inf,label=sub('chr','',as.character(newdata[st,2][1])),size=4,vjust=-1)
    for(i in 3:length(chr)){
        end<-st
        st<-st+chr[[i]]
        pos<-(end+st)/2
        p3=p3+geom_vline(xintercept = st,linetype=2)
        p3=p3+annotate("text",x=pos,y=-Inf,label=sub('chr','',as.character(newdata[st,2][1])),size=4,vjust=-1)
    
    }
    
    ggsave(p3,filename=paste0(out_dir,"/sup/",name,".tiff"), width=174, height=100, units='mm', dpi=500)
    return(p3)
}

# Y625 14; Y626 15; Y627 16; Y641 18; Y645 21; Y648 24
sup_dir = "f:/Cornell/experiment/c.glabrata/manuscript/analysis/scripts/wig/"

Y625=draw.cnv(sup_dir = sup_dir,out_dir = out_dir,note = 14,name='Y625')
Y626=draw.cnv(sup_dir = sup_dir,out_dir = out_dir,note = 15,name='Y626')
Y627=draw.cnv(sup_dir = sup_dir,out_dir = out_dir,note = 16,name='Y627')
Y641=draw.cnv(sup_dir = sup_dir,out_dir = out_dir,note = 18,name='Y641')
Y645=draw.cnv(sup_dir = sup_dir,out_dir = out_dir,note = 21,name='Y645')
Y648=draw.cnv(sup_dir = sup_dir,out_dir = out_dir,note = 24,name='Y648')

plotFS1=grid_arrange_shared_legend(Y625,Y626,Y627,Y641,Y645,Y648, ncol = 1, nrow = 6)
ggsave(plotFS1,filename=paste0(out_dir,"/sup/","FigS1.tiff"), width=174, height=300, units='mm', dpi=500)


# Fig S2
draw.chr = function(sup_dir,out_dir,note,name,chr){
    column.types <- c("character", "numeric", "numeric", "numeric")
    all.data <- read.table(paste0(sup_dir,note,"_cnv.seg"), header=FALSE, colClasses=column.types)
    colnames(all.data)<-c("Chrom","Start","End","CopyNum")
    chr.data <- subset(all.data, Chrom==chr)
    
    p=ggplot(data=chr.data, aes(x= Start, y= CopyNum))+ geom_point(aes(color=CopyNum))+
        scale_x_continuous(labels = scales::unit_format("Kb", 1e-3), breaks=scales::pretty_breaks(n=3))+
        ylim(c(-1.5,2))+
        scale_colour_gradient2(name='log2(Copy Number)',limits=c(-2,2),low = 'purple', mid = "blue", high = "red")+ggtitle(name)+
        ylab("log2(Copy Number)")	+ xlab("") +theme(plot.title = element_text(hjust = 0.5))
    return(p)
}


# chrD 4: Y622 11; Y625 14; Y1640 41; Y1641 42; Y1642 43; Y1644 45
Y622=draw.chr(sup_dir = sup_dir,out_dir = out_dir,note=11,name='Y622',chr="chr4")
Y625=draw.chr(sup_dir = sup_dir,out_dir = out_dir,note=14,name='Y625',chr="chr4")
Y1640=draw.chr(sup_dir = sup_dir,out_dir = out_dir,note=41,name='Y1640',chr="chr4")
Y1641=draw.chr(sup_dir = sup_dir,out_dir = out_dir,note=42,name='Y1641',chr="chr4")
Y1642=draw.chr(sup_dir = sup_dir,out_dir = out_dir,note=43,name='Y1642',chr="chr4")
Y1644=draw.chr(sup_dir = sup_dir,out_dir = out_dir,note=45,name='Y1644',chr="chr4")
plotA=grid_arrange_shared_legend(Y622,Y625,Y1640,Y1641,Y1642,Y1644, ncol = 2, nrow = 3)
ggsave(plotA,filename=paste0(out_dir,"/sup/","FigS2A.tiff"), width=174, height=200, units='mm', dpi=500)

# chrE 5: Y621 10; Y622 11; Y624 13; Y1641 42
Y621=draw.chr(sup_dir = sup_dir,out_dir = out_dir,note=10,name='Y621',chr="chr5")
Y622=draw.chr(sup_dir = sup_dir,out_dir = out_dir,note=11,name='Y622',chr="chr5")
Y624=draw.chr(sup_dir = sup_dir,out_dir = out_dir,note=13,name='Y624',chr="chr5")
Y1641=draw.chr(sup_dir = sup_dir,out_dir = out_dir,note=42,name='Y1641',chr="chr5")
plotA=grid_arrange_shared_legend(Y621,Y622,Y624,Y1641, ncol = 2, nrow = 2)
ggsave(plotA,filename=paste0(out_dir,"/sup/","FigS2B.tiff"), width=174, height=150, units='mm', dpi=500)


# chrL 12: Y622 11; Y645 21; Y656 31
Y622=draw.chr(sup_dir = sup_dir,out_dir = out_dir,note=11,name='Y622',chr="chr12")
Y645=draw.chr(sup_dir = sup_dir,out_dir = out_dir,note=21,name='Y645',chr="chr12")
Y656=draw.chr(sup_dir = sup_dir,out_dir = out_dir,note=31,name='Y656',chr="chr12")
plotA=grid_arrange_shared_legend(Y622,Y645,Y656, ncol = 2, nrow = 2)
ggsave(plotA,filename=paste0(out_dir,"/sup/","FigS2C.tiff"), width=174, height=150, units='mm', dpi=500)



