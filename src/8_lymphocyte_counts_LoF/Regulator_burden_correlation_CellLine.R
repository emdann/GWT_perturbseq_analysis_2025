args   <- commandArgs(trailingOnly = T)  
options("scipen"=10) 
FILE<-as.character( args[1] ) ##exp. "Backman_2021_86.per_gene_estimates.tsv"
FILE_perturb<-as.character( args[2] ) ##exp. "K562GW_limma_logFC_sum.txt"

LOF<-read.table(paste0("input/", FILE), sep="\t", quote="", header=T, stringsAsFactor=F)
LOF$post_mean[is.element(LOF$post_mean, "Inf")]<-max(LOF$post_mean[!is.infinite(LOF$post_mean)])
LOF$post_mean[is.element(LOF$post_mean, "-Inf")]<-min(LOF$post_mean[!is.infinite(LOF$post_mean)])

library(data.table)
df1<-fread(paste0(FILE2), data.table=F, header=T)
row.names(df1)<-df1[,1]
df1<-df1[,-1]

corresp<-read.table("input/gencode_v41_gname_gid_ALL_sorted", header=F, stringsAsFactor=F)
corresp[,1]<-sapply(strsplit(as.character(corresp[,1]), "\\."), function(x){x[1]})
corresp<-corresp[!duplicated(corresp[,1]),]
corresp<-corresp[!duplicated(corresp[,2]),]
row.names(corresp)<-corresp[,2]
library("boot")

shet<-read.table("input/shet_10bins.txt", header=T, stringsAsFactor=F)
shet<-shet[is.element(shet$ensg, LOF$ensg),]

summary<-data.frame()
for(i in 1:nrow(df1)){

tmp<-df1[i,,drop=F]
tmp<-t(tmp)
tmp<-data.frame(gene=row.names(tmp), tmp)
colnames(tmp)[2]<-"perturb_beta"
tmp<-data.frame(tmp, ensg=corresp[tmp$gene, 1])

df<-merge(tmp, LOF, by="ensg")

target=row.names(df1)[i]
df<-df[!is.element(df$ensg, target),]

##with shet regression
df<-merge(df, shet, by="ensg")
df$post_mean<-scale(df$post_mean)
df$perturb_beta<-scale(df$perturb_beta)

fit<- lm(post_mean~perturb_beta + shet, data=df)
P_withShet<-summary(fit)$coefficients[2,4]
beta_withShet<-summary(fit)$coefficients[2,1]
betaSE_withShet<-summary(fit)$coefficients[2,2]


P_pearson<-cor.test(df$perturb_beta, df$post_mean, method="pearson")$p.value
R_pearson<-cor.test(df$perturb_beta, df$post_mean, method="pearson")$estimate
R_pearson_CIlower<-cor.test(df$perturb_beta, df$post_mean, method="pearson")$conf.int[1]
R_pearson_CIupper<-cor.test(df$perturb_beta, df$post_mean, method="pearson")$conf.int[2]


hoge<-data.frame(ensg=row.names(df1)[i], P_withShet=P_withShet, beta_withShet=beta_withShet, betaSE_withShet=betaSE_withShet, P_pearson=P_pearson, R_pearson=R_pearson, R_pearson_CIlower=R_pearson_CIlower,R_pearson_CIupper=R_pearson_CIupper )
summary<-rbind(summary, hoge)
}

corresp<-read.table("input/gencode_v41_gname_gid_ALL_sorted", header=F, stringsAsFactor=F)
corresp[,1]<-sapply(strsplit(as.character(corresp[,1]), "\\."), function(x){x[1]})
corresp<-corresp[!duplicated(corresp[,1]),]
corresp<-corresp[!duplicated(corresp[,2]),]
row.names(corresp)<-corresp[,1]

summary<-data.frame(gene=corresp[as.character(summary$ensg), 2], summary)

write.table(summary, paste0("BurdenRegCor/GeneLevel/", FILE2, "_", FILE, "_geneRegulation_correlation.txt"), row.names=F, sep="\t", quote=F)
