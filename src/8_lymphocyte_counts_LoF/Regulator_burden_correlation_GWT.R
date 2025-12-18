args   <- commandArgs(trailingOnly = T)  
options("scipen"=10) 
FILE<-as.character( args[1] ) ##exp. "Backman_2021_86.per_gene_estimates.tsv"

LOF<-read.table(paste0("input/", FILE), sep="\t", quote="", header=T, stringsAsFactor=F)
LOF$post_mean[is.element(LOF$post_mean, "Inf")]<-max(LOF$post_mean[!is.infinite(LOF$post_mean)])
LOF$post_mean[is.element(LOF$post_mean, "-Inf")]<-min(LOF$post_mean[!is.infinite(LOF$post_mean)])

library(data.table)
df<-fread(paste0("CD4i_final_merged.DE_pseudobulk_logFC.csv"), data.table=F, header=T)
df<-df[,-1]
df<-t(df)

pert<-read.csv("CD4i_final_merged.DE_pseudobulk_metadata.csv", header=T, stringsAsFactor=F)
colnames(df)<-paste0(pert$culture_condition, "_", pert$target_contrast_gene_name)

genes<-read.csv("CD4i_final_merged.DE_pseudobulk_genes.csv", header=T, stringsAsFactor=F)
row.names(df)<-genes$gene_ids

library("boot")

shet<-read.table("input/shet_10bins.txt", header=T, stringsAsFactor=F)
shet<-shet[is.element(shet$ensg, LOF$ensg),]


for(COND in c("Rest",  "Stim8hr", "Stim48hr")){
df1<-df[,grep(COND, colnames(df))]
colnames(df1)<-sapply(strsplit(colnames(df1), "_"), function(x){x[2]})
corresp<-read.table("input/gencode_v41_gname_gid_ALL_sorted", header=F, stringsAsFactor=F)
corresp<-corresp[!duplicated(corresp[,1]),]
corresp<-corresp[!duplicated(corresp[,2]),]
row.names(corresp)<-corresp[,2]

summary<-data.frame()
for(i in 1:nrow(df1)){

tmp<-df1[i,,drop=F]
tmp<-t(tmp)
tmp<-data.frame(gene=row.names(tmp), tmp)
colnames(tmp)[2]<-"perturb_beta"
tmp<-data.frame(tmp, ensg=corresp[tmp$gene, 1])

df2<-merge(tmp, LOF, by="ensg")

target=row.names(df1)[i]
df2<-df2[!is.element(df2$ensg, target),]
df2<-df2[!is.na(df2$perturb_beta),]

if(nrow(df2)>100){
##with shet regression
df2<-merge(df2, shet, by="ensg")
df2$post_mean<-scale(df2$post_mean)
df2$perturb_beta<-scale(df2$perturb_beta)

fit<- lm(post_mean~perturb_beta + shet, data=df2)
P_withShet<-summary(fit)$coefficients[2,4]
beta_withShet<-summary(fit)$coefficients[2,1]
betaSE_withShet<-summary(fit)$coefficients[2,2]


P_pearson<-cor.test(df2$perturb_beta, df2$post_mean, method="pearson")$p.value
R_pearson<-cor.test(df2$perturb_beta, df2$post_mean, method="pearson")$estimate
R_pearson_CIlower<-cor.test(df2$perturb_beta, df2$post_mean, method="pearson")$conf.int[1]
R_pearson_CIupper<-cor.test(df2$perturb_beta, df2$post_mean, method="pearson")$conf.int[2]


hoge<-data.frame(ensg=row.names(df1)[i], P_withShet=P_withShet, beta_withShet=beta_withShet, betaSE_withShet=betaSE_withShet, P_pearson=P_pearson, R_pearson=R_pearson, R_pearson_CIlower=R_pearson_CIlower,R_pearson_CIupper=R_pearson_CIupper )
summary<-rbind(summary, hoge)
}
}


corresp<-read.table("input/gencode_v41_gname_gid_ALL_sorted", header=F, stringsAsFactor=F)
corresp<-corresp[!duplicated(corresp[,1]),]
corresp<-corresp[!duplicated(corresp[,2]),]
row.names(corresp)<-corresp[,1]

summary<-data.frame(gene=corresp[as.character(summary$ensg), 2], summary)

write.table(summary, paste0("BurdenRegCor/GeneLevel/", COND, "_", FILE, "_geneRegulation_correlation.txt"), row.names=F, sep="\t", quote=F)
}
