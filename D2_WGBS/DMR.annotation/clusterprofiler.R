
library("ChIPseeker")
library("tidyverse")
library("DOSE")
library("clusterProfiler") 
library("org.Mm.eg.db") 
library(TxDb.Mmusculus.UCSC.mm10.knownGene) 
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene 

dat_M_hyper <- readPeakFile("dmrs_D2_M_WTvsCre.delta0.p1e-3.0705.2023.txt.bedgraph.hyper.bedgraph")
dat_M_hypo <- readPeakFile("dmrs_D2_M_WTvsCre.delta0.p1e-3.0705.2023.txt.bedgraph.hypo.bedgraph")

dat_F_hyper <- readPeakFile("dmrs_D2_F_WTvsCre_p1e-3.delta0.all.txt.bedgraph.hyper.bedgraph")
dat_F_hypo <- readPeakFile("dmrs_D2_F_WTvsCre_p1e-3.delta0.all.txt.bedgraph.hypo.bedgraph")
dat_WT_hyper <- readPeakFile("dmrs_D2_WT_MvsF.delta0.p1e-3.0705.2023.txt.bedgraph.hyper.bedgraph")
dat_WT_hypo <- readPeakFile("dmrs_D2_WT_MvsF.delta0.p1e-3.0705.2023.txt.bedgraph.hypo.bedgraph")
dat_KO_hyper <- readPeakFile("dmrs_D2_KO_MvsF.delta0.p1e-3.all.txt.bedgraph.hyper.bedgraph")
dat_KO_hypo <- readPeakFile("dmrs_D2_KO_MvsF.delta0.p1e-3.all.txt.bedgraph.hypo.bedgraph")



peaks <- list(M_hyper=dat_M_hyper, M_hypo=dat_M_hypo, F_hyper=dat_F_hyper, F_hypo=dat_F_hypo, WT_hyper=dat_WT_hyper, WT_hypo=dat_WT_hypo,  KO_hyper=dat_KO_hyper, KO_hypo = dat_KO_hypo )
peakAnnoList <- lapply(peaks, annotatePeak, TxDb=txdb,tssRegion=c(-1000, 100), verbose=FALSE,addFlankGeneInfo=TRUE, flankDistance=5000,annoDb="org.Mm.eg.db")
dat_gene_list = lapply(peakAnnoList, function(i) as.data.frame(i)$geneId)


dat_geneName_list = lapply(peakAnnoList, function(i) as.data.frame(i)$SYMBOL)
all_name  =dat_geneName_list %>%flatten()%>%as_vector()%>%unique()

write.table(all_name,"all_DMR_gene.name.txt",quote=F,row.names = F)

library(writexl)
dat_geneName_list = lapply(peakAnnoList, function(i) as.data.frame(i))
write_xlsx(dat_geneName_list, "my_list.xlsx")

#GO_1 <- enrichGO(gene = dat_gene_list$M_DMR, keyType = "ENTREZID", OrgDb = org.Mm.eg.db, ont = "BP", pAdjustMethod = "BH", qvalueCutoff = 0.05, readable = TRUE)
#GO_2 <- enrichGO(gene = dat_gene_list$F_DMR, keyType = "ENTREZID", OrgDb = org.Mm.eg.db, ont = "BP", pAdjustMethod = "BH", qvalueCutoff = 0.05, readable = TRUE)
#GO_1_simp <- simplify(GO_1, cutoff = 0.4, by = "p.adjust",select_fun = min, measure = "Wang", semData = NULL)
#GO_2_simp <- simplify(GO_2, cutoff = 0.4, by = "p.adjust",select_fun = min, measure = "Wang", semData = NULL)



#GO_1_simp <- simplify(GO_1, cutoff = 0.7, by = "p.adjust",select_fun = min, measure = "Wang", semData = NULL)



compGO_BP <- compareCluster(dat_gene_list,fun="enrichGO",OrgDb = org.Mm.eg.db, ont = "BP",pvalueCutoff=0.05,qvalueCutoff=0.05)
ComB <- function(x){paste(bitr(unlist(strsplit(x,"/")), fromType = "ENTREZID",toType = c("ENSEMBL","SYMBOL"), OrgDb = org.Mm.eg.db)$SYMBOL,collapse = "/")}
GO_cluster_summary <- as.data.frame(compGO_BP)
GO_cluster_summary$Symbol <- sapply(GO_cluster_summary$geneID,FUN=ComB)
write.table(GO_cluster_summary, "all_DMR_GO.BP.txt",sep="\t",quote=F,row.names = F)

compGO_CC <- compareCluster(dat_gene_list,fun="enrichGO",OrgDb = org.Mm.eg.db, ont = "CC",pvalueCutoff=0.05,qvalueCutoff=0.05)
compGO_MF <- compareCluster(dat_gene_list,fun="enrichGO",OrgDb = org.Mm.eg.db, ont = "MF",pvalueCutoff=0.05,qvalueCutoff=0.05)

#compGO_simp_BP <- simplify(compGO_BP, cutoff = 0.5, by = "p.adjust", select_fun = min, measure = "Wang", semData = NULL)#
#compGO_simp_CC <- simplify(compGO_CC, cutoff = 0.5, by = "p.adjust", select_fun = min, measure = "Wang", semData = NULL)
#compGO_simp_MF <- simplify(compGO_MF, cutoff = 0.5, by = "p.adjust", select_fun = min, measure = "Wang", semData = NULL)

#p5.1 <- dotplot(compGO_simp,showCategory = 50, title = "GO Analysis")
p6.1 <- dotplot(compGO_BP,size= "p.adjust", showCategory = 10, includeAll=T)+ scale_color_gradient(high="gold", low="red")

##### manually plot
BP_res = compGO_BP@compareClusterResult
Top5_terms <- BP_res %>%
  group_by(Cluster) %>%
  slice_min(order_by = pvalue, n = 5) %>%
  ungroup() %>% dplyr::select(Description)

Top5_res <- BP_res %>% filter(Description %in% Top5_terms$Description) 

Top5_res$DMR = sapply(strsplit(as.character(Top5_res$Cluster), split ="_"), `[`, 2)


Top5_res <- Top5_res %>%
  arrange(Cluster, p.adjust) %>%
  mutate(Description = factor(Description, levels = unique(Description)))

Top5_res$Description <- factor(Top5_res$Description, levels = rev(levels(Top5_res$Description)))


p_tf2 <- ggplot(Top5_res, aes(x=Cluster,y=Description, size= -log(p.adjust))) +
geom_point(aes(fill=DMR), color="black", pch=21, stroke=0.5 ) + 
scale_size(range=c(1,2)) +
xlab("") + 
# coord_flip() +
scale_fill_manual(values=c("hyper"="red", "hypo"="blue")) 


p_tf.2 = p_tf2 + theme_bw(base_size = 5)+ theme(plot.margin = unit(c(0,0,0,0), "cm"),legend.title= element_blank(),
                                                     axis.text.x =element_text(size = 5,angle = 90, vjust=0),
                                                      #axis.text.x=element_blank(),
                                                     #legend.position =c(0.3,0.7),
                                                    axis.text.y =element_text(size = 5, vjust=0),
                                                     legend.position ="top", 
                                                     legend.text = element_text(size=5),
                                                     legend.key.size = unit(0.2, 'cm'),
                                                     legend.key.height = unit(0.2, 'cm'),
                                                     legend.key.width = unit(0.2, 'cm'),
                                                     legend.box = "horizontal") 

ggsave("DMR.GO.BP.2.svg",plot=p_tf.2, device="svg",width=8, height=7,unit="cm") 

source("D:/yli_R_lib/plot_top_n_terms.r")

p1 = create_top_terms_plot(BP_res,"BP_Comp", 5)

p6.1 <- dotplot(compGO_BP,size= "p.adjust", showCategory = 3, includeAll=T)


p7.1 <- dotplot(compGO_CC,size="p.adjust", showCategory = 10,includeAll=T)+ scale_color_gradient(high="gold", low="red")
CC_res = compGO_CC@compareClusterResult
p2 = create_top_terms_plot(CC_res,"CC_Comp", 5)

p8.1 <- dotplot(compGO_MF,size="p.adjust", showCategory = 10, includeAll=T)+ scale_color_gradient(high="gold", low="red")

p3 = create_top_terms_plot(compGO_MF@compareClusterResult,"MF_Comp", 5)

ggsave("MF.DMR.GO_BP.top5.1.svg",plot=p1, device="svg",width=4, height=4)
ggsave("MF.DMR.GO_CC.top5.1.svg",plot=p2, device="svg",width=4, height=4)
ggsave("MF.DMR.GO_MF.top5.1.svg",plot=p3, device="svg",width=4, height=4)

#compGO_BP <- compareCluster(dat_gene_list,fun="enrichGO",OrgDb = org.Mm.eg.db, ont = "BP",pvalueCutoff=0.1,qvalueCutoff=0.1)
#p6.11<- dotplot(compGO_BP,showCategory = 5,includeAll=T)+ scale_color_gradient(high="gold", low="red")
#ggsave("MF.DMR.GO_BP_100.top10.p0.1.svg",plot=p6.11, device="svg", width=8, height=8)

#p1.1 <- dotplot(GO_2_simp, showCategory = 50, title = "GO Analysis")
#p2.1 <- dotplot(GO_3_simp,showCategory = 50, title = "GO Analysis")
#p3.1 <- dotplot(GO_4_simp,showCategory = 50, title = "GO Analysis")

compKEGG <- compareCluster(geneCluster = dat_gene_list, fun = "enrichKEGG", organism = "mouse", pvalueCutoff = 0.05, pAdjustMethod = "BH")
p4.1 <- dotplot(compKEGG, showCategory = 20, title = "KEGG Pathway Enrichment Analysis",includeAll=T)+ scale_color_gradient(high="gold", low="red")
ggsave("MF.DMR.KEGG.svg",plot=p4.1, device="svg",width=12, height=10)
p4 = create_top_terms_plot(compKEGG@compareClusterResult,"MF_Comp", 5)
ggsave("MF.DMR.KEGG.top5.1.svg",plot=p4, device="svg",width=4, height=4)




ComB <- function(x){paste(bitr(unlist(strsplit(x,"/")), fromType = "ENTREZID",toType = c("ENSEMBL","SYMBOL"), OrgDb = org.Mm.eg.db)$SYMBOL,collapse = "/")}
GO_cluster_summary <- as.data.frame(compGO_BP)
GO_cluster_summary$Symbol <- sapply(GO_cluster_summary$geneID,FUN=ComB)
write.table(GO_cluster_summary, "MF.DMR.GO.BP_simplified_0.2.txt",sep="\t",quote=F,row.names = F)

GO_cluster_summary <- as.data.frame(compGO_CC)
GO_cluster_summary$Symbol <- sapply(GO_cluster_summary$geneID,FUN=ComB)
write.table(GO_cluster_summary, "MF.DMR.GO.CC_simplified_0.2.txt",sep="\t",quote=F,row.names = F)

GO_cluster_summary <- as.data.frame(compGO_MF)
GO_cluster_summary$Symbol <- sapply(GO_cluster_summary$geneID,FUN=ComB)
write.table(GO_cluster_summary, "MF.DMR.GO.MF_simplified_0.2.txt",sep="\t",quote=F,row.names = F)

KEGG_cluster_summary <- as.data.frame(compKEGG)
KEGG_cluster_summary$Symbol <- sapply(KEGG_cluster_summary$geneID,FUN=ComB)
write.table(KEGG_cluster_summary, "MF.DMR.KEGG.txt",sep="\t",quote=F,row.names = F)




# plot manually

#p5 = ggplot(dat_tf_2, aes(x=Group,y=TF,size=(20-rank), col=-log(p.value))) +
##  geom_point() +scale_color_gradient(low = "blue", high = "red")+ 
#  theme(axis.text.x = element_text(face="bold", color="black", angle=45,hjust=1,size=15), axis.text.y=element_text(size=15))+
#  labs(x="",y="",title="UMR-enriched TF ")
#p_test = ggplot(compKEGG@compareClusterResult,  aes(x=Cluster,y=Description, size=GeneRatio, col=-log(pvalue))) + geom_point() +scale_color_gradient(low = "blue", high = "red")


#df = compKEGG@compareClusterResult
#dim(df)

##df[nrow(df)+1,] <-df[20,]
#df[nrow(df)+1,] <-df[20,]
#df[21,1] <- "F_DMR"
#df[22,1] <- "KO_DMR"

# idx <- order(df[["qvalue"]], decreasing = F)
#    df$Description <- factor(df$Description,
#                          levels=rev(unique(df$Description[idx])))
    #
#   g_test =  ggplot(df, aes(x=Cluster, y=Description, size=GeneRatio, color=qvalue)) +
#        geom_point() + scale_color_continuous(low="red", high="blue", name ="color", guide=guide_colorbar(reverse=TRUE))   + ylab(NULL)   +
#        scale_size(range=c(1, 2))+ scale_y_discrete()
#  g_test.1 = g_test + theme_classic(base_size = 5)
#ggsave("KEGG_enrichment.2.svg",plot=g_test.1, device="svg",width=7, height=9,units = "cm")

p6.2 = p6.1+ scale_size(range=c(1, 3)) + theme(text = element_text(size = 5),
                                               axis.text.x = element_text(angle = 90, hjust = 1, size=5,color="black"),
                                               axis.text.y = element_text( size=5,color="black"),
                                              legend.title=element_text(size=5), 
                                              panel.grid.major=element_line(colour=NA),
                                               legend.text = element_text(size=5),
                                               legend.key.size = unit(0.2, 'cm'),
                                                legend.key.height = unit(0.2, 'cm'),
                                                 legend.key.width = unit(0.2, 'cm')
  
)

ggsave("KEGG_enrichment.GO.BP.5.svg",plot=p6.2, device="svg",width=8, height=8,units = "cm")

p7.2 = p7.1 + scale_size(range=c(1, 2)) + theme(text = element_text(size = 5),
                                               axis.text.x = element_text(angle = 90, hjust = 1, size=5,color="black"),
                                               axis.text.y = element_text( size=5,color="black"),
                                              legend.title=element_text(size=5), 
                                              panel.grid.major=element_line(colour=NA),
                                               legend.text = element_text(size=5),
                                               legend.key.size = unit(0.2, 'cm'),
                                                legend.key.height = unit(0.2, 'cm'),
                                                 legend.key.width = unit(0.2, 'cm')
  
)

ggsave("KEGG_enrichment.GO.CC.2.svg",plot=p7.2, device="svg",width=10, height=10,units = "cm")


p8.2 = p8.1 + scale_size(range=c(1, 2)) + theme(text = element_text(size = 5),
                                               axis.text.x = element_text(angle = 90, hjust = 1, size=5,color="black"),
                                               axis.text.y = element_text( size=5,color="black"),
                                              legend.title=element_text(size=5), 
                                              panel.grid.major=element_line(colour=NA),
                                               legend.text = element_text(size=5),
                                               legend.key.size = unit(0.2, 'cm'),
                                                legend.key.height = unit(0.2, 'cm'),
                                                 legend.key.width = unit(0.2, 'cm')
)

ggsave("KEGG_enrichment.GO.MF.8.svg",plot=p8.2, device="svg",width=8, height=8, units = "cm")






library(pathview)
        
hsa04110 <- pathview(gene.data  = dat_gene_list$M_DMR,
                     pathway.id = "hsa04110",
                     species    = "hsa")



#### reactome
library(ReactomePA)

compGO_RA <- compareCluster(dat_gene_list,fun="enrichPathway",organism = "mouse",pvalueCutoff=0.05,qvalueCutoff=0.05)
WT <- enrichPathway(gene=dat_gene_list$M_DMR,organism = "mouse",pvalueCutoff = 0.05, readable=TRUE)WT
all =dat_gene_list %>%flatten()%>%as_vector()%>%unique()
compGO_RA_all <- enrichPathway(gene = all,organism = "mouse",pvalueCutoff = 0.05, readable=TRUE)
res=as.data.frame(compGO_RA_all@result)
res$Description = str_replace_all(res$Description,pattern = "\r",replacement = "")

write.table(res, "all_DMR.reactome_simplified_0.2.txt",sep="\t", quote=F,row.names = F)

dat_DMR_cDMR = read.table("DMR_cDMR.genes.txt",header=F)
dat_DMR_cDMR_ID=bitr(dat_DMR_cDMR$V1,fromType = "SYMBOL",toType = "ENTREZID", OrgDb = "org.Mm.eg.db",drop=F )

compGO_RA_all_2 <- enrichPathway(gene = dat_DMR_cDMR_ID$ENTREZID,organism = "mouse",pvalueCutoff = 0.05, readable=TRUE)
res=as.data.frame(compGO_RA_all_2@result)
res$Description = str_replace_all(res$Description,pattern = "\r",replacement = "")

write.table(res, "all_DMR_cDMR.reactome_simplified_0.2.txt",sep="\t", quote=F,row.names = F)


#### save the results
save(compGO_BP, compGO_CC, compGO_MF, compKEGG ,compGO_RA , file="GO_KEGG.rData")



#GO simplify
library("simplifyEnrichment")
  

#GO_clr1 = clr1_anno$result$term_id[grep("GO",clr1_anno$result$term_id)]
#all_GO = compGO_BP@compareClusterResult$ID
set.seed(888)
#go_id = random_GO(500)
#mat = GO_similarity(all_GO, ont = "BP")

#cluster_by_kmeans(mat, max_k = max(2, min(round(nrow(mat)/5), 100)))

#svg("all_DMR_GO_enrich.cluster.svg",width = 20,height = 10)
##plot_binary_cut(mat, value_fun = area_above_ecdf, cutoff = 0.7,partition_fun = partition_by_pam, dend = NULL, dend_width = unit(3, "cm"), depth = NULL, show_heatmap_legend = TRUE)
#dev.off()


GO_clr2 = clr2_anno$result$term_id[grep("GO",clr2_anno$result$term_id)]
set.seed(888)
#go_id = random_GO(500)
mat = GO_similarity(GO_clr2, ont = "BP")
cluster_by_kmeans(mat, max_k = max(2, min(round(nrow(mat)/5), 100)))
svg("GO_enrich.cluster.svg",width = 20,height = 10)
plot_binary_cut(mat, value_fun = area_above_ecdf, cutoff = 0.85,partition_fun = partition_by_pam, dend = NULL, dend_width = unit(3, "cm"), depth = NULL, show_heatmap_legend = TRUE)
dev.off()

svg("GO_enrich.cluster.2.svg",width = 30,height = 10)
df = simplifyGO(mat)
dev.off()

#David_4 <- enrichDAVID(gene = dat_gene_list$UMRwithoutCGI_intergenic, pAdjustMethod = "BH", qvalueCutoff = 0.01)
