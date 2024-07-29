library("AnnotationHub")
library("ensembldb")

release <- 101

anno <- query(AnnotationHub(), pattern=c("Mus musculus", "EnsDb", release))[[1]]

proms <- promoters(genes(anno), upstream=1000, downstream=100)

library(GenomicRanges)

# Access the metadata columns using mcols()
metadata <- mcols(proms)

# Create a logical vector where 'entrezid' is not NA
non_na_indices <- !is.na(metadata$entrezid)

# Subset 'proms' using this logical vector
filtered_proms <- proms[non_na_indices]

promoters_filtered = as.data.frame(filtered_proms)
output.M = data.frame("chr"=promoters_filtered$seqnames, "s"=promoters_filtered$start, "e"=promoters_filtered$end, "gene"=promoters_filtered$gene_name)

write.table(output.M,"mm10.promoter.bed", quote=F, row.names = F, col.names = F)

#proms <- trim(proms)[, "gene_id"]
#proms
