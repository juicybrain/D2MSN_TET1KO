    library(hictoolsr)
    library("dbscan")
        loop_25k_files = c("file1", "file2")
        mergeBedpe(bedpeFiles = loop_25k_files, res=25000, selectCol = 12, dist_method = "manhattan", minPts =2  )
        write.table( dat_loop_25k, "dat_loop_25k.bedpe", quote=F, sep="\t", row.names=F)
       
