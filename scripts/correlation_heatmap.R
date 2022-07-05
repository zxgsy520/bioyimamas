#!/usr/bin/env Rscript

library(psych)
library(pheatmap)
library(reshape2)
options(bitmapType='cairo') #关闭服务器与界面的互动响应


run_data <- function(data){

    data <- read.table(data, header=T, row.names=1, sep="\t", comment.char="",quote = "")
    data <- data[which(rowSums(data) > 0),] #过滤所有丰度都为0得行
    data <- t(data)

    return(data) 
}


mark_pvalue <- function(pvalue){

    if(!is.null(pvalue)){
        site1 <- pvalue<0.01
        pvalue[site1] <- "**"
        site2 <- pvalue > 0.01& pvalue < 0.05
        pvalue[site2] <- "*"
        pvalue[!site2&!site1]<- ""
    } else{
        pvalue <- F
    }
    return(pvalue)
}


filter_pvalue <- function(pvalue, rvalue){

    indexs <- which(apply(pvalue, 1, min, na.rm=TRUE) < 0.05)
    pvalue <- pvalue[indexs,]
    rvalue <- rvalue[indexs,]

    r <- list(pvalue=pvalue, rvalue=rvalue)

    return(r)
}


correlation_analysis <- function(data1, data2, prefix){

    data1 <- run_data(data1)
    data2 <- run_data(data2)

    result <- corr.test(data1, data2, method="pearson", adjust="none")
    rvalue <- result$r
    pvalue <- result$p
    rvalue <- t(rvalue)
    #print(rvalue)
    pvalue <- t(pvalue)

    r <- filter_pvalue(pvalue, rvalue)
    rvalue <- r$rvalue
    pvalue <- r$pvalue

    result <- melt(rvalue, value.name="cor")
    result$pvalue <- as.vector(pvalue)
    write.table(result, file=paste(prefix, ".correlation.csv", sep=""), sep= ",")

    pvalue <- mark_pvalue(pvalue)


    mycol <- colorRampPalette(c("blue","white","tomato"))(800)
    pheatmap(rvalue, scale="none", cluster_row=T, cluster_col=T,
             border=NA, display_numbers=pvalue, fontsize_number=12,
             number_color = "white", cellwidth=12, cellheight=15, 
             color=mycol, filename=paste(prefix, ".correlation_heatmap.pdf", sep=""))
    pheatmap(rvalue, scale="none", cluster_row=T, cluster_col=T,
             border=NA, display_numbers=pvalue, fontsize_number=12,
             number_color = "white", cellwidth=12, cellheight=15,
             color=mycol, filename=paste(prefix, ".correlation_heatmap.png", sep=""))
}


add_help_args <- function(args){

    if(length(args) != 3) {
        cat("Version: v1.0.1\n")
        cat("Author:Xingguo Zhang\n")
        cat("Email:invicoun@foxmail.com\n")
        cat("Function:Metabolism and Species Association Analysis.\n")
        cat("Example:correlation_heatmap.R msms_Intensity.txt diff_taxa_table_L5.txt  prefix\n")
        cat("Input file format:\n")
        cat("Tax Id\tsample1\tsample2\tsample3\n")
        cat("OTU1\t10\t9\t9")
        quit()
    }
}


args <- commandArgs(trailingOnly=TRUE)
add_help_args(args)
correlation_analysis(args[1], args[2], args[3])
