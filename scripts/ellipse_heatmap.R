#!/usr/bin/env Rscript

library(psych)
library(corrplot)
library(reshape2)
options(bitmapType='cairo') #关闭服务器与界面的互动响应


run_data <- function(data){

    data <- read.table(data, header=T, row.names=1, sep="\t", comment.char="",quote = "")
    data <- data[which(rowSums(data) > 0),] #过滤所有丰度都为0得行
    data <- t(data)

    return(data)
}


ellipse_heatmap <- function(data1, data2, prefix){

    data1 <- run_data(data1)
    data2 <- run_data(data2)

    result <- corr.test(data1, data2, method="pearson", adjust="none")
    rvalue <- result$r
    rvalue <- t(rvalue)
    #print(rvalue)

    rows <- nrow(rvalue)
    columns <- length(rvalue[1,])
    if(rows<=20) {
        heights <- 6
    }else if(rows<=50) {
        heights <- 14
    }else {
        heights <- 18
    }

    if(columns<=4) {
        widths <- 6
    }else if(columns<=8) {
        widths <- 10
    }else {
        widths <- 14
    }

    pdf(paste(prefix, ".ellipse_heatmap.pdf", sep=""), width=widths, height=heights)
    a <- dev.cur()
    png(paste(prefix, ".ellipse_heatmap.png", sep=""), width=widths*60, height=650)
    dev.control("enable")
    corrplot(rvalue, method="ellipse")
    dev.copy(which=a)
    dev.off()
    dev.off()
    
}


add_help_args <- function(args){

    if(length(args) != 3) {
        cat("Version: v1.0.0\n")
        cat("Author:Xingguo Zhang\n")
        cat("Email:invicoun@foxmail.com\n")
        cat("Function:Mapping a heatmap of microbial metabolism correlations.\n")
        cat("Example:ellipse_heatmap.R msms_Intensity.txt diff_taxa_table_L5.txt  prefix\n")
        cat("Input file format:\n")
        cat("Tax Id\tsample1\tsample2\tsample3\n")
        cat("OTU1\t10\t9\t9")
        quit()
    }
}


args <- commandArgs(trailingOnly=TRUE)
add_help_args(args)
ellipse_heatmap(args[1], args[2], args[3])
