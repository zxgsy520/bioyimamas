#!/usr/bin/env Rscript

library(psych)
library(vegan)
library(reshape2)
library(ComplexHeatmap)
library(circlize)
options(bitmapType='cairo') #关闭服务器与界面的互动响应


read_data <- function(data){
  
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
}#筛选标注哪些小于0.01 小于0.05并标注


filter_pvalue <- function(pvalue, rvalue){
  
    indexs <- which(apply(pvalue, 1, min, na.rm=TRUE) < 0.05)
    pvalue <- pvalue[indexs,]
    rvalue <- rvalue[indexs,]
  
    r <- list(pvalue=pvalue, rvalue=rvalue)
  
    return(r)
}


correlation_analysis <- function(data1, data2, prefix){
  
    data1 <- read_data(data1)
    data2 <- read_data(data2)
  
    result <- corr.test(data1, data2, method="pearson", adjust="none") #计算相关性
    rvalue <- result$r
    pvalue <- result$p
    rvalue <- t(rvalue)
    pvalue <- t(pvalue)
  
    r <- filter_pvalue(pvalue, rvalue)
    rvalue <- r$rvalue
    pvalue <- r$pvalue
  
    result <- melt(rvalue, value.name="cor")
    result$pvalue <- as.vector(pvalue)
    write.table(result, file=paste(prefix, ".correlation_pvalue.csv", sep=""), sep= ",")
    data.1<-scale(data1)
    data.2<-scale(data2)
    rvalue<-t(rvalue)
    pvalue <- mark_pvalue(pvalue)
    pvalue<-t(pvalue)
    data.2 <- data.2[,match(colnames(rvalue), colnames(data.2))]
    data.1<-t(data.1)
    #p3 <- Heatmap(data.1,cluster_columns = F, cluster_rows = F)
    col_fun = colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))
    #ha <- rowAnnotation(empty = anno_empty(border = FALSE),Microbiome = data.1,annotation_legend_param = list(Microbiome = list(title = "Microbiome")), col = list(Microbiome = col_fun))
    ha <- rowAnnotation(empty=anno_empty(border=FALSE), 
        Microbiome=data.1, col=list(Microbiome=col_fun), 
        show_legend=F, annotation_label="Microbiome")
    #right_annotation=ha
    p1 <- Heatmap(rvalue, name="Correlation", 
        column_title="",
        column_title_gp=gpar(fontface='bold', fontsize=20),
        cluster_columns=T, cluster_rows=T, show_row_names=F,
        show_heatmap_legend=TRUE, width=unit(90, "cm"), height=unit(5, "cm"),
        column_dend_height=unit(3,"cm"),
        cell_fun = function(j, i, x, y, width, height, fill) {
            grid.text(sprintf("%s",pvalue[i, j]), x, y, gp = gpar(fontsize = 10))
        }, right_annotation=ha)

    p2 <- Heatmap(data.2, name="Metabolome", column_title="Metabolome",
        column_title_side="bottom", cluster_columns=F,
        cluster_rows=F, width=unit(90, "cm"), height=unit(5, "cm"))

    lgd.list = Legend(col_fun=col_fun, title="Microbiome")
    ht = p1%v%p2
    pdf(paste(prefix, ".combine_heatmap.pdf", sep=""), width=45, height=10)
    ptemp <- dev.cur()
    png(paste(prefix, ".combine_heatmap.png", sep=""), width=4000, height=900)
    dev.control("enable")

    draw(ht, annotation_legend_list=lgd.list, column_title="Metabolome",
        column_title_side="bottom", column_title_gp=gpar(fontface='bold', fontsize=20),
        padding=unit(c(2, 2, 10, 2), "mm"))
    decorate_annotation("Microbiome", { 
        grid.text("Microbiome", gp=gpar(fontface='bold', fontsize=20) ,
        y=unit(1, "npc") + unit(2, "mm"), just = "bottom")}
    )

    dev.copy(which=ptemp)
    dev.off()
    dev.off()
}


add_help_args <- function(args){
  
  if(length(args) != 3) {
    cat("Version: v1.0.0\n")
    cat("Author:BoYaXu\n")
    cat("Email:invicoun@foxmail.com\n")
    cat("Function:Plotting a heatmap of metabolic and microbial correlations.\n")
    cat("Example:combine_heatmap.R genus.tsv differential_metabolites.tsv prefix\n")
    quit()
  }
}


args <- commandArgs(trailingOnly=TRUE)
add_help_args(args)
correlation_analysis(args[1], args[2], args[3])
