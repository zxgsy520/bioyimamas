#!/usr/bin/env Rscript

library(vegan)
library(psych)
library(ggplot2)
library(ggforce)
options(bitmapType='cairo') #关闭服务器与界面的互动响应


read_data <- function(data){

    data <- read.table(data, header=T, row.names=1, sep="\t", comment.char="", quote="")
    data <- data[which(rowSums(data) > 0),] #过滤所有丰度都为0得行
    data <- t(data)

    return(data)
}


filter_pvalue <- function(pvalue, rvalue){

    indexs <- which(apply(pvalue, 1, min, na.rm=TRUE) < 0.05)
    pvalue <- pvalue[indexs,]
    rvalue <- rvalue[indexs,]
    indexs <- which(apply(rvalue, 1, max, na.rm=TRUE) >= 0.8 | apply(rvalue, 1, min, na.rm=TRUE) <= -0.8)
    pvalue <- pvalue[indexs,]
    rvalue <- rvalue[indexs,]
    r <- list(pvalue=pvalue, rvalue=rvalue)

    return(r)
}


cca_plot<-function(data1, data2, prefix){

    data1 <- read_data(data1)
    data2 <- read_data(data2)

    result <- corr.test(data1, data2, method="spearman", adjust="none")
    rvalue <- result$r
    pvalue <- result$p
    rvalue <- t(rvalue)
    pvalue <- t(pvalue)

    r <- filter_pvalue(pvalue, rvalue)
    rvalue <- r$r
    pvalue <- r$p
    a <- c(rownames(rvalue))
    data2 <- t(data2)
    select_data <- data2[a,]
    select_data <- t(select_data)
    rownames(select_data) <- rownames(data1)
    #envdat <- envdat_raw[match(row.names(select_data), row.names(data1)),]
    decorana(select_data)
    data1 <- data.frame(data1)
    otu_cca <- cca(select_data ~ ., data1)
    envdat_raw <- data1
    res <- otu_cca
    xxxx <- summary(res)
    aa <- xxxx$concont$importance
    aa <- round(aa, 4)
    pdat <- res$CCA
    samples <- data.frame(sample=row.names(pdat$u), CCA1=pdat$u[,1], CCA2=pdat$u[,2])
    species <- data.frame(spece=row.names(pdat$v), CCA1=pdat$v[,1], CCA2=pdat$v[,2])
    envi <- data.frame(en=row.names(pdat$biplot), CCA1=pdat$biplot[,1], CCA2=pdat$biplot[,2])
    r <- max(max(species$CCA1),max(species$CCA2))#规定半径

    p <- ggplot() +
        geom_hline(aes(yintercept = 0)) +
        geom_vline(aes(xintercept = 0))  +
        geom_circle(aes(x0=0,y0=0,r=r*1.2)) +
        geom_circle(aes(x0=0,y0=0,r=r*0.8))+
        geom_point(data=species, aes(x=CCA1, y=CCA2, colour='Metabolites'), size=1.5, show.legend=T) +
        geom_point(data=envi, aes(x=CCA1, y=CCA2, colour='Microbe'), size=1.5, show.legend=T)  +
        theme_bw() + theme(panel.grid.major=element_line(colour=NA), panel.grid.minor=element_blank()) +
        xlab(paste('CCA1 (', aa[2,1]*100, '%)', sep ='')) + ylab(paste('CCA2 (', aa[2,2]*100, '%)', sep =''))+
        theme(legend.title=element_blank())
        ggsave(paste(prefix, ".plot_cca.pdf", sep=""), p, width=135, height=100, units="mm")
        ggsave(paste(prefix, ".plot_cca.png", sep=""), p, width=135, height=100, units="mm")
}


add_help_args <- function(args){

  if(length(args) != 3) {
    cat("Version: v1.0.0\n")
    cat("Author:BoYaXu\n")
    cat("Email:invicoun@foxmail.com\n")
    cat("Function:CCA Analysis.\n")
    cat("Example:cca_analysis.R genus.tsv differential_metabolites.tsv prefix\n")
    quit()
  }
}


args <- commandArgs(trailingOnly=TRUE)
add_help_args(args)
cca_plot(args[1], args[2], args[3])
