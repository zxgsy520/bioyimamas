#!/usr/bin/env Rscript

library(igraph)
library(RColorBrewer)
options(bitmapType='cairo') #关闭服务器与界面的互动响应


correlation_network <- function(data, prefix){

    data <- read.table(data, header=T, sep="\t")
    colnames(data) <- c("Source", "Target", "Strand", "value")

    dsource <- unique(data$Source)
    dtarget <- unique(data$Target)
    dnode <- c()
    for (i in dsource){
        dnode <- c(dnode, "pink")
    }
    for (i in dtarget){
        dnode <- c(dnode, "green")
    }
    node <- data.frame(name=c(dsource, dtarget), carac=dnode)
    network <- graph_from_data_frame(d=data, vertices=node, directed=F)

    node_color = c()
    for (i in 1:length(V(network)$carac)){
        if (V(network)$carac[i] == "pink"){
            node_color[i] = "pink"
        }
        if (V(network)$carac[i] == "green"){
            node_color[i] = "green"
        }
    }

    link_color = c()
    for (i in 1:length(E(network)$Strand)){
        if (E(network)$Strand[i] == "+"){
            link_color[i] = "red"
        }
        if (E(network)$Strand[i] == "-"){
            link_color[i] = "blue"
        }
    }

    pdf(paste(prefix, ".network.pdf", sep=""))
    a <- dev.cur()
    png(paste(prefix, ".network.png", sep=""), width=1600, height=1600, res=72*3)
    dev.control('enable')
    plot(network, layout=layout.fruchterman.reingold, vertex.color=node_color, vertex.size=5,
        vertex.label.cex=0.5, vertex.label.color="black", edge.color=link_color)
    dev.copy(which=a)#复制png的图形给pdf
    dev.off()
    dev.off()

}


add_help_args <- function(args){

    if(length(args) != 2) {
        cat("Author:XuChang\n")
        cat("Email:xc@bioyigene.com\n")
        cat("Function:Plot a network.\n")
        cat("Example:correlation_network.R link.tsv prefix\n")
        cat("Input file format:\n")
        cat("Source\tTarget\tStrand\tvalue\n")
        cat("MW0104052\tRhodobacter\t-\t0.870586")
        quit()
    }

}


args <- commandArgs(trailingOnly=TRUE)
add_help_args(args)
correlation_network(args[1], args[2])
