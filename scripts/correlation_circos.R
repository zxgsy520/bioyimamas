#!/usr/bin/env Rscript

library(circlize)
options(bitmapType='cairo') #关闭服务器与界面的互动响应


correlation_circos <- function(data, prefix){

    data <- read.table(data, header=T, sep="\t")
    colnames(data) <- c("Source", "Target", "Strand", "value")

    link_color <- c()
    for (i in 1:length(data$Strand)){
    	  if (data$Strand[i] == "+"){
    		    link_color[i] = "pink"
        }
    	  if (data$Strand[i] == "-"){
            link_color[i] = "blue"
        }
    }

    data <- data.frame(data$Target, data$Source, data$value, link_color)
    pdf(paste(prefix, ".circos.pdf", sep=""))
    a <- dev.cur()
    png(paste(prefix, ".circos.png", sep=""), width=1700, height = 1700, res=72*4)
    circos.par(start.degree=360, clock.wise=FALSE)
    chordDiagram(data, annotationTrack="grid", preAllocateTracks=list(track.height=0.1), col=data$link_color)
                 #grid.col=link_color, transparency=0)
    circos.track(track.index=1, panel.fun=function(x, y){
                circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,
                    facing="clockwise",
                    niceFacing=TRUE,
                    adj=c(0, 0.5),
                    cex=0.5
                    )
                }, bg.border = NA)

    circos.clear()
    dev.copy(which=a)
    dev.off()
    dev.off()
}


add_help_args <- function(args){

    if(length(args) != 2) {
        cat("Author:XuChang\n")
        cat("Email:xc@bioyigene.com\n")
        cat("Function:Plot the correlation circle\n")
        cat("Example:Rscript metabolism_cca.R link.tsv prefix\n")
        cat("Input file format:\n")
        cat("Source\tTarget\tStrand\tvalue\n")
        cat("MW0104052\tRhodobacter\t-\t0.870586")
        quit()
    }

}

args <- commandArgs(trailingOnly=TRUE)
add_help_args(args)
correlation_circos(args[1], args[2])
