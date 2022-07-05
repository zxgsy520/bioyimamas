#!/usr/bin/env Rscript

library(vegan)
library(ggrepel) 
library(ggplot2)
options(bitmapType='cairo') #关闭服务器与界面的互动响应


metabolism_cca <- function(data1, data2, group, prefix){

    data1 <- read.table(data1, header=T, row.names=1, sep="\t", comment.char="", quote="")
    data1 <- data.frame(t(data1))
    data2 <- read.table(data2, header=T, row.names=1, sep="\t", comment.char="", quote="")
    data2 <- data.frame(t(data2))
    
    otu_cca <- cca(data1~., data2)
    cca_coef <- coef(otu_cca)
    r2 <- RsquareAdj(otu_cca)
    otu_cca_noadj <- r2$r.squared   #原始 R2
    otu_cca_adj <- r2$adj.r.squared #校正后的 R2
    #关于约束轴承载的特征值或解释率，应当在 R2 校正后重新计算
    otu_cca_exp_adj <- otu_cca_adj * otu_cca$CCA$eig/sum(otu_cca$CCA$eig)
    otu_cca_eig_adj <- otu_cca_exp_adj * otu_cca$tot.chi

    #所有约束轴的置换检验，即全局检验，基于 999 次置换，详情 ?anova.cca
    otu_cca_test <- anova.cca(otu_cca, permutations=999)
    #各约束轴逐一检验，基于 999 次置换
    otu_cca_test_axis <- anova.cca(otu_cca, by ='axis', permutations=999)

    #p 值校正（Bonferroni 为例）
    otu_cca_test_axis$`Pr(>F)` <- p.adjust(otu_cca_test_axis$`Pr(>F)`, method = 'bonferroni')

    #计算方差膨胀因子，详情 ?vif.cca
    vif.cca(otu_cca)

    #前向选择，这里以 ordiR2step() 的方法为例，基于 999 次置换检验，详情 ?ordiR2step
    otu_cca_forward_pr <- ordiR2step(cca(data1~1, data2), scope=formula(otu_cca), R2scope=TRUE, direction='forward', permutations=999)
    #计算校正 R2 后的约束轴解释率
    exp_adj <- RsquareAdj(otu_cca_forward_pr)$adj.r.squared * otu_cca_forward_pr$CCA$eig/sum(otu_cca_forward_pr$CCA$eig)
    cca1_exp <- paste('CCA1:', round(exp_adj[1]*100, 2), '%')
    cca2_exp <- paste('CCA2:', round(exp_adj[2]*100, 2), '%')

    #下面是 ggplot2 方法
    #提取样方和环境因子排序坐标，前两轴，I 型标尺
    otu_cca_forward_pr.scaling1 <- summary(otu_cca_forward_pr, scaling = 1)
    otu_cca_forward_pr.site <- data.frame(otu_cca_forward_pr.scaling1$sites)[1:2]
    otu_cca_forward_pr.env <- data.frame(otu_cca_forward_pr.scaling1$biplot)[1:2]

    #添加分组
    otu_cca_forward_pr.env$name <- rownames(otu_cca_forward_pr.env)
    #读取分组文件按
    map <- read.delim(group, row.names = 1, sep = '\t')
    colnames(map) <- c("group")
    otu_cca_forward_pr.site$name <- rownames(otu_cca_forward_pr.site)
    otu_cca_forward_pr.site$group <- map$group
    #merged2<-merge(map,otu,by="row.names",all.x=TRUE)

    color=c( "#3C5488B2","#00A087B2",
           "#F39B7FB2","#91D1C2B2",
           "#8491B4B2", "#DC0000B2",
           "#7E6148B2","yellow",
           "darkolivegreen1", "lightskyblue",
           "darkgreen", "deeppink", "khaki2",
           "firebrick", "brown1", "darkorange1",
           "cyan1", "royalblue4", "darksalmon",
           "darkgoldenrod1", "darkseagreen", "darkorchid")

    p <- ggplot(otu_cca_forward_pr.site, aes(CCA1, CCA2)) +
    geom_point(size=1,aes(color = group,shape = group)) +
    stat_ellipse(aes(color = group), level = 0.95, show.legend = FALSE, linetype = 2) +
    scale_color_manual(values = color[1:length(unique(map$group))]) +
    theme(panel.grid.major = element_line(color = 'gray', size = 0.1), panel.background = element_rect(color = 'black', fill = 'transparent'),
      legend.title = element_blank(), legend.key = element_rect(fill = 'transparent'), plot.title = element_text(hjust = 0.5)) +
    labs(x = cca1_exp, y = cca2_exp) +
    geom_vline(xintercept = 0, color = 'gray', size = 0.5) +
    geom_hline(yintercept = 0, color = 'gray', size = 0.5) +
    geom_segment(data = otu_cca_forward_pr.env, aes(x = 0, y = 0, xend = CCA1, yend = CCA2), arrow = arrow(length = unit(0.2, 'cm')), size = 0.3, color = 'blue') +
    geom_text(data = otu_cca_forward_pr.env, aes(CCA1 * 1.2, CCA2 * 1.2, label = name), color = 'blue', size = 3) +
    geom_label_repel(aes(label = name, color = group), size = 3, box.padding = unit(0, 'lines'), show.legend = FALSE)

    ggsave(paste(prefix, ".cca.pdf", sep=""), p, width = 5.5, height = 5.5)
    ggsave(paste(prefix, ".cca.png", sep=""), p, width = 5.5, height = 5.5)
   

}


add_help_args <- function(args){

    if(length(args) != 4) {
        cat("Version: v1.0.0\n")
        cat("Author:Xingguo Zhang\n")
        cat("Email:invicoun@foxmail.com\n")
        cat("Function:Metabolism and Species Association Analysis.\n")
        cat("Example:metabolism_association.R msms_Intensity.txt diff_taxa_table_L5.txt  group.list prefix\n")
        cat("Input file format:\n")
        cat("Tax Id\tsample1\tsample2\tsample3\n")
        cat("OTU1\t10\t9\t9")
        quit()
    }
}


args <- commandArgs(trailingOnly=TRUE)
add_help_args(args)
metabolism_cca(args[1], args[2], args[3], args[4])
