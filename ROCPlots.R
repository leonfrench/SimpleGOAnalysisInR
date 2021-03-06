library(dplyr)
library(ggplot2)
library(tidyr)
library(plotROC)

#input 
#sorted genes, gene ID names, tmod object for getting genes <-> group
createPlots <- function(sortedGenes, groupIDs, tmodSets, customNames=NULL) {
  
  ranking <- tbl_df(data.frame(gene_symbol = rev(sortedGenes), rank = 1: length(sortedGenes), stringsAsFactors = F))
  geneToClass <- NULL
  for(groupID in groupIDs) {
    geneToClass <- bind_rows(geneToClass, tbl_df(data.frame(gene_symbol = unlist(tmodSets$MODULES2GENES[groupID]), group = tmodSets$MODULES[groupID,]$Title, stringsAsFactors = F)))  
  }
  geneToClassAUC <- left_join(ranking, geneToClass, by = "gene_symbol") %>% spread(key=group, value=group)
  geneToClassAUC %<>% gather(group, present, -gene_symbol, -rank) %>% filter(group != "<NA>")
  geneToClassAUC <- geneToClassAUC %>% mutate(present = if_else(is.na(present), 0, 1))
  
  #order groups by direction of AUC
  forOrder <- tmodUtest(c(sortedGenes), mset=tmodSets, qval = 1, filter = T)
  forOrder <- tbl_df(forOrder)
  forOrder %<>% filter(ID %in% groupIDs ) %>% arrange(desc(AUC))
  #forOrder %<>% filter(ID %in% groupIDs ) %>% arrange(desc(sign(AUC-0.5)* P.Value))
  geneToClassAUC$group <- factor(geneToClassAUC$group, levels= as.character(forOrder$Title))
  
  (AUCPlot <- ggplot(geneToClassAUC, aes(d = present, m = rank, color=group)) + ylab("") + 
    geom_roc(n.cuts=0) + 
    style_roc() + coord_cartesian(expand=F) +
    theme(legend.position = c(1,0), legend.justification = c(1, 0), legend.background= element_rect(fill = "transparent", colour = "transparent"), plot.margin=unit(c(.5,.5,.5,.5),"cm")) + 
    labs(color='Gene Group')  + 
    theme(strip.background = element_blank(), strip.placement = "inside", strip.text = element_blank()) )

  forOrder$labelWithAUC <- paste0(tmodSets$MODULES[groupID,]$Title, " (AUC=", signif(forOrder[groupID, "AUC"],digits=2), ")")
  forOrder %<>% mutate(labelWithAUC = paste0(Title, " (AUC=", signif(AUC,digits=2), ")")) %>% dplyr::select(group = Title, labelWithAUC)
  forOrder$group <- as.character(forOrder$group)
  geneToClassAUC$group<- as.character(geneToClassAUC$group)
  geneToClassAUC <- inner_join(geneToClassAUC, forOrder, by="group") %>% dplyr::select(-group) %>% dplyr::rename(group = labelWithAUC)

  geneToClassAUC$group <- factor(geneToClassAUC$group, levels= as.character(forOrder$labelWithAUC))

  geneToClassAUC$rank <- -1*geneToClassAUC$rank
  (rasterPlot <- ggplot(geneToClassAUC, aes(x = rank, y = present, color= group)) + 
    geom_blank() + 
    geom_vline(data = filter(geneToClassAUC, present == 1), aes(xintercept=rank, color=group)) + #,color="black") + #, size=0.07) + 
    theme_bw()+coord_cartesian(expand=F) +
    ylab("Transcriptomic cell type") + 
    facet_wrap(~group, strip.position="top",ncol=1) + #, switch = "both"
    theme(strip.background = element_blank(), strip.placement = "inside") + #, strip.text.y = element_text(angle = 180)) +
    theme(axis.title.y = element_blank(),  axis.text.y=element_blank(), axis.ticks.y=element_blank(),axis.ticks.x=element_blank()) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + guides(color=FALSE) +
    scale_x_continuous(name = paste0("Gene ranking (",length(unique(geneToClassAUC$gene_symbol))," genes)"), breaks= c(min(geneToClassAUC$rank)+700, max(geneToClassAUC$rank)-700), labels = c("Top Rank", "Bottom Rank")))
  returnPlots = list()
  returnPlots[["AUCPlot"]] <- AUCPlot
  returnPlots[["rasterPlot"]] <- rasterPlot
  returnPlots
}

