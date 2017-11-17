library(cowplot)
library(readr)
library(dplyr)
library(ggplot2)
library(tidyr)
library(magrittr)
library(AnnotationDbi)
library(annotate)
library(GO.db)
library(tmod)


#parameters to set:
species <- "human" #human or mouse
maxGOgroupSize <- 200
minGOgroupSize <- 10
methodToUse <- "hyper" #hyper (two lists) or AUC (sorted list)

print(paste("GO Database date:", tbl_df(GO_dbInfo()) %>% filter(name == "GOSOURCEDATE") %>% .$value))

if (species == "human") {
  library(org.Hs.eg.db) #switch to org.Mm.eg.db for mouse
  goSource <- 'org.Hs.eg'
} else if (species == "mouse") {
  library(org.Mm.eg.db) 
  goSource <- 'org.Mm.eg'
}

if (methodToUse == "AUC") {
  #for demostration purposes - genes sorted alphbetically, human genes
  #this could be replaced by your genes that have been sorted by their p-value multipled by the sign of the direction of effect
  sortedGenes <- sort(unique(names(as.list(org.Hs.egALIAS2EG))))
}
if (methodToUse == "hyper") {
  #the hypergeometric method just compares the hit list or foreground genes against the remaining genes tested
  #again, this is just for demonstration purposes
  sortedGenes <- sort(unique(names(as.list(org.Hs.egALIAS2EG))))
  hitListGenes <- head(sortedGenes, n=5000)
}


#select GO groups
if (exists("geneSetsGO") && length(geneSetsGO$MODULES2GENES) > 1000 ) { #assume it's already loaded - needs a fix to see if the variable is declared
} else {
  if (species == "human") {
    go_object <- as.list(org.Hs.egGO2ALLEGS)
  } else if (species == "mouse") {
    go_object <- as.list(org.Mm.egGO2ALLEGS)
  }
  
  symbolsInGO <- getSYMBOL(unique(unlist(go_object)), data=goSource)
  
  #build GO sets for tmod -slow
  tmodNames <- data.frame()
  modules2genes <- list()
  goGroupName <- names(go_object)[1]
  showMethods(Term)
  
  goCount <- length(go_object)
  count <- 1
  for(goGroupName in names(go_object)) {
    if (count %% 1000 == 0) print(paste(count, "of", goCount))
    count <- count + 1
    
    goGroup <- go_object[goGroupName]
    geneIDs <- unique(unlist(goGroup, use.names=F))  #discard evidence codes
    genesymbols <- unique(getSYMBOL(geneIDs, data=goSource))
    
    genesymbols <- intersect(genesymbols, sortedGenes) #get size after intersecting with our full gene set
    if (!(length(genesymbols) >= minGOgroupSize & length(genesymbols) <= maxGOgroupSize)) next();
    
    modules2genes[goGroupName] <- list(genesymbols)
    
    tmodNames <- rbind(tmodNames, data.frame(ID=goGroupName, Title = Term(goGroupName)))
  }
  geneSetsGO <- makeTmod(modules = tmodNames, modules2genes = modules2genes)
}

#sorted gene list method
if (method == "AUC") {
  result <- tbl_df(tmodUtest(c(sortedGenes), mset=geneSetsGO, qval = 1, filter = T))
  result %<>% rowwise() %>% mutate(P.Value = P.Value * 2) %>% ungroup() %>% mutate(adj.P.Val=p.adjust(P.Value)) #tmod runs one-sided tests
  result %<>% rowwise() %>% mutate(aspect = Ontology(ID)) #add the source ontology (could be filterd for just biological process)
  
  #collapse genesets that have the exact same set of genes
  (result %<>% group_by(U, N1, AUC, P.Value,adj.P.Val) %>% summarize(MainTitle = first(Title),  ID=paste(ID, collapse=","), aspect= first(aspect), allNames = if_else(n() > 1, paste(Title[2:length(Title)], collapse=","), "")))
  result %<>% arrange(P.Value)
  result$rank <- 1:nrow(result)
  result %<>% dplyr::select(MainTitle, N1, AUC, P.Value, adj.P.Val, everything()) %>% arrange(adj.P.Val)

  #print top 20 at top and bottom
  print(head(filter(result, AUC > 0.5) %>% dplyr::select(-ID), n=20))
  print(head(filter(result, AUC < 0.5) %>% dplyr::select(-ID), n=20))

  source("./ROCPlots.R")
  
  #plot the top and bottom GO Groups
  plots <- createPlots(sortedGenes, c("GO:0043492", "GO:0015293"), geneSetsGO)
  plots$AUCPlot
  plots$rasterPlot
  (bothPlots <- plot_grid(plots$AUCPlot, plots$rasterPlot,  nrow = 2, align = "v", rel_heights=c(1,0.2),scale = 0.95))
}

#this just compares two lists
if (method == "hyper") {
  result <- tbl_df(tmodHGtest(fg = hitListGenes, bg = sortedGenes, mset=geneSetsGO, qval = 1.01, filter = T))
  
  result %<>% rowwise() %>% mutate(aspect = Ontology(ID)) #add the source ontology (could be filterd for just biological process)
  result$rank <- 1:nrow(result)
  result %<>% dplyr::select(Title, overlap = b, setSize=B, hitListSize = n, genesWithGOAnnotation = N, P.Value, adj.P.Val, everything()) %>% arrange(adj.P.Val)
  
  #print top 20 
  print(head(result %>% dplyr::select(-ID, -E, -rank), n=20))
}


