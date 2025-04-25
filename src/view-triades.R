########################PLOT###################

for (i in 1:nrow(allchildren)){
  child <- as.character(allchildren[i,1])
  parent1 <- as.character(allchildren[i,2])
  parent2 <- as.character(allchildren[i,3])
  offspringDir <- paste0(outDir, "/offsprings/", nodes$parent[nodes$ID==child])
  if (!dir.exists(offspringDir)){
    dir.create(offspringDir, recursive = TRUE)
  }
  siblingsDir <- paste0(outDir, "/siblings/")
  if (!dir.exists(siblingsDir)){
    dir.create(siblingsDir, recursive = TRUE)
  }
  
  gg <- parentsmap(chrombarcodes, 
                   child, 
                   parent1, parent2,
                   chlen = chromlength)
  png(paste0(offspringDir, "/",
             child, "_from_", parent2, "-", parent1, ".png"), 
      width = 8, height = 5, res = 200, units = "in")
  plot(gg)
  dev.off()
  gg <- parentsmap(chrombarcodes, 
                   child, 
                   parent2, parent1,
                   chlen = chromlength)
  png(paste0(offspringDir, "/",
             child, "_from_", parent1, "-", parent2, ".png"), 
      width = 8, height = 5, res = 200, units = "in")
  plot(gg)
  dev.off()
  allnextchildren <- allchildren[-c(1:i),1]
  next
  for (nextchild in allnextchildren){
    gg <- parentsmap(chrombarcodes, 
                     child, 
                     nextchild, nextchild,
                     chlen = chromlength, mode = 1)
    gg
    ggsave(paste0(siblingsDir, "/", child, "_with_", nextchild, ".png"), 
           plot = gg, width = 8, height = 5)
  }
}

########################################
