####################READING####################
nodes <- read.table(IBD_nodes_file, sep = "\t", h = T)[,c(1:4)]
colnames(nodes) <- c("ID", "individual", "date", "location")

edges <- read.table(IBD_edges_file, sep = "\t", h = T)[,c(1:3)]
colnames(edges) <- c("ID1", "ID2", "fract_sites_IBD")

ibdsegments <- read_delim(file = hmmIBD_file, delim = "\t")[,c(1:7)]
colnames(ibdsegments) <- c("ID1", "ID2", "chr", "start", "end", "different", "Nsnp")
########################################

#############FIND PARENTS#############
edges.filtered <- edges[edges$fract_sites_IBD>=parentalmin & edges$fract_sites_IBD<=parentalmax,]
# Get rid of loops and ensure right naming of vertices
edges.topology <- edges.filtered[,c("ID1", "ID2")]
edges.topology <- simplify(graph_from_data_frame(edges.topology[order(edges.topology[[1]]),],directed = FALSE))

# Find all components
parent <- components(edges.topology)
parent.df <- data.frame(parent$membership)
parent.df <- cbind(row.names(parent.df), parent.df)
colnames(parent.df) <- c("ID", "parent")
nodes <- merge(nodes, parent.df, all.x = T)
nodes$parent[which(is.na(nodes$parent))] <- c(1:length(which(is.na(nodes$parent))))+
  max(nodes$parent, na.rm = TRUE)

nodes$parent <- sprintf(paste0("%0",max(nchar(nodes$parent)),"d"), nodes$parent)
##########################

#################FIND IDENTICAL SAMPLES#################
edges.filtered <- edges[edges$fract_sites_IBD>=identical,]
# Get rid of loops and ensure right naming of vertices
edges.topology <- edges.filtered[,c("ID1", "ID2")]
edges.topology <- simplify(graph_from_data_frame(edges.topology[order(edges.topology[[1]]),],directed = FALSE))

# Find all components
strain <- components(edges.topology)
strain.df <- data.frame(strain$membership)
strain.df <- cbind(row.names(strain.df), strain.df)
colnames(strain.df) <- c("ID", "strain")
nodes <- merge(nodes, strain.df, all.x = T)
nodes$strain[which(is.na(nodes$strain))] <- c(1:length(which(is.na(nodes$strain))))+
  max(nodes$strain, na.rm = TRUE)

nodes$strain <- sprintf(paste0("%0",max(nchar(nodes$strain)),"d"), nodes$strain)
##########################

####################FIND TRIADES####################
edges <- merge(edges, nodes, by.x = "ID1", by.y = "ID",
               suffixes = c("1", "2"))
edges <- merge(edges, nodes, by.x = "ID2", by.y = "ID",
               suffixes = c("1", "2"))

siblings <- edges[edges$parent1 == edges$parent2 & !edges$strain1 == edges$strain2 & 
                    edges$fract_sites_IBD >= 0,]
families <- siblings%>%
  group_by(parent1)%>%
  summarise(IBDmax = max(fract_sites_IBD), IBDmin = min(fract_sites_IBD))

allchildren <- data.frame("ID" = NA, "parent1" = NA, "parent2" = NA, 
                          "parentsIBD" = NA, "misbases" = NA,
                          "parent1IBD" = NA, "parent2IBD" = NA, 
                          "parent1bases" = NA, "parent2bases" = NA, 
                          "samebases" = NA, "family" = NA)

for (family in families$parent1){
  siblings.family <- siblings[siblings$parent1==family & siblings$fract_sites_IBD < parentalmax,]
  siblings.family.unrelated <- siblings.family[siblings.family$fract_sites_IBD < different,]
  siblings.family.unrelated[,c("ID1", "ID2")] <- data.frame(pmin(siblings.family.unrelated$ID1, 
                                                              siblings.family.unrelated$ID2),
                                                         pmax(siblings.family.unrelated$ID1,
                                                              siblings.family.unrelated$ID2)) 
  
  siblings.family.related <- siblings.family[siblings.family$fract_sites_IBD > parentalmin,]
  
  wayone <- siblings.family.related%>%
    group_by(ID1)%>%
    reframe(expand.grid(ID2, ID2, stringsAsFactors = F))
  waytwo <- siblings.family.related%>%
    group_by(ID2)%>%
    reframe(expand.grid(ID1, ID1, stringsAsFactors = F))
  
  colnames(wayone) <- c("ID", "parent1", "parent2")
  colnames(waytwo) <- c("ID", "parent1", "parent2")
  
  children <- rbind(wayone, waytwo)
  
  children[,c("parent1", "parent2")] <- data.frame(pmin(children$parent1, 
                                                        children$parent2),
                                                   pmax(children$parent1,
                                                        children$parent2))
  children$parentsIBD <- apply(children, 1,
                        FUN = function(x){
                          value <- siblings.family.unrelated$fract_sites_IBD[siblings.family.unrelated$ID1 == x[2] & siblings.family.unrelated$ID2 == x[3]]
                          if (length(value) > 0){
                            return(value)
                          }else{
                            return(NA)
                          }
                        })
  children <- children[!is.na(children$parentsIBD),]
  children <- children[!duplicated(children),]
  if(nrow(children)==0){
    next
  }
  children$misbases <- apply(children,1,
                                FUN = function(x){
                                  dfunkown <- parentsdiff(ibdsegments, x[1], x[2], x[3])
                                  return(sum(dfunkown$end - dfunkown$start))
                                })
  children <- children[order(children$ID, children$misbases),]
  children.best <- children[!duplicated(children$ID),]
  children.best[order(children.best$parent1, children.best$parent2),]
  children.best$parent1IBD <- apply(children.best, 1,
                               FUN = function(x){
                                 value <- siblings.family$fract_sites_IBD[siblings.family$ID1 == x[1] & siblings.family$ID2 == x[2]]
                                 value2 <- siblings.family$fract_sites_IBD[siblings.family$ID1 == x[2] & siblings.family$ID2 == x[1]]
                                 if (length(value) == 1){
                                   return(value)
                                 }else{
                                   return(value2)
                                 }
                               })
  children.best$parent2IBD <- apply(children.best, 1,
                               FUN = function(x){
                                 value <- siblings.family$fract_sites_IBD[siblings.family$ID1 == x[1] & siblings.family$ID2 == x[3]]
                                 value2 <- siblings.family$fract_sites_IBD[siblings.family$ID1 == x[3] & siblings.family$ID2 == x[1]]
                                 if (length(value) == 1){
                                   return(value)
                                 }else{
                                   return(value2)
                                 }
                               })
  children.best$parent1bases <- apply(children.best,1,
                                      FUN = function(x){
                                        subdata <- subduo(ibdsegments, x[1], x[2])
                                        return(sum(subdata$end[subdata$different==0] - subdata$start[subdata$different==0]))
                                      })
  children.best$parent2bases <- apply(children.best,1,
                                      FUN = function(x){
                                        subdata <- subduo(ibdsegments, x[1], x[3])
                                        return(sum(subdata$end[subdata$different==0] - subdata$start[subdata$different==0]))
                                      })
  children.best$samebases <- apply(children.best,1,
                               FUN = function(x){
                                 dfsame <- parentssame(ibdsegments, x[1], x[2], x[3])
                                 return(sum(dfsame$end - dfsame$start))
                               })
  children.best$family <- family
  allchildren <- rbind(allchildren, children.best)
}

allchildren <- allchildren[!is.na(allchildren$ID),]
########################################

########################SAVE DATA###################
write.table(allchildren, paste0(outDir, "/triades.tsv"),
            sep ="\t", quote = F, row.names = F)
########################################

