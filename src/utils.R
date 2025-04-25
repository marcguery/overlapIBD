
##################FUNCTIONS###############
sets <-function(start, end, group, overlap=length(unique(group))) {
  dd<-rbind(data.frame(pos=start, event=1), data.frame(pos=end, event=-1))
  dd<-aggregate(event~pos, dd, sum)
  dd<-dd[order(dd$pos),]
  dd$open <- cumsum(dd$event)
  r<-rle(dd$open>=overlap)
  ex<-cumsum(r$lengths-1 + rep(1, length(r$lengths))) 
  sx<-ex-r$lengths+1
  data.frame(start = dd$pos[sx[r$values]],
             end = dd$pos[ex[r$values]+1])
} 
parentsdiff <- function(ibddata, child, parent1, parent2){
  subdata.1 <- subduo(ibddata, child, parent1)
  subdata.1$id <- 1
  subdata.1 <- subdata.1[subdata.1$different==1,]
  subdata.2 <- subduo(ibddata, child, parent2)
  subdata.2$id <- 2
  subdata.2 <- subdata.2[subdata.2$different==1,]
  subdata <- rbind(subdata.1, subdata.2)
  if(nrow(subdata)==0){
    return(subdata[,c("chr", "start", "end")])
  }
  subdata%>%
    group_by(chr)%>%
    do(sets(.$start, .$end, .$id, 2))
}
parentssame <- function(ibddata, child, parent1, parent2){
  subdata.1 <- subduo(ibddata, child, parent1)
  subdata.1$id <- 1
  subdata.1 <- subdata.1[subdata.1$different==0,]
  subdata.2 <- subduo(ibddata, child, parent2)
  subdata.2$id <- 2
  subdata.2 <- subdata.2[subdata.2$different==0,]
  subdata <- rbind(subdata.1, subdata.2)
  if(nrow(subdata)==0){
    return(subdata[,c("chr", "start", "end")])
  }
  subdata%>%
    group_by(chr)%>%
    do(sets(.$start, .$end, .$id, 2))
}

subduo <- function(ibddata, above, below){
  subdata <- ibddata[ibddata$ID1==above,]
  subdata <- subdata[subdata$ID2==below,]
  if (dim(subdata)[1]==0){
    subdata <- ibddata[ibddata$ID2==above,]
    subdata <- subdata[subdata$ID1==below,]
  }
  if(dim(subdata)[1]==0){
    stop("There is no correlation between ", above, " and ", below)
  }
  return(subdata)
  
}

parentsmap <- function(ibddata, child, parent1, parent2, chlen, mode=0){
  print(paste(child, parent1, parent2))
  subdata.1 <- subduo(ibddata, child, parent1)
  subdata.2 <- subduo(ibddata, child, parent2)
  bothdiff <- parentsdiff(ibddata, child, parent1, parent2)
  bothsame <- parentssame(ibddata, child, parent1, parent2)
  
  chs <- chlen
  if(mode==0){
    gg <- ggplot()+
      geom_rect(data=chs,
                aes(xmin = 0, xmax=length, 
                    ymin=as.numeric(chr)*2.75-0.5, ymax=as.numeric(chr)*2.75+0.5), 
                color=NA, fill="grey80")+
      geom_text(data=chs,
                aes(x = 0-100000, y=as.numeric(chr)*2.75, 
                    label=chr), 
                color="grey50", size=3)+
      geom_rect(data=subdata.2[subdata.2$different==0,],
                aes(xmin=start, xmax=end, 
                    ymin=as.numeric(chr)*2.75-0.5, ymax=as.numeric(chr)*2.75+0.5),
                color = NA, fill = "gold1", linewidth = 0.1)+
      geom_rect(data=subdata.1[subdata.1$different==0,],
                aes(xmin=start, xmax=end, 
                    ymin=as.numeric(chr)*2.75-0.5, ymax=as.numeric(chr)*2.75+0.5), 
                fill = "blue", color = NA, linewidth = 0.1)+
      geom_rect_pattern(data=bothsame,
                        aes(xmin=start, xmax=end, 
                            ymin=as.numeric(chr)*2.75-0.5, ymax=as.numeric(chr)*2.75+0.5), 
                        fill = "blue", pattern_fill ="gold1", pattern_density = 0.5, pattern_spacing = 0.015,
                        color = NA, pattern_color = NA, linewidth = 0.1)+
      geom_rect(data=bothdiff,
                aes(xmin=start, xmax=end, 
                    ymin=as.numeric(chr)*2.75-0.5, ymax=as.numeric(chr)*2.75+0.5),
                color = NA, fill = "pink", linewidth = 0.1)+
      labs(title = paste(child, "from", parent1, "(blue) and", parent2, "(yellow)"))+
      scale_x_continuous(breaks = seq(0,3.25*10^6,0.25*10^6), 
                         labels = seq(0,3.25*10^6,0.25*10^6)/10^6,
                         minor_breaks = seq(0,3.25*10^6,0.125*10^6),
                         expand = c(0.01,0.01))+
      scale_y_discrete(expand = c(0,0))+
      theme(plot.background = element_blank(),
            panel.grid.major.y = element_blank(),
            panel.grid.minor.y = element_blank(),
            panel.grid.major.x = element_line(color = "grey60", linetype = 2, linewidth = 0.5),
            panel.grid.minor.x = element_line(color = "grey80", linetype = 2, linewidth = 0.5),
            panel.background = element_rect(fill = "grey95"),
            legend.position = "none",
            axis.title = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.ticks.x = element_line(color = "grey60", linewidth = 0.5))
    return(gg)
  }
  
  if(mode==1){
    gg <- ggplot()+
      geom_rect(data=chs,
                aes(xmin = 0, xmax=length, 
                    ymin=as.numeric(chr)*2.75-0.5, ymax=as.numeric(chr)*2.75+0.5), 
                color=NA, fill="grey80")+
      geom_text(data=chs,
                aes(x = 0-100000, y=as.numeric(chr)*2.75, 
                    label=chr), 
                color="grey50", size=3)+
      geom_rect(data=subdata.1[subdata.1$different!=1,],
                aes(xmin=start, xmax=end, 
                    ymin=as.numeric(chr)*2.75-0.5, ymax=as.numeric(chr)*2.75+0.5),
                color = NA, fill = "lightgreen", linewidth = 0.1)+
      geom_rect(data=subdata.1[subdata.1$different==1,],
                aes(xmin=start, xmax=end, 
                    ymin=as.numeric(chr)*2.75-0.5, ymax=as.numeric(chr)*2.75+0.5),
                color = NA, fill = "pink", linewidth = 0.1)+
      labs(title = paste(child, "with", parent1))+
      scale_x_continuous(breaks = seq(0,3.25*10^6,0.25*10^6), 
                         labels = seq(0,3.25*10^6,0.25*10^6)/10^6,
                         minor_breaks = seq(0,3.25*10^6,0.125*10^6),
                         expand = c(0.01,0.01))+
      scale_y_discrete(expand = c(0,0))+
      theme(plot.background = element_blank(),
            panel.grid.major.y = element_blank(),
            panel.grid.minor.y = element_blank(),
            panel.grid.major.x = element_line(color = "grey60", linetype = 2, linewidth = 0.5),
            panel.grid.minor.x = element_line(color = "grey80", linetype = 2, linewidth = 0.5),
            panel.background = element_rect(fill = "grey95"),
            legend.position = "none",
            axis.title = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.ticks.x = element_line(color = "grey60", linewidth = 0.5))
    return(gg)
    
  }
}

familymap <- function(ibddata, genome, versus, chlen){
  maxother = min(8, length(versus))
  coordsplot <- c(-0.5,c(-0.5)+c(1:maxother)/maxother)
  colorlist <- brewer.pal(n = 8, name = "Set1")
  chs <- chlen
  
  gg <- ggplot()+
    geom_rect(data=chs,
              aes(xmin = 0, xmax=length, 
                  ymin=as.numeric(chr)*2.75-0.5, ymax=as.numeric(chr)*2.75+0.5), 
              color=NA, fill="grey80")+
    geom_text(data=chs,
              aes(x = 0-100000, y=as.numeric(chr)*2.75, 
                  label=chr), 
              color="grey50", size=3)
  
  for (i in c(1:maxother)){
    othergenome <- versus[i]
    if (othergenome == genome){
      selfcomp <- chs
      selfcomp$coord1 <- coordsplot[i]
      selfcomp$coord2 <- coordsplot[i+1]
      gg <- gg+
        geom_rect(data=selfcomp,
                  aes(xmin=0, xmax=length, 
                      ymin=as.numeric(chr)*2.75+coord1, 
                      ymax=as.numeric(chr)*2.75+coord2),
                  color = NA, fill = colorlist[i], linewidth = 0.1, 
                  alpha = 0.5)
      next
    }
    
    subdata <- subduo(ibddata, genome, othergenome)
    subdata$coord1 <- coordsplot[i]
    subdata$coord2 <- coordsplot[i+1]
    
    gg <- gg+
      geom_rect(data=subdata[subdata$different==0,],
                aes(xmin=start, xmax=end, 
                    ymin=as.numeric(chr)*2.75+coord1, 
                    ymax=as.numeric(chr)*2.75+coord2),
                color = NA, fill = colorlist[i], linewidth = 0.1)
    
  }
  
  gg <- gg+
    scale_x_continuous(breaks = seq(0,3.25*10^6,0.25*10^6), 
                       labels = seq(0,3.25*10^6,0.25*10^6)/10^6,
                       minor_breaks = seq(0,3.25*10^6,0.125*10^6),
                       expand = c(0.01,0.01))+
    scale_y_discrete(expand = c(0,0))+
    theme(plot.background = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.grid.major.x = element_line(color = "grey60", linetype = 2, linewidth = 0.5),
          panel.grid.minor.x = element_line(color = "grey80", linetype = 2, linewidth = 0.5),
          panel.background = element_rect(fill = "grey95"),
          legend.position = "none",
          axis.title = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.ticks.x = element_line(color = "grey60", linewidth = 0.5))
  return(gg)
}

getShareCoverage <- function(chrombarcodes){
  subj <-
    with(chrombarcodes, GenomicRanges::GRanges(chr, IRanges(start, end)))
  chrombarcodes$coverage <- GenomicRanges::countOverlaps(subj, subj, type = "within")
  chrombarcodes.unique <- chrombarcodes[!duplicated(paste(chrombarcodes$chr, chrombarcodes$start, 
                                                          chrombarcodes$end, chrombarcodes$coverage)),]
  return(chrombarcodes.unique)
}

overlapIBD <- function(chrombarcodes){
  
  chrombarcodes.dt <- data.table::data.table(id=c(1:nrow(chrombarcodes)),
                                             chr = chrombarcodes$chr,
                                             start = chrombarcodes$start,
                                             end = chrombarcodes$end,
                                             coverage = chrombarcodes$coverage)
  data.table::setkey(chrombarcodes.dt, chr, start, end)
  
  #perform overlap-join
  chrom.overlaps <- data.table::foverlaps( chrombarcodes.dt, 
                                           chrombarcodes.dt,
                                           type = "any", 
                                           mult = "all", 
                                           nomatch = NULL)
  
  colnames(chrom.overlaps) <- c("chr1", "id1", "start1","end1", "coverage1",
                                "id2", "start2", "end2", "coverage2")
  
  chrom.overlaps <- chrom.overlaps[chrom.overlaps$id1 != chrom.overlaps$id2,]
  
  chrom.partlycovered.1 <- chrom.overlaps[chrom.overlaps$coverage1 <= chrom.overlaps$coverage2,]%>%
    group_by(id1)%>%
    reframe(as.matrix(interval_difference(Intervals_full(unique(matrix(c(start1, end1), ncol = 2))),
                                          Intervals_full(matrix(c(start2, end2), ncol = 2)),
                                          check_valid = TRUE)),
            chr = unique(chr1),
            coverage = unique(coverage1))
  
  chrom.partlycovered.2 <- chrom.overlaps[chrom.overlaps$coverage2 <= chrom.overlaps$coverage1,]%>%
    group_by(id2)%>%
    reframe(as.matrix(interval_difference(Intervals_full(unique(matrix(c(start2, end2), ncol = 2))),
                                          Intervals_full(matrix(c(start1, end1), ncol = 2)),
                                          check_valid = TRUE)),
            chr = unique(chr1),
            coverage = unique(coverage2))%>%
    unnest(everything())
  
  if (nrow(chrom.partlycovered.1)==0){
    chrom.partlycovered.1 <- data.frame(matrix(nrow = 0, ncol = 5))
  }
  if (nrow(chrom.partlycovered.2)==0){
    chrom.partlycovered.2 <- data.frame(matrix(nrow = 0, ncol = 5))
  }
  chrom.partlycovered.1 <- as.data.frame(as.matrix(chrom.partlycovered.1))
  chrom.partlycovered.2 <- as.data.frame(as.matrix(chrom.partlycovered.2))
  colnames(chrom.partlycovered.1) <- c("id", "start", "end", "chr", "coverage")
  colnames(chrom.partlycovered.2) <- c("id", "start", "end", "chr", "coverage")
  chrom.partlycovered <- rbind(chrom.partlycovered.1,chrom.partlycovered.2)
  
  chrom.covered.ids <- unique(c(chrom.overlaps$id2[chrom.overlaps$coverage2 < chrom.overlaps$coverage1],
                                chrom.overlaps$id1[chrom.overlaps$coverage1 < chrom.overlaps$coverage2]))
  
  chrom.overlaps.density <- chrombarcodes.dt[!chrombarcodes.dt$id %in% chrom.covered.ids,]
  chrom.overlaps.density <- rbind(chrom.overlaps.density,chrom.partlycovered)
  chrom.overlaps.density <- chrom.overlaps.density[!duplicated(chrom.overlaps.density[,c("chr", "start", "end")]),]
  
  chrom.overlaps.density <- chrom.overlaps.density[with(chrom.overlaps.density, order(chr, start)),]
  chrom.overlaps.density <- chrom.overlaps.density[,-1]
  return(chrom.overlaps.density)
  
}

#################################

