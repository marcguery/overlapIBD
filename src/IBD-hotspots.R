#############READ#############
edges <- read.table(IBD_edges_file, sep = "\t", h = T)[,c(1:3)]
colnames(edges) <- c("ID1tmp", "ID2tmp", "fract_sites_IBD")
edges$ID1 <- pmin(edges$ID1tmp, edges$ID2tmp)
edges$ID2 <- pmax(edges$ID1tmp, edges$ID2tmp)
edges <- edges[,c("ID1", "ID2", "fract_sites_IBD")]

ibdsegments <- read_delim(file = hmmIBD_file, delim = "\t")[,c(1:7)]
colnames(ibdsegments) <- c("ID1tmp", "ID2tmp", "chr", "start", "end", "different", "Nsnp")
ibdsegments$ID1 <- pmin(ibdsegments$ID1tmp, ibdsegments$ID2tmp)
ibdsegments$ID2 <- pmax(ibdsegments$ID1tmp, ibdsegments$ID2tmp)
ibdsegments <- ibdsegments[,c("ID1", "ID2", "chr", "start", "end", "different", "Nsnp")]

chromlength <- read.table(chromosomes_file)[,c(1:2)]
colnames(chromlength) <- c("chr", "length")
chromlength$chr <- as.numeric(chromlength$chr)
chromlength <- chromlength[!is.na(chromlength$chr),]

if (!is.na(highlight_regions_file)){
  highlight_regions <- read.table(highlight_regions_file, h = T)[,c(1:3)]
  colnames(highlight_regions) <- c("chr", "start", "end")
  highlight_regions$chr <- as.numeric(highlight_regions$chr)
  highlight_regions <- highlight_regions[!is.na(highlight_regions$chr),]
  highlight_regions$start <- as.numeric(highlight_regions$start)
  highlight_regions$end <- as.numeric(highlight_regions$end)
}

if (!is.na(highlight_genes_file)){
  highlight_genes <- read.table(highlight_genes_file, h = T)
  colnames(highlight_genes) <- c("chr", "start", "end", "direction", "ID", "gene", "protein")
  highlight_genes$chr <- as.numeric(highlight_genes$chr)
  highlight_genes <- highlight_genes[!is.na(highlight_genes$chr),]
  highlight_genes$start <- as.numeric(highlight_genes$start)
  highlight_genes$end <- as.numeric(highlight_genes$end)
}
##########################

#############Non-identical genomes statistics#############
edges.notsimilar <- edges[edges$fract_sites_IBD < similar,]
ibdsegments.notsimilar <- ibdsegments[paste0(ibdsegments$ID1, ibdsegments$ID2) %in% paste0(edges.notsimilar$ID1, edges.notsimilar$ID2),]
ibdsegments.notsimilar.same <- ibdsegments.notsimilar[ibdsegments.notsimilar$different == 0,]
ibdsegments.notsimilar.same.unique <- getShareCoverage(ibdsegments.notsimilar.same)
quantile.cov <- quantile(ibdsegments.notsimilar.same.unique$coverage, probs = c(0,0.01,0.05,0.1,0.25,0.5,0.75,0.9,0.95,0.99,1))

compnum <- length(which(!duplicated(ibdsegments.notsimilar[,c("ID1", "ID2")])))
print(compnum)
samenum <- length(which(!duplicated(ibdsegments.notsimilar.same[,c("ID1", "ID2")])))
##########################
#############Keeping highest percentage of genomes in IBD per interval#############
if (is.na(density_table)){
  chrom.overlaps.density <- overlapIBD(ibdsegments.notsimilar.same.unique)
  
  chrom.overlaps.connect.density <- chrom.overlaps.density %>%
    group_by(chr)%>%
    reframe(id = NA, chr = chr, end = lead(start), start = end,
              startcoverage = coverage, endcoverage = lead(coverage))
  chrom.overlaps.connect.density <- chrom.overlaps.connect.density[with(chrom.overlaps.connect.density, order(chr, start)),]
  chrom.overlaps.connect.density <- chrom.overlaps.connect.density[!is.na(chrom.overlaps.connect.density$end),]
  chrom.overlaps.connect.density <- chrom.overlaps.connect.density[,c(1,4,3,5,6)]
  write.csv(chrom.overlaps.density, paste0(outDir, "/IBD-density.csv"), quote = F, row.names = F)
}else{
  chrom.overlaps.density <- read.csv(density_table)
  
  chrom.overlaps.connect.density <- chrom.overlaps.density %>%
    group_by(chr)%>%
    reframe(id = NA, chr = chr, end = lead(start), start = end,
              startcoverage = coverage, endcoverage = lead(coverage))
  chrom.overlaps.connect.density <-chrom.overlaps.connect.density[with(chrom.overlaps.connect.density, order(chr, start)),]
  chrom.overlaps.connect.density <- chrom.overlaps.connect.density[!is.na(chrom.overlaps.connect.density$end),]
  chrom.overlaps.connect.density <- chrom.overlaps.connect.density[,c(1,4,3,5,6)]
}

##########################
#############Most covered regions#############
similar_cutoff <- quantile.cov[9]
percentile <- as.numeric(sub("%", "", names(similar_cutoff)))
100*similar_cutoff/compnum
bestcov <- chrom.overlaps.density[chrom.overlaps.density$coverage > similar_cutoff,]
curr <- bestcov[1,]
cov <- c(bestcov[1,4])
curr$maxcov <- c(bestcov[1,4])
bestcov.continuous <- curr
for (i in c(2:nrow(bestcov))){
  if(bestcov[i,1] == curr[1,1] & as.numeric(bestcov[i,2]) <= as.numeric(curr[1,3])+10000){
    curr[1,3] <- bestcov[i,3]
    cov <- c(cov, as.numeric(bestcov[i,4]))
    curr[1,4] <- round(mean(as.numeric(cov)),1)
    curr[1,5] <- max(as.numeric(cov))
  }else{
    bestcov.continuous <- rbind(bestcov.continuous, curr)
    curr <- bestcov[i,]
    curr$maxcov <- as.numeric(bestcov[i,4])
    cov <- c(as.numeric(bestcov[i,4]))}
}
bestcov.continuous <- rbind(bestcov.continuous, curr)
bestcov.continuous <- bestcov.continuous[-1,]
bestcov.continuous <- bestcov.continuous[bestcov.continuous$end - bestcov.continuous$start>=1000,]
bestcov.continuous$prctibd <- round(100*bestcov.continuous$coverage/compnum,1)
bestcov.continuous$maxprctibd <- round(100*bestcov.continuous$maxcov/compnum,1)
bestcov.continuous
write.csv(bestcov.continuous, paste0(outDir, "/", 100-percentile, "prct-most-covered.csv"), quote = F, row.names = F)

chrom.overlaps.density.mostcovered <- chrom.overlaps.density[1,]
chrom.overlaps.density.mostcovered$id <- 0
chrom.overlaps.connect.density.mostcovered <- chrom.overlaps.connect.density[1,]
chrom.overlaps.connect.density.mostcovered$id <- 0
for (i in c(1:nrow(bestcov.continuous))){
  curr <- bestcov.continuous[i,]
  mostcovered <- chrom.overlaps.density[chrom.overlaps.density$chr == curr$chr & 
                                          chrom.overlaps.density$end >= as.numeric(curr$start)-25000 &
                                          chrom.overlaps.density$start <= as.numeric(curr$end)+25000,]
  mostcovered$id <- i
  chrom.overlaps.density.mostcovered <- rbind(chrom.overlaps.density.mostcovered, mostcovered)
  mostcovered <- chrom.overlaps.connect.density[chrom.overlaps.connect.density$chr == curr$chr & 
                                                  chrom.overlaps.connect.density$end >= as.numeric(curr$start)-25000 &
                                                  chrom.overlaps.connect.density$start <= as.numeric(curr$end)+25000,]
  mostcovered$id <- i
  chrom.overlaps.connect.density.mostcovered <- rbind(chrom.overlaps.connect.density.mostcovered, mostcovered)
}
chrom.overlaps.density.mostcovered <- chrom.overlaps.density.mostcovered[-1,]
chrom.overlaps.connect.density.mostcovered <- chrom.overlaps.connect.density.mostcovered[-1,]

##########################
#############Plots#############
upper_lim <- 1.12*max(bestcov.continuous$maxprctibd)/100
gg <- ggplot(chrom.overlaps.density.mostcovered)+
  geom_rect(data = bestcov.continuous,
            aes(xmin = start, xmax = end,
                ymin = 0, ymax = 1.01*maxprctibd/100), 
            fill = "grey30", color = NA,
            alpha = 0.2)+
  geom_hline(data = data.frame(1),
             aes(color = "green", yintercept = similar_cutoff/compnum),
             linetype = 2, alpha = 0.6, show.legend = T)+
  geom_hline(data = data.frame(1),
             aes(color = "red", yintercept = quantile.cov[8]/compnum),
             linetype = 2, alpha = 0.6, show.legend = T)+
  geom_segment(data = highlight_genes[highlight_genes$gene%in%c("aat1", "crt"),],
               aes(x = (as.numeric(start)+as.numeric(end))/2, xend = (as.numeric(start)+as.numeric(end))/2,
                   y = 0, yend = 1.025*max(chrom.overlaps.density$coverage)/compnum),
               linetype = 2 , color = "blue", size = 0.5)+
  geom_text(data=highlight_genes[highlight_genes$gene%in%c("aat1", "crt"),],
            aes(x = as.numeric(start)-2000,
                y=1.03*max(chrom.overlaps.density$coverage)/compnum, label = gene),
            color = "purple", size = 6, fontface = "bold.italic",hjust = 1)+
  geom_point(data = highlight_genes[highlight_genes$gene%in%c("aat1", "crt"),],
             aes(x = (as.numeric(start)+as.numeric(end))/2,
                 y = 1.025*max(chrom.overlaps.density$coverage)/compnum),
             color = "blue")+
  geom_segment(aes(x = start,
                   xend = end,
                   y = coverage/compnum,
                   yend = coverage/compnum))+
  geom_segment(data = chrom.overlaps.connect.density.mostcovered,
               aes(x = start,
                   xend = end,
                   y = startcoverage/compnum,
                   yend = endcoverage/compnum))+
  scale_y_continuous(labels = function(x){x*100},
                     n.breaks = 5, expand = c(0,0),
                     limits = c(0,upper_lim))+
  scale_x_continuous(labels = function(x){x/1000},
                     n.breaks = 20, expand = c(0.01,0.01))+
  scale_color_manual("Most covered regions",
                     breaks = c("green", "red"),
                     labels = c("Top 5%", "Top 10%"),
                     values = c("green", "red"))+
  xlab("Chromosomal position (kb)")+
  ylab("Percentage of pairs of genomes in IBD")+
  theme(plot.background = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_line(color = "grey60", linetype = 3, size = 0.5),
        panel.grid.minor.x = element_line(color = "grey80", linetype = 3, size = 0.5),
        panel.background = element_rect(fill = "white"),
        legend.position = "top",
        legend.background = element_rect(color = "black"), 
        legend.direction = "horizontal",
        text = element_text(size = 17),
        axis.ticks.x = element_line(color = "grey60", size = 0.5),
        axis.ticks.y = element_line(color = "grey60", size = 0.5),
        strip.background=element_rect(fill="grey90"))+
  facet_wrap(vars(chr), scales = "free_x", ncol = 1)
gg

ggsave(paste0(outDir, "/IBD-most-covered.png"), 
       plot = gg, width = 12, height = 10, dpi = 400)


chr.maxend <- chrom.overlaps.density %>%
  group_by(chr)%>%
  summarise(maxvalue = max(end))
for (i in c(1:14)){
  if(!i%in%chr.maxend$chr){
    chr.maxend <- rbind(chr.maxend, c(i, 0))
  }
}

chr.maxend$shiftvalue <- chr.maxend$maxvalue+25000
chr.maxend <- rbind(chr.maxend, 0)
chr.maxend <- chr.maxend[with(chr.maxend, order(chr)),]
chr.maxend$shiftvalue <- cumsum(chr.maxend$shiftvalue )

chrom.overlaps.density.linedup <- data.frame(t(apply(chrom.overlaps.density,1,
                                                     FUN = function(x) {
                                                       x[2] <- x[2] + chr.maxend$shiftvalue[chr.maxend$chr == (x[1]-1)]
                                                       x[3] <- x[3] + chr.maxend$shiftvalue[chr.maxend$chr == (x[1]-1)]
                                                       return(x)})))

chrom.overlaps.connect.density.linedup <- data.frame(t(apply(chrom.overlaps.connect.density,1,
                                                             FUN = function(x) {
                                                               x[2] <- x[2] + chr.maxend$shiftvalue[chr.maxend$chr == (x[1]-1)]
                                                               x[3] <- x[3] + chr.maxend$shiftvalue[chr.maxend$chr == (x[1]-1)]
                                                               return(x)})))
highlight_genes.linedup <- data.frame(t(apply(highlight_genes,1,
                                        FUN = function(x) {
                                          x[2] <- as.numeric(x[2]) + chr.maxend$shiftvalue[chr.maxend$chr == (as.numeric(x[1])-1)]
                                          x[3] <- as.numeric(x[3]) + chr.maxend$shiftvalue[chr.maxend$chr == (as.numeric(x[1])-1)]
                                          return(x)})))
highlight_genes.linedup.mostcovered <- highlight_genes.linedup[highlight_genes.linedup$gene%in%c("aat1", "crt"),]
highlight_genes.linedup.mostcovered$chr <- as.numeric(highlight_genes.linedup.mostcovered$chr)
splitchr <- list("firsthaflfchr" = c(0:9),
                 "secondhalfchr" = c(9:14))
for (split in splitchr){
  gg <- ggplot(chrom.overlaps.density.linedup[chrom.overlaps.density.linedup$chr%in%split[-1],])+
    geom_vline(data = chr.maxend[chr.maxend$chr%in%split,][-1,],
               aes(xintercept = shiftvalue+25000/2),
               linetype = 2, color = "grey70", size = 0.5)+
    geom_segment(data = highlight_genes.linedup[as.numeric(highlight_genes.linedup$chr)%in%split[-1],],
                 aes(x = (as.numeric(start)+as.numeric(end))/2, xend = (as.numeric(start)+as.numeric(end))/2,
                     y = 0, yend = 1.05*max(chrom.overlaps.density$coverage)/compnum),
                 linetype = 3 , color = "blue", size = 0.3)+
    geom_text(data=highlight_genes.linedup[as.numeric(highlight_genes.linedup$chr)%in%split[-1],],
              aes(x = as.numeric(start)-7500,
                  y=1.065*max(chrom.overlaps.density$coverage)/compnum, label = gene),
              color = "purple", size = 3.5, fontface = "bold.italic",hjust = 1)+
    geom_point(data = highlight_genes.linedup[as.numeric(highlight_genes.linedup$chr)%in%split[-1],],
               aes(x = (as.numeric(start)+as.numeric(end))/2,
                   y = 1.05*max(chrom.overlaps.density$coverage)/compnum),
               color = "blue")+
    geom_segment(aes(x = start, xend = end,
                     y = coverage/compnum, yend = coverage/compnum,
                     color = coverage > quantile.cov[9]),
                 size = 0.25
    )+
    geom_segment(data = chrom.overlaps.connect.density.linedup[chrom.overlaps.connect.density.linedup$chr%in%split[-1],],
                 aes(x = start, xend = end,
                     y = startcoverage/compnum, yend = endcoverage/compnum,
                     color = startcoverage > quantile.cov[9] & endcoverage > quantile.cov[9]),
                 size = 0.25)+
    scale_color_manual("Coverage",
                       breaks = c("TRUE", "FALSE"),
                       values = c("black", "black"))+
    scale_y_continuous(labels = function(x){x*100},
                       n.breaks = 10, expand = c(0,0),
                       limits = c(0,upper_lim))+
    scale_x_continuous(breaks = (chr.maxend$shiftvalue[chr.maxend$chr%in%split]+
                                   lag(chr.maxend$shiftvalue[chr.maxend$chr%in%split])-125000)/2,
                       labels = chr.maxend$chr[chr.maxend$chr%in%split],
                       expand = c(0,0))+
    xlab("Chromosome")+
    ylab("Percentage of pairs of genome in IBD")+
    theme(panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          legend.position = "none")
  gg
  ggsave(paste0(outDir, "/IBD-coverage-", min(split)+1, "_", max(split), ".png"),
         plot = gg, width = 16, height = 8, dpi = 250)
}

gg <- ggplot(chrom.overlaps.density)+
  geom_segment(aes(x = start,
                   xend = end,
                   y = coverage/compnum,
                   yend = coverage/compnum))+
  geom_segment(data = chrom.overlaps.connect.density,
               aes(x = start,
                   xend = end,
                   y = startcoverage/compnum,
                   yend = endcoverage/compnum))+
  scale_y_continuous(labels = function(x){x*100})+
  scale_x_continuous(labels = function(x){x/1000},
                     n.breaks = 10)+
  xlab("Chromosomal position (kb)")+
  ylab("Percentage of pairs of genome in IBD")+
  facet_wrap(vars(chr), scales = "free", ncol = 3)
gg

ggsave(paste0(outDir, "/IBD-coverage-split.png"), 
       plot = gg, width = 16, height = 12, dpi = 250)

chrom.overlaps.density.capped <- chrom.overlaps.density
chrom.overlaps.density.capped$coverage[chrom.overlaps.density.capped$coverage/compnum > uppercov/100] <- (uppercov/100)*compnum
chrom.overlaps.density.capped$coverage[chrom.overlaps.density.capped$coverage/compnum < lowercov/100] <- 0

gg <- ggplot()+
  geom_blank(aes(fill = factor(c(lowercov, c(ceiling(lowercov):uppercov)))))+
  geom_rect(data=chromlength,
            aes(xmin = 0, xmax=length, 
                ymin=as.numeric(chr)*2.75-0.5, ymax=as.numeric(chr)*2.75+0.5), 
            color=NA, fill="grey80")+
  geom_text(data=chromlength,
            aes(x = 0-100000, y=as.numeric(chr)*2.75, 
                label=chr), 
            color="grey50", size=5)+
  geom_rect(data=chrom.overlaps.density.capped[chrom.overlaps.density.capped$coverage > 0,],
            aes(xmin=start, xmax=end-1, 
                ymin=as.numeric(chr)*2.75-0.5, ymax=as.numeric(chr)*2.75+0.5,
                fill = factor(pmax(lowercov,floor(coverage*100/compnum)))),
            linewidth = 0, color = NA)
if (!is.na(highlight_genes_file)){
  gg <- gg+
    geom_text(data=highlight_genes,
            aes(x = start-10000,
                y=chr*2.75+1.3, label = gene),
            color = "purple", size = 4, fontface = "bold.italic",hjust = 1)+
    geom_point(data=highlight_genes,
             aes(x = (end+start)/2,
                 y=chr*2.75+0.75),
             color = "purple", size = 1.5)
  }

if (!is.na(highlight_regions_file)){
  gg <- gg+
    geom_rect(data=highlight_regions,
              aes(xmin = start, xmax = end,
                  ymin=chr*2.75-0.57, ymax=chr*2.75+0.57),
              color = "blue", fill = NA, linetype = 5, size = 0.4)
}

gg <- gg+
  scale_fill_brewer("Percentage of pair of genomes in IBD",
                    type = "div", palette = "RdYlGn",
                    labels = function(x){
                      numbers <- as.numeric(x)
                      paste(numbers,
                            ifelse(is.na(lead(numbers)),"+",
                                   paste("-",
                                         lead(numbers))))})+
  scale_x_continuous(breaks = seq(0,3.25*10^6,0.1*10^6), 
                     labels = seq(0,3.25*10^6,0.1*10^6)/10^6,
                     minor_breaks = seq(0,3.25*10^6,0.025*10^6),
                     expand = c(0.01,0.05))+
  scale_y_discrete(expand = c(0.01,0.01))+
  theme(plot.background = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_line(color = "grey60", linetype = 3, size = 0.5),
        panel.grid.minor.x = element_line(color = "grey80", linetype = 3, size = 0.5),
        panel.background = element_rect(fill = "white"),
        legend.position = "inside",
        legend.position.inside = c(0.7,0.15),
        legend.background = element_rect(color = "black"), 
        legend.direction = "horizontal",
        axis.title = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_line(color = "grey60", size = 0.5))
gg
ggsave(paste0(outDir, "/conserved-regions.png"),
         plot = gg, width = 11, height = 6.5, dpi = 800)

##########################

