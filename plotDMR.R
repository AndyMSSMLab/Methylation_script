library(tidyverse)
library(data.table)
library(reshape2)

plotDMR <- function(data, samples, startPos, endPos, chrCol, startCol,title = ""){
  data = data %>% slice(startPos:endPos)
  epivarSamples = data %>% select(Sign_individuals_t0.1_n3_w1k_BOTH) %>%
    filter(!is.na(Sign_individuals_t0.1_n3_w1k_BOTH)) %>%
    separate_rows(Sign_individuals_t0.1_n3_w1k_BOTH, sep = ",") %>%
    distinct %>% pull
  
  data = data %>% 
    melt(id.vars = setdiff(colnames(data), samples), variable.name = "sampleID", value.name = "Beta") %>%
    mutate(status = ifelse(sampleID %in% epivarSamples, "DMR-carrier", "normal"))
  
  g = ggplot(data =  data, aes_string(x=startCol, y="Beta")) +
    geom_line(aes(group=sampleID,size=status, linetype = status, colour = status)) +
    scale_colour_manual(values = c("DMR-carrier" = "red","normal" = "black")) +
    scale_size_manual(values = c("DMR-carrier"=0.9, "normal" = 0.2 )) +
    scale_linetype_manual(values = c("DMR-carrier"="solid", "normal" = "dashed" )) +
    coord_cartesian(ylim = c(0,1)) +
    xlab(paste0("Genomic Position at ",data[1,chrCol])) +
    ylab("Methylation Value") +
    labs(title = title) +
    theme_bw()
  g
}

data = fread("example.epi.out.txt", sep = "\t", header = T) %>% as.data.frame %>%
  filter(Sign_window_t0.1_n3_w1k_BOTH ==1)

dmr = findDMR(data)

##List of all samples
samples = c("Sample1","Sample2","Sample3","Sample4","Sample5")

###Total number of DMRs can be found with length(dmr$dmrStart)
i = 1   ## First DMR
plotDMR(data, samples, dmr$dmrStart[i], dmr$dmrEnd[i], "Chr_CpG", "Start_CpG")

