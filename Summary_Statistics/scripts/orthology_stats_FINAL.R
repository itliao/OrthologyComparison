work_dir <-"/path_to_directory/"
setwd (work_dir)

library(ggplot2)
library(lattice)
library(car)

library(dplyr)
library(tidyr)
library(tibble)
library(tidyverse)

library(openxlsx)
library(RColorBrewer)

########################
# directories and path #
########################
fullDir = paste0(work_dir, "counts/fullBrass/")
diploidDir = paste0(work_dir, "counts/diploidBrass/")
listDirectory = list(c("BR","broccoli"),
                     c("OF_blast","orthofinder"),
                     c("OF_diamond","orthofinder"),
                     c("OF_mmseqs","orthofinder"),
                     c("ON", "orthnet"),
                     c("SP_diamond","sonicparanoid"),
                     c("SP_mmseqs","sonicparanoid"))

########################################################
# diploid set - reformat files for downstream analyses #
########################################################

numSpOG1 <- list()
numGenesSpOG1 <- list()
numGenesSpOG2 <- list()

for (i in listDirectory){
  category <- i[1]
  program <- i[2]
  
  allCountsPath <- paste0(diploidDir, category, "/", program, "_allCounts.tsv")
  binaryPath <- paste0(diploidDir, category, "/", program, "_binary.tsv")
  tenMaxPath <- paste0(diploidDir, category, "/", program, "_tenMax.tsv")
  
  #read in files
  allCounts <- read.table(allCountsPath, header=T)
  binary <- read.table(binaryPath, header=T)
  tenMax <- read.table(tenMaxPath, header=T)
  
  #create new column - binary - for number of species/OG
  binary$numPerOG = apply(binary[,c(2:6)], 1, sum)
  
  if (program == "orthofinder" || program == "orthnet"){
    onePlus = binary[binary$numPerOG != 1, ]
    onePlus$program <- category
    
    df_2column <- onePlus[,c('program','numPerOG')]
    numSpOG1 <- append(numSpOG1, list(df_2column))
  } else {
    binary$program <- category
    
    df_2column <- binary[,c('program','numPerOG')]
    numSpOG1 <- append(numSpOG1, list(df_2column))
  }
    
  # for tenMax - for number of genes/species per OG
  if (program == "orthofinder" || program == "orthnet"){
    onePlus = binary[binary$numPerOG != 1, ]
    
    if(program == "orthofinder"){
      mergedDF <- merge(onePlus,tenMax,by="HOG")
    }else{
      mergedDF <- merge(onePlus,tenMax,by="groupNumber")
    }
    
    # remove columns from onePlus
    OGtenMaxOnePlus <- mergedDF[-c(2:6)]
    
    # rename columns
    names(OGtenMaxOnePlus)[names(OGtenMaxOnePlus) == "Ath.y"] <- "Ath"
    names(OGtenMaxOnePlus)[names(OGtenMaxOnePlus) == "Chi.y"] <- "Chi"
    names(OGtenMaxOnePlus)[names(OGtenMaxOnePlus) == "Cru.y"] <- "Cru"
    names(OGtenMaxOnePlus)[names(OGtenMaxOnePlus) == "Tar.y"] <- "Tar"
    names(OGtenMaxOnePlus)[names(OGtenMaxOnePlus) == "Aar.y"] <- "Aar"
    
    #table format
    OGtenMaxOnePlus$program <- category
    numGenesSpOG2 <- append(numGenesSpOG2, list(OGtenMaxOnePlus))
    
    OGtenMaxOnePlus %>% 
      select("Ath", "Chi", "Cru", "Tar", "Aar") %>%
      gather(key = species , value = numGenes) %>%
      group_by(species) %>%
      count(numGenes) %>%
      ungroup %>%
      spread(species, n) -> geneSum1
    
    # graphical format
    OGtenMaxOnePlus %>% 
      select("Ath", "Chi", "Cru", "Tar", "Aar") %>%
      gather(key = species , value = numGenes) %>%
      group_by(species) %>%
      count(numGenes) -> geneSum2
    
    names(geneSum2)[3] = category
    numGenesSpOG1 <- append(numGenesSpOG1, list(geneSum2))
    
  } else {
    # table format
    tenMax$program <- category
    numGenesSpOG2 <- append(numGenesSpOG2, list(tenMax))
    
    tenMax %>% 
      select("Ath", "Chi", "Cru", "Tar", "Aar") %>%
      gather(key = species , value = numGenes) %>%
      group_by(species) %>%
      count(numGenes) %>%
      ungroup %>%
      spread(species, n) -> geneSum1
    
    # graphical format
    tenMax %>% 
      select("Ath", "Chi", "Cru", "Tar", "Aar") %>%
      gather(key = species , value = numGenes) %>%
      group_by(species) %>%
      count(numGenes) -> geneSum2
    
    names(geneSum2)[3] = category
    numGenesSpOG1 <- append(numGenesSpOG1, list(geneSum2))
  }
}

merged1 <- bind_rows(numSpOG1) 
merged2 <- bind_cols(numGenesSpOG1) 
merged3 <- bind_rows(numGenesSpOG2) 

########################################################################
# diploid - stacked bar plot for the distributions of all the algorithms #
########################################################################

merged1 %>%
  gather(key = program, value = numPerOG) %>% 
  group_by(program) %>%
  count(numPerOG) -> sumSpOG

sumSpOG$program <- factor(sumSpOG$program, 
                          levels = c("BR", "OF_blast", "OF_diamond", "OF_mmseqs",
                                     "SP_diamond", "SP_mmseqs", "ON"))


numSpOG_pdf <- paste0(work_dir, "comparison/230928_combine/diploid_numSpPerOG_230928.pdf") 
pdf(numSpOG_pdf)

p1 <- ggplot(sumSpOG, aes(fill=factor(numPerOG), y=n, x=program)) + 
  geom_bar(position="stack", stat="identity") +
  scale_fill_brewer(palette="BrBG")+
  ggtitle("number of species in OG - stacked bars") +
  xlab("Program") +
  ylab("number of species per OG") +
  theme_bw()

p2 <- ggplot(sumSpOG, aes(fill=factor(numPerOG), y=n, x=program)) + 
  geom_bar(position="fill", stat="identity") +
  scale_fill_brewer(palette="BrBG")+
  ggtitle("number of species in OG - stacked bars") +
  xlab("Program") +
  ylab("porportion of the number of species per OG") +
  theme_bw()

print(p1)
print(p2)

dev.off()

# summary statistic
merged1 %>%
  group_by(program) %>%
  summarise(
    n = n(),
    mean = mean(numPerOG),
    median = median(numPerOG),
    sd   = sd(numPerOG),
    se = sd/sqrt(n),
    mean_se_high = mean + se,
    mean_se_low = mean - se) -> numSpPerOG_sum

numSpPerOG_sum
numSpPerOG_df <- as.data.frame(numSpPerOG_sum)
numSpOG_df_out <- paste0(work_dir, "comparison/230928_combine/diploid_numSpPerOG_sumstats_230928.csv") 
write.csv(numSpPerOG_df,numSpOG_df_out)

######################################################################
# diploid - violin plot with average number of species in orthogroup #
######################################################################
brewer.pal(9, "Purples")
brewer.pal(9, "Blues")

#"#FCFBFD" "#EFEDF5" "#DADAEB" "#BCBDDC" "#9E9AC8" "#807DBA" "#6A51A3" "#54278F" "#3F007D"
#"#F7FBFF" "#DEEBF7" "#C6DBEF" "#9ECAE1" "#6BAED6" "#4292C6" "#2171B5" "#08519C" "#08306B"

violin_plot <- paste0(work_dir,"comparison/230928_combine/diploid_numSpPerOG_violin_230928.pdf")
pdf(violin_plot)
merged1$program <- factor(merged1$program, 
                          levels = c("BR", "OF_blast", "OF_diamond", "OF_mmseqs",
                                     "SP_diamond", "SP_mmseqs", "ON"))
p2 <- ggplot(merged1, aes(x=program, y=numPerOG, color = program)) + ylim(2,5) +
  geom_violin(width=0.5) +
  geom_pointrange(aes(x = program, y = mean, ymin = mean_se_low, ymax = mean_se_high , color=program),
                  data = numSpPerOG_sum,
                  position = position_nudge(x=0.25)) +
  scale_color_manual(values = c("#08519C","#2171B5","#4292C6","#807DBA","#6A51A3","#54278F","#3F007D")) +
  labs(x = "program", y = "number of species per orthogroup") +
  theme_bw()

print(p2)
dev.off()

##################################################################################################
# diploid -  statistics (Kruskal-Wallis, non-parametric) for numbers of species in an orthogroup #
##################################################################################################
kwtest_out <- paste0(work_dir,"comparison/230928_combine/diploid_numSpPerOG_kwTest_230928.txt")
sink(kwtest_out)
kruskal.test(numPerOG ~ program, data = merged1)
pairwise.wilcox.test(merged1$numPerOG, merged1$program,
                     p.adjust.method = "fdr")

merged1_noON <- merged1[!(merged1$program=="ON"),]
kruskal.test(numPerOG ~ program, data = merged1_noON)
pairwise.wilcox.test(merged1_noON$numPerOG, merged1_noON$program,
                     p.adjust.method = "fdr")
sink()

######################################################################################
# diploid - Rearrange bar plots for number of genes/species in OG - group by species #
######################################################################################

geneSpOG_pdf <- paste0(work_dir,"comparison/230928_combine/diploid_numGenesPerSpOG_230928.pdf")
pdf(geneSpOG_pdf)
listOfSpecies <- c("Ath", "Chi", "Cru", "Tar", "Aar")
for (species2 in listOfSpecies){
  targetSp <- subset(merged2, species...1==species2)
  
  numGenesPerSp <- c(0:11)
  targetSp %>%
    select("BR", "OF_blast", "OF_diamond", "OF_mmseqs", "SP_diamond", "SP_mmseqs", "ON") %>%
    gather(key = program, value = numGenes) %>% 
    group_by(program) %>%
    mutate(numGenesSp = numGenesPerSp[numGenesPerSp])-> geneMatrix1
  
  geneMatrix1$program <- factor(geneMatrix1$program,
                                levels = c("BR", "OF_blast", "OF_diamond", "OF_mmseqs", 
                                           "SP_diamond", "SP_mmseqs", "ON"))
  print(
    ggplot(geneMatrix1, aes(fill=factor(numGenesSp), y=numGenes, x=program)) + 
      geom_bar(position="stack", stat="identity") +
      scale_fill_brewer(palette="PRGn") +
      ggtitle(paste0("number of genes per species in OG:", species2)) +
      xlab("Program") +
      ylab("number of genes per species per OG") +
      theme_bw()
  )
  
  print(
    ggplot(geneMatrix1, aes(fill=factor(numGenesSp), y=numGenes, x=program)) + 
      geom_bar(position="fill", stat="identity") +
      scale_fill_brewer(palette="PRGn") +
      ggtitle(paste0("number of genes per species in OG:", species2)) +
      xlab("Program") +
      ylab("number of genes per species per OG") +
      theme_bw()
  )
}

dev.off()

##############################################################################################
# diploid - statistics (Kruskal-Wallis) for number of genes/species in OG - group by species #
##############################################################################################

#stats - Ath
merged3 %>%
  select("program", "Ath") %>%
  gather(key = program, value = Ath) %>% 
  group_by(program) %>%
  count(Ath) -> sumMatrix_ath

ktest_ath <- kruskal.test(Ath ~ program, data = merged3)
wilcox_ath <- pairwise.wilcox.test(merged3$Ath, merged3$program, p.adjust.method = "fdr")
wilcox_ath_df <- as.data.frame(wilcox_ath$p.value)

#stats - Chi
merged3 %>%
  select("program", "Chi") %>%
  gather(key = program, value = Chi) %>% 
  group_by(program) %>%
  count(Chi) -> sumMatrix_chi

ktest_chi <- kruskal.test(Chi ~ program, data = merged3)
wilcox_chi <- pairwise.wilcox.test(merged3$Chi, merged3$program, p.adjust.method = "fdr")
wilcox_chi_df <- as.data.frame(wilcox_chi$p.value)

#stats - Cru
merged3 %>%
  select("program", "Cru") %>%
  gather(key = program, value = Cru) %>% 
  group_by(program) %>%
  count(Cru) -> sumMatrix_cru

ktest_cru <- kruskal.test(Cru ~ program, data = merged3)
wilcox_cru <- pairwise.wilcox.test(merged3$Cru, merged3$program, p.adjust.method = "fdr")
wilcox_cru_df <- as.data.frame(wilcox_cru$p.value)

#stats - Tar
merged3 %>%
  select("program", "Tar") %>%
  gather(key = program, value = Tar) %>% 
  group_by(program) %>%
  count(Tar) -> sumMatrix_tar

ktest_tar <- kruskal.test(Tar ~ program, data = merged3)
wilcox_tar <- pairwise.wilcox.test(merged3$Tar, merged3$program, p.adjust.method = "fdr")
wilcox_tar_df <- as.data.frame(wilcox_tar$p.value)

#stats - Aar
merged3 %>%
  select("program", "Aar") %>%
  gather(key = program, value = Aar) %>% 
  group_by(program) %>%
  count(Aar) -> sumMatrix_aar

ktest_aar <- kruskal.test(Aar ~ program, data = merged3)
wilcox_aar <- pairwise.wilcox.test(merged3$Aar, merged3$program, p.adjust.method = "fdr")
wilcox_aar_df <- as.data.frame(wilcox_aar$p.value)

geneSpOG <- c(sumMatrix_ath, sumMatrix_chi, sumMatrix_cru, sumMatrix_tar, sumMatrix_aar)
geneSpOG_merge <- bind_cols(geneSpOG)
wilcox_list <- c(wilcox_ath_df, wilcox_chi_df,wilcox_cru_df, wilcox_tar_df, wilcox_aar_df)
wilcox_all <- bind_cols(wilcox_list)

geneSpOG_out <- paste0(work_dir,"comparison/230928_combine/diploid_numGenesPerSpOG_sumStats_230928.csv")
write.csv(geneSpOG_merge, geneSpOG_out)
wilcox_out <- paste0(work_dir,"comparison/230928_combine/diploid_numGenesPerSpOG_wilcox_230928.csv")
write.csv(wilcox_all, wilcox_out)

kwtest_out <- paste0(work_dir,"comparison/230928_combine/diploid_numGenesPerSpOG_kwTest_230928.txt")
sink(kwtest_out)
print(ktest_ath)
print(ktest_chi)
print(ktest_cru)
print(ktest_tar)
print(ktest_aar)
sink()

#####################################################################
# higher ploidy "full" set - reformat files for downstream analyses #
#####################################################################
numSpOG2 <- list()
numGenesSpOG3 <- list()
numGenesSpOG4 <- list()

for (i in listDirectory){
  category <- i[1]
  program <- i[2]
    
  allCountsPath <- paste0(fullDir, category, "/", program, "_allCounts.tsv")
  binaryPath <- paste0(fullDir, category, "/", program, "_binary.tsv")
  tenMaxPath <- paste0(fullDir, category, "/", program, "_tenMax.tsv")
    
  #read in files
  allCounts <- read.table(allCountsPath, header=T)
  binary <- read.table(binaryPath, header=T)
  tenMax <- read.table(tenMaxPath, header=T)
    
  #create new column - binary - for number of species/OG
  binary$numPerOG = apply(binary[,c(2:9)], 1, sum)
    
  if (program == "orthofinder" || program == "orthnet"){
    onePlus = binary[binary$numPerOG != 1, ]
    onePlus$program <- category
    
    df_2column <- onePlus[,c('program','numPerOG')]
    numSpOG2 <- append(numSpOG2, list(df_2column))
    } else {
      binary$program <- category
      
      df_2column <- binary[,c('program','numPerOG')]
      numSpOG2 <- append(numSpOG2, list(df_2column))
    }
  
  # for tenMax - for number of genes/species per OG
  if (program == "orthofinder" || program == "orthnet"){
    onePlus = binary[binary$numPerOG != 1, ]
    
    if(program == "orthofinder"){
      mergedDF <- merge(onePlus,tenMax,by="HOG")
    }else{
      mergedDF <- merge(onePlus,tenMax,by="groupNumber")
    }
    
    # remove columns from onePlus
    OGtenMaxOnePlus <- mergedDF[-c(2:9)]
    
    # rename columns
    names(OGtenMaxOnePlus)[names(OGtenMaxOnePlus) == "Ath.y"] <- "Ath"
    names(OGtenMaxOnePlus)[names(OGtenMaxOnePlus) == "Bra.y"] <- "Bra"
    names(OGtenMaxOnePlus)[names(OGtenMaxOnePlus) == "Chi.y"] <- "Chi"
    names(OGtenMaxOnePlus)[names(OGtenMaxOnePlus) == "Cru.y"] <- "Cru"
    names(OGtenMaxOnePlus)[names(OGtenMaxOnePlus) == "Csa.y"] <- "Csa"
    names(OGtenMaxOnePlus)[names(OGtenMaxOnePlus) == "Sal.y"] <- "Sal"
    names(OGtenMaxOnePlus)[names(OGtenMaxOnePlus) == "Tar.y"] <- "Tar"
    names(OGtenMaxOnePlus)[names(OGtenMaxOnePlus) == "Aar.y"] <- "Aar"
    
    #table format
    OGtenMaxOnePlus$program <- category
    numGenesSpOG4 <- append(numGenesSpOG4, list(OGtenMaxOnePlus))

    OGtenMaxOnePlus %>% 
      select("Ath", "Chi", "Cru", "Tar", "Bra", "Csa", "Sal", "Aar") %>%
      gather(key = species , value = numGenes) %>%
      group_by(species) %>%
      count(numGenes) %>%
      ungroup %>%
      spread(species, n) -> geneSum1
    
    geneSum1$program <- category
    numGenesSpOG1 <- append(numGenesSpOG1, list(geneSum1))
    
    # graphical format
    OGtenMaxOnePlus %>% 
      select("Ath", "Chi", "Cru", "Tar", "Bra", "Csa", "Sal", "Aar") %>%
      gather(key = species , value = numGenes) %>%
      group_by(species) %>%
      count(numGenes) -> geneSum2
    
    names(geneSum2)[3] = category
    numGenesSpOG3 <- append(numGenesSpOG3, list(geneSum2))
    
  } else {
    # table format
    tenMax$program <- category
    numGenesSpOG4 <- append(numGenesSpOG4, list(tenMax))
    
    tenMax %>% 
      select("Ath", "Chi", "Cru", "Tar", "Bra", "Csa", "Sal", "Aar") %>%
      gather(key = species , value = numGenes) %>%
      group_by(species) %>%
      count(numGenes) %>%
      ungroup %>%
      spread(species, n) -> geneSum1
    
    geneSum1$program <- category
    numGenesSpOG1 <- append(numGenesSpOG1, list(geneSum1))
    
    # graphical format
    tenMax %>% 
      select("Ath", "Chi", "Cru", "Tar", "Bra", "Csa", "Sal", "Aar") %>%
      gather(key = species , value = numGenes) %>%
      group_by(species) %>%
      count(numGenes) -> geneSum2
    
    names(geneSum2)[3] = category
    numGenesSpOG3 <- append(numGenesSpOG3, list(geneSum2))
  }
}

merged4 <- bind_cols(numGenesSpOG3) 
merged5 <- bind_rows(numSpOG2) 
merged6 <- bind_rows(numGenesSpOG4)

#######################################################################
# full - stacked bar plot for the distributions of all the algorithms #
#######################################################################

merged5 %>%
  gather(key = program, value = numPerOG) %>% 
  group_by(program) %>%
  count(numPerOG) -> sumSpOG

sumSpOG$program <- factor(sumSpOG$program, 
                          levels = c("BR", "OF_blast", "OF_diamond", "OF_mmseqs",
                                     "SP_diamond", "SP_mmseqs", "ON"))

numSpOG_pdf <- paste0(work_dir, "comparison/230928_combine/full_numSpPerOG_230928.pdf") 
pdf(numSpOG_pdf)

p1 <- ggplot(sumSpOG, aes(fill=factor(numPerOG), y=n, x=program)) + 
  geom_bar(position="stack", stat="identity") +
  scale_fill_brewer(palette="BrBG")+
  ggtitle("number of species in OG - stacked bars") +
  xlab("Program") +
  ylab("number of species per OG") +
  theme_bw()

p2 <- ggplot(sumSpOG, aes(fill=factor(numPerOG), y=n, x=program)) + 
  geom_bar(position="fill", stat="identity") +
  scale_fill_brewer(palette="BrBG")+
  ggtitle("number of species in OG - stacked bars") +
  xlab("Program") +
  ylab("proportion of the number of species per OG") +
  theme_bw()

print(p1)
print(p2)
dev.off()

# summary statistic
merged5 %>%
  group_by(program) %>%
  summarise(
    n = n(),
    mean = mean(numPerOG),
    median = median(numPerOG),
    sd   = sd(numPerOG),
    se = sd/sqrt(n),
    mean_se_high = mean + se,
    mean_se_low = mean - se) -> numSpPerOG_sum

numSpPerOG_sum
numSpPerOG_df <- as.data.frame(numSpPerOG_sum)
numSpOG_df_out <- paste0(work_dir, "comparison/230928_combine/full_numSpPerOG_sumstats_230928.csv")
write.csv(numSpPerOG_df,numSpOG_df_out)

###################################################################
# full - violin plot with average number of species in orthogroup #
###################################################################
brewer.pal(9, "Purples")
brewer.pal(9, "Blues")

#"#FCFBFD" "#EFEDF5" "#DADAEB" "#BCBDDC" "#9E9AC8" "#807DBA" "#6A51A3" "#54278F" "#3F007D"
#"#F7FBFF" "#DEEBF7" "#C6DBEF" "#9ECAE1" "#6BAED6" "#4292C6" "#2171B5" "#08519C" "#08306B"

violin_plot <- paste0(work_dir,"comparison/230928_combine/full_numSpPerOG_violin_230928.pdf")
pdf(violin_plot)
merged5$program <- factor(merged5$program, 
                          levels = c("BR", "OF_blast", "OF_diamond", "OF_mmseqs",
                                     "SP_diamond", "SP_mmseqs", "ON"))
p2 <- ggplot(merged5, aes(x=program, y=numPerOG, color = program)) + ylim(2,8) +
  geom_violin(width=0.5) +
  geom_pointrange(aes(x = program, y = mean, ymin = mean_se_low, ymax = mean_se_high , color=program),
                  data = numSpPerOG_sum,
                  position = position_nudge(x=0.25)) +
  scale_color_manual(values = c("#08519C","#2171B5","#4292C6","#807DBA","#6A51A3","#54278F","#3F007D")) +
  labs(x = "program", y = "number of species per orthogroup") +
  theme_bw()

print(p2)
dev.off()

##############################################################################################
# full - statistics (Kruskal-Wallis, non-parametric) for numbers of species in an orthogroup #
##############################################################################################
kwtest_out <- paste0(work_dir,"comparison/230928_combine/full_numSpPerOG_kwTest_230928.txt")
sink(kwtest_out)
kruskal.test(numPerOG ~ program, data = merged5)
pairwise.wilcox.test(merged5$numPerOG, merged5$program,
                     p.adjust.method = "fdr")

merged5_noON <- merged5[!(merged5$program=="ON"),]
kruskal.test(numPerOG ~ program, data = merged5_noON)
pairwise.wilcox.test(merged5_noON$numPerOG, merged5_noON$program,
                     p.adjust.method = "fdr")

sink()
###################################################################################
# full - Rearrange bar plots for number of genes/species in OG - group by species #
###################################################################################

geneSpOG_pdf <- paste0(work_dir,"comparison/230928_combine/full_numGenesPerSpOG_230928.pdf")
pdf(geneSpOG_pdf)
listOfSpecies <- c("Aar", "Ath", "Bra", "Chi", "Cru", "Csa", "Sal", "Tar")
for (species2 in listOfSpecies){
  targetSp <- subset(merged4, species...1==species2)
  
  numGenesPerSp <- c(0:11)
  targetSp %>%
    select("BR", "OF_blast", "OF_diamond", "OF_mmseqs", "SP_diamond", "SP_mmseqs", "ON") %>%
    gather(key = program, value = numGenes) %>% 
    group_by(program) %>%
    mutate(numGenesSp = numGenesPerSp[numGenesPerSp])-> geneMatrix1
  
  geneMatrix1$program <- factor(geneMatrix1$program,
                                levels = c("BR", "OF_blast", "OF_diamond", "OF_mmseqs", 
                                           "SP_diamond", "SP_mmseqs", "ON"))
  print(
    ggplot(geneMatrix1, aes(fill=factor(numGenesSp), y=numGenes, x=program)) + 
      geom_bar(position="stack", stat="identity") +
      scale_fill_brewer(palette="PRGn") +
      ggtitle(paste0("number of genes per species in OG:", species2)) +
      xlab("Program") +
      ylab("number of genes per species per OG") +
      theme_bw()
  )
  
  print(
    ggplot(geneMatrix1, aes(fill=factor(numGenesSp), y=numGenes, x=program)) + 
      geom_bar(position="fill", stat="identity") +
      scale_fill_brewer(palette="PRGn") +
      ggtitle(paste0("number of genes per species in OG:", species2)) +
      xlab("Program") +
      ylab("number of genes per species per OG") +
      theme_bw()
  )
}

dev.off()
###########################################################################################
# full - statistics (Kruskal-Wallis) for number of genes/species in OG - group by species #
###########################################################################################

#stats - Aar
merged6 %>%
  select("program", "Aar") %>%
  gather(key = program, value = Aar) %>% 
  group_by(program) %>%
  count(Aar) -> sumMatrix_aar

ktest_aar <- kruskal.test(Aar ~ program, data = merged6)
wilcox_aar <- pairwise.wilcox.test(merged6$Aar, merged6$program, p.adjust.method = "fdr")
wilcox_aar_df <- as.data.frame(wilcox_aar$p.value)

#stats - Ath
merged6 %>%
  select("program", "Ath") %>%
  gather(key = program, value = Ath) %>% 
  group_by(program) %>%
  count(Ath) -> sumMatrix_ath

ktest_ath <- kruskal.test(Ath ~ program, data = merged6)
wilcox_ath <- pairwise.wilcox.test(merged6$Ath, merged6$program, p.adjust.method = "fdr")
wilcox_ath_df <- as.data.frame(wilcox_ath$p.value)

#stats - Bra
merged6 %>%
  select("program", "Bra") %>%
  gather(key = program, value = Bra) %>% 
  group_by(program) %>%
  count(Bra) -> sumMatrix_bra

ktest_bra <- kruskal.test(Bra ~ program, data = merged6)
wilcox_bra <- pairwise.wilcox.test(merged6$Bra, merged6$program, p.adjust.method = "fdr")
wilcox_bra_df <- as.data.frame(wilcox_bra$p.value)

#stats - Chi
merged6 %>%
  select("program", "Chi") %>%
  gather(key = program, value = Chi) %>% 
  group_by(program) %>%
  count(Chi) -> sumMatrix_chi

ktest_chi <- kruskal.test(Chi ~ program, data = merged6)
wilcox_chi <- pairwise.wilcox.test(merged6$Chi, merged6$program, p.adjust.method = "fdr")
wilcox_chi_df <- as.data.frame(wilcox_chi$p.value)

#stats - Cru
merged6 %>%
  select("program", "Cru") %>%
  gather(key = program, value = Cru) %>% 
  group_by(program) %>%
  count(Cru) -> sumMatrix_cru

ktest_cru <- kruskal.test(Cru ~ program, data = merged6)
wilcox_cru <- pairwise.wilcox.test(merged6$Cru, merged6$program, p.adjust.method = "fdr")
wilcox_cru_df <- as.data.frame(wilcox_cru$p.value)

#stats - Csa
merged6 %>%
  select("program", "Csa") %>%
  gather(key = program, value = Csa) %>% 
  group_by(program) %>%
  count(Csa) -> sumMatrix_csa

ktest_csa <- kruskal.test(Csa ~ program, data = merged6)
wilcox_csa <- pairwise.wilcox.test(merged6$Csa, merged6$program, p.adjust.method = "fdr")
wilcox_csa_df <- as.data.frame(wilcox_csa$p.value)

#stats - Sal
merged6 %>%
  select("program", "Sal") %>%
  gather(key = program, value = Sal) %>% 
  group_by(program) %>%
  count(Sal) -> sumMatrix_sal

ktest_sal <- kruskal.test(Sal ~ program, data = merged6)
wilcox_sal <- pairwise.wilcox.test(merged6$Sal, merged6$program, p.adjust.method = "fdr")
wilcox_sal_df <- as.data.frame(wilcox_sal$p.value)

#stats - Tar
merged6 %>%
  select("program", "Tar") %>%
  gather(key = program, value = Tar) %>% 
  group_by(program) %>%
  count(Tar) -> sumMatrix_tar

ktest_tar <- kruskal.test(Tar ~ program, data = merged6)
wilcox_tar <- pairwise.wilcox.test(merged6$Tar, merged6$program, p.adjust.method = "fdr")
wilcox_tar_df <- as.data.frame(wilcox_tar$p.value)

geneSpOG <- c(sumMatrix_aar, sumMatrix_ath, sumMatrix_bra, sumMatrix_chi, sumMatrix_cru,
              sumMatrix_csa, sumMatrix_sal, sumMatrix_tar)
geneSpOG_merge <- bind_cols(geneSpOG)
wilcox_list <- c(wilcox_aar_df, wilcox_ath_df, wilcox_bra_df, wilcox_chi_df, 
                 wilcox_cru_df, wilcox_csa_df, wilcox_sal_df, wilcox_tar_df)
wilcox_all <- bind_cols(wilcox_list)

geneSpOG_out <- paste0(work_dir,"comparison/230928_combine/full_numGenesPerSpOG_sumStats_230928.csv")
write.csv(geneSpOG_merge, geneSpOG_out)
wilcox_out <- paste0(work_dir,"comparison/230928_combine/full_numGenesPerSpOG_wilcox_230928.csv")
write.csv(wilcox_all, wilcox_out)

kwtest_out <- paste0(work_dir,"comparison/230928_combine/full_numGenesPerSpOG_kwTest_230928.txt")
sink(kwtest_out)
print(ktest_aar)
print(ktest_ath)
print(ktest_bra)
print(ktest_chi)
print(ktest_cru)
print(ktest_csa)
print(ktest_sal)
print(ktest_tar)
sink()

#################################################################################
# orthology counts - comparing programs with a baseline search of species pairs #
#################################################################################

# hand-formatted table from information in TableS7
ortho_path <- paste0(work_dir, "orthologCounts/combinedOrthologCount_reformat_230928.csv")
orthoCount <- read.csv(ortho_path)

plots_orthoCount <- paste0(work_dir, "comparison/230928_combine/orthologyCountsPlots_230928.pdf")
pdf(plots_orthoCount)

diploid_list <- c("AthChi", "AthCru","AthTar", "AthAar",
                  "ChiCru", "ChiTar", "ChiAar", 
                  "CruTar", "CruAar", "TarAar")

full_list <- c("AthChi","AthCru","AthTar","AthSal","AthBra","AthCsa","AthAar","ChiCru",
               "ChiTar","ChiSal","ChiBra","ChiCsa","ChiAar","CruTar","CruSal","CruBra",
               "CruCsa","CruAar","TarSal","TarBra","TarCsa","TarAar","SalBra","SalCsa",
               "SalAar","BraCsa","BraAar","CsaAar")

for (speciesPair in diploid_list){
  subset_df1 <- subset(orthoCount, ploidy=="diploid" & SpPair==speciesPair)
  subset_df1$program <- factor(subset_df1$program,
                              levels = c("OFb_base", "BR", "OFb", "OFd", "OFm",
                                         "SPd", "SPm", "ON"))
  p1 <- ggplot(subset_df1, aes(fill=category, y=n, x=program)) + 
    geom_bar(position="stack", stat="identity") +
    scale_fill_brewer(palette="RdPu")+
    ggtitle(paste0("number of counts in categories, diploid: ", speciesPair)) +
    xlab("Program") +
    ylab("number in each category") +
    theme_bw()
  p2 <- ggplot(subset_df1, aes(fill=category, y=n, x=program)) + 
    geom_bar(position="fill", stat="identity") +
    scale_fill_brewer(palette="RdPu")+
    ggtitle(paste0("number of counts in categories, diploid: ", speciesPair)) +
    xlab("Program") +
    ylab("proportion of the number in each category") +
    theme_bw()
  
  print(p1)
  print(p2)
}

for (speciesPair in full_list){
  subset_df <- subset(orthoCount, SpPair==speciesPair & ploidy=="full")
  subset_df$program <- factor(subset_df$program,
                              levels = c("OFb_base", "BR", "OFb", "OFd", "OFm",
                                         "SPd", "SPm", "ON"))
  p1 <- ggplot(subset_df, aes(fill=category, y=n, x=program)) + 
    geom_bar(position="stack", stat="identity") +
    scale_fill_brewer(palette="OrRd")+
    ggtitle(paste0("number of counts in categories, full: ", speciesPair)) +
    xlab("Program") +
    ylab("number in each category") +
    theme_bw()
  p2 <- ggplot(subset_df, aes(fill=category, y=n, x=program)) + 
    geom_bar(position="fill", stat="identity") +
    scale_fill_brewer(palette="OrRd")+
    ggtitle(paste0("number of counts in categories, full: ", speciesPair)) +
    xlab("Program") +
    ylab("proportion of the number in each category") +
    theme_bw()
  
  print(p1)
  print(p2)
}
dev.off()
#############################################################
# pairwise orthogroup composition comparisons - new metrics #
#############################################################
metrics_dip_path <- paste0(work_dir, "comparison/230928_diploid/diploid_pairwiseMetrics_231003.csv")
metrics_mix_path <- paste0(work_dir, "comparison/230928_full/full_pairwiseMetrics_231004.csv")

metrics_dip <- read.csv(metrics_dip_path)
metrics_mix <- read.csv(metrics_mix_path)

## summary stats - diploid ##                                     
metrics_dip %>%
  group_by(program1, program2) %>% 
  filter(!grepl('na', RS)) %>%
  mutate_at("RS", as.numeric) %>%
  count(RS) %>%
  mutate(sum = sum(n), prop = n/sum(n)) %>%
  filter(RS==1) -> metrics_dip_R1

metrics_dip %>%
  group_by(program1, program2) %>%
  filter(!grepl('na', RS)) %>%
  mutate_at("RS", as.numeric) %>%
  summarise(
    n = n(),
    mean = mean(RS),
    median = median(RS),
    sd   = sd(RS),
    se = sd/sqrt(n),
    mean_se_high = mean + se,
    mean_se_low = mean - se) -> metrics_dip_R2

metrics_dip %>%
  group_by(program1, program2) %>%
  filter(!grepl('na', ARS)) %>%
  mutate_at("ARS", as.numeric) %>%
  count(ARS) %>%
  mutate(sum = sum(n), prop = n/sum(n)) %>%
  filter(ARS==1) -> metrics_dip_A1

metrics_dip %>%
  group_by(program1, program2) %>%
  filter(!grepl('na', ARS)) %>%
  mutate_at("ARS", as.numeric) %>%
  summarise(
    n = n(),
    mean = mean(ARS),
    median = median(ARS),
    sd   = sd(ARS),
    se = sd/sqrt(n),
    mean_se_high = mean + se,
    mean_se_low = mean - se) -> metrics_dip_A2

metrics_dip %>%
  group_by(program1, program2) %>%
  filter(!grepl('na', Jaccard)) %>%
  mutate_at("Jaccard", as.numeric) %>%
  count(Jaccard) %>%
  mutate(sum = sum(n), prop = n/sum(n)) %>%
  filter(Jaccard==1) -> metrics_dip_J1

metrics_dip %>%
  group_by(program1, program2) %>%
  filter(!grepl('na', Jaccard)) %>%
  mutate_at("Jaccard", as.numeric) %>%
  summarise(
    n = n(),
    mean = mean(Jaccard),
    median = median(Jaccard),
    sd   = sd(Jaccard),
    se = sd/sqrt(n),
    mean_se_high = mean + se,
    mean_se_low = mean - se) -> metrics_dip_J2
 
## summary stats - mix ##
metrics_mix %>%
  group_by(program1, program2) %>%
  filter(!grepl('na', RS)) %>%
  mutate_at("RS", as.numeric) %>%
  count(RS) %>%
  mutate(sum = sum(n), prop = n/sum(n)) %>%
  filter(RS==1) -> metrics_mix_R1

metrics_mix %>%
  group_by(program1, program2) %>%
  filter(!grepl('na', RS)) %>%
  mutate_at("RS", as.numeric) %>%
  summarise(
    n = n(),
    mean = mean(RS),
    median = median(RS),
    sd   = sd(RS),
    se = sd/sqrt(n),
    mean_se_high = mean + se,
    mean_se_low = mean - se) -> metrics_mix_R2

metrics_mix %>%
  group_by(program1, program2) %>%
  filter(!grepl('na', ARS)) %>%
  mutate_at("ARS", as.numeric) %>%
  count(ARS) %>%
  mutate(sum = sum(n), prop = n/sum(n)) %>%
  filter(ARS==1) -> metrics_mix_A1

metrics_mix %>%
  group_by(program1, program2) %>%
  filter(!grepl('na', ARS)) %>%
  mutate_at("ARS", as.numeric) %>%
  summarise(
    n = n(),
    mean = mean(ARS),
    median = median(ARS),
    sd   = sd(ARS),
    se = sd/sqrt(n),
    mean_se_high = mean + se,
    mean_se_low = mean - se) -> metrics_mix_A2

metrics_mix %>%
  group_by(program1, program2) %>%
  filter(!grepl('na', Jaccard)) %>%
  mutate_at("Jaccard", as.numeric) %>%
  count(Jaccard) %>%
  mutate(sum = sum(n), prop = n/sum(n)) %>%
  filter(Jaccard==1) -> metrics_mix_J1

metrics_mix %>%
  group_by(program1, program2) %>%
  filter(!grepl('na', Jaccard)) %>%
  mutate_at("Jaccard", as.numeric) %>%
  summarise(
    n = n(),
    mean = mean(Jaccard),
    median = median(Jaccard),
    sd   = sd(Jaccard),
    se = sd/sqrt(n),
    mean_se_high = mean + se,
    mean_se_low = mean - se) -> metrics_mix_J2

# Create a blank workbook
metrics_wb <- createWorkbook()
# Add some sheets to the workbook
addWorksheet(metrics_wb, "dipRS1")
addWorksheet(metrics_wb, "dipRS2")
addWorksheet(metrics_wb, "dipARS1")
addWorksheet(metrics_wb, "dipARS2")
addWorksheet(metrics_wb, "dipJI1")
addWorksheet(metrics_wb, "dipJI2")
addWorksheet(metrics_wb, "mixRS1")
addWorksheet(metrics_wb, "mixRS2")
addWorksheet(metrics_wb, "mixARS1")
addWorksheet(metrics_wb, "mixARS2")
addWorksheet(metrics_wb, "mixJI1")
addWorksheet(metrics_wb, "mixJI2")
# Write the data to the sheets
writeData(metrics_wb, sheet = "dipRS1", x = metrics_dip_R1)
writeData(metrics_wb, sheet = "dipRS2", x = metrics_dip_R2)
writeData(metrics_wb, sheet = "dipARS1", x = metrics_dip_A1)
writeData(metrics_wb, sheet = "dipARS2", x = metrics_dip_A2)
writeData(metrics_wb, sheet = "dipJI1", x = metrics_dip_J1)
writeData(metrics_wb, sheet = "dipJI2", x = metrics_dip_J2)
writeData(metrics_wb, sheet = "mixRS1", x = metrics_mix_R1)
writeData(metrics_wb, sheet = "mixRS2", x = metrics_mix_R2)
writeData(metrics_wb, sheet = "mixARS1", x = metrics_mix_A1)
writeData(metrics_wb, sheet = "mixARS2", x = metrics_mix_A2)
writeData(metrics_wb, sheet = "mixJI1", x = metrics_mix_J1)
writeData(metrics_wb, sheet = "mixJI2", x = metrics_mix_J2)
# Export the file
metrics_out <- paste0(work_dir,"comparison/230928_combine/metricsSummary_DipMix_231004.xlsx")
saveWorkbook(metrics_wb, metrics_out)


## heatmap - diploid ##
metrics_dip_R1$program2 <- factor(metrics_dip_R1$program2,
                                  levels = c("BR", "OFb", "OFd", "OFm", "SPd", "SPm", "ON"))
m1 <- ggplot(metrics_dip_R1, aes(x=program1, y=program2, fill= prop)) + 
  geom_tile() +
  scale_fill_distiller(palette = "Reds", trans = "reverse", limits = c(1, 0)) +
  geom_text(aes(label = paste(n, paste("(",round(prop, digits=3),")",sep=""),
                              sep = "\n"))) + 
  ggtitle("number of orthogroups that are exactly the same, diploid - RS") +
  xlab("method 1") +
  ylab("method 2") +
  theme_bw()

metrics_dip_R2$program1 <- factor(metrics_dip_R2$program1,
                                  levels = c("BR", "OFb", "OFd", "OFm", "SPd", "SPm", "ON"))
metrics_dip_R2$program2 <- factor(metrics_dip_R2$program2,
                                  levels = c("BR", "OFb", "OFd", "OFm", "SPd", "SPm", "ON"))
m2 <- ggplot(metrics_dip_R2, aes(x=program2, y=program1, fill= mean)) + 
  geom_tile() +
  scale_fill_distiller(palette = "Greys", trans = "reverse", limits = c(1, 0)) +
  geom_text(aes(label = round(mean, digits=3)), color = "white") + 
  ggtitle("mean RS value, diploid") +
  xlab("method 1") +
  ylab("method 2") +
  theme_bw()

metrics_dip_A1$program2 <- factor(metrics_dip_A1$program2,
  levels = c("BR", "OFb", "OFd", "OFm", "SPd", "SPm", "ON"))
m3 <- ggplot(metrics_dip_A1, aes(x=program1, y=program2, fill= prop)) + 
        geom_tile() +
        scale_fill_distiller(palette = "Reds", trans = "reverse", limits = c(1, 0)) +
        geom_text(aes(label = paste(n, paste("(",round(prop, digits=3),")",sep=""),
                                    sep = "\n"))) + 
        ggtitle("number of orthogroups that are exactly the same, diploid - ARS") +
        xlab("method 1") +
        ylab("method 2") +
        theme_bw()

metrics_dip_A2$program1 <- factor(metrics_dip_A2$program1,
                                  levels = c("BR", "OFb", "OFd", "OFm", "SPd", "SPm", "ON"))
metrics_dip_A2$program2 <- factor(metrics_dip_A2$program2,
                                  levels = c("BR", "OFb", "OFd", "OFm", "SPd", "SPm", "ON"))
m4 <- ggplot(metrics_dip_A2, aes(x=program2, y=program1, fill= mean)) + 
        geom_tile() +
        scale_fill_distiller(palette = "Greys", trans = "reverse", limits = c(1, 0)) +
        geom_text(aes(label = round(mean, digits=3)), color = "white") + 
        ggtitle("mean ARS value, diploid") +
        xlab("method 1") +
        ylab("method 2") +
        theme_bw()

metrics_dip_J1$program2 <- factor(metrics_dip_J1$program2,
                                  levels = c("BR", "OFb", "OFd", "OFm", "SPd", "SPm", "ON"))
m5 <- ggplot(metrics_dip_J1, aes(x=program1, y=program2, fill= prop)) + 
  geom_tile() +
  scale_fill_distiller(palette = "Reds", trans = "reverse", limits = c(1, 0)) +
  geom_text(aes(label = paste(n, paste("(",round(prop, digits=3),")",sep=""),
                              sep = "\n"))) + 
  ggtitle("number of orthogroups that are exactly the same, diploid - JI") +
  xlab("method 1") +
  ylab("method 2") +
  theme_bw()

metrics_dip_J2$program1 <- factor(metrics_dip_J2$program1,
                                  levels = c("BR", "OFb", "OFd", "OFm", "SPd", "SPm", "ON"))
metrics_dip_J2$program2 <- factor(metrics_dip_J2$program2,
                                  levels = c("BR", "OFb", "OFd", "OFm", "SPd", "SPm", "ON"))
m6 <- ggplot(metrics_dip_J2, aes(x=program2, y=program1, fill= mean)) + 
  geom_tile() +
  scale_fill_distiller(palette = "Greys", trans = "reverse", limits = c(1, 0)) +
  geom_text(aes(label = round(mean, digits=3)), color = "white") +  
  ggtitle("mean JI value, diploid") +
  xlab("method 1") +
  ylab("method 2") +
  theme_bw()

## heatmap - mixed ##
metrics_mix_R1$program2 <- factor(metrics_mix_R1$program2,
                                  levels = c("BR", "OFb", "OFd", "OFm", "SPd", "SPm", "ON"))
m7 <- ggplot(metrics_mix_R1, aes(x=program1, y=program2, fill= prop)) + 
  geom_tile() +
  scale_fill_distiller(palette = "Reds", trans = "reverse", limits = c(1, 0)) +
  geom_text(aes(label = paste(n, paste("(",round(prop, digits=3),")",sep=""),
                              sep = "\n"))) + 
  ggtitle("number of orthogroups that are exactly the same, mixed - RS") +
  xlab("method 1") +
  ylab("method 2") +
  theme_bw()

metrics_mix_R2$program1 <- factor(metrics_mix_R2$program1,
                                  levels = c("BR", "OFb", "OFd", "OFm", "SPd", "SPm", "ON"))
metrics_mix_R2$program2 <- factor(metrics_mix_R2$program2,
                                  levels = c("BR", "OFb", "OFd", "OFm", "SPd", "SPm", "ON"))
m8 <- ggplot(metrics_mix_R2, aes(x=program2, y=program1, fill= mean)) + 
  geom_tile() +
  scale_fill_distiller(palette = "Greys", trans = "reverse", limits = c(1, 0)) +
  geom_text(aes(label = round(mean, digits=3)), color = "white") + 
  ggtitle("mean RS value, mixed") +
  xlab("method 1") +
  ylab("method 2") +
  theme_bw()

metrics_mix_A1$program2 <- factor(metrics_mix_A1$program2,
                                  levels = c("BR", "OFb", "OFd", "OFm", "SPd", "SPm", "ON"))
m9 <- ggplot(metrics_mix_A1, aes(x=program1, y=program2, fill= prop)) + 
  geom_tile() +
  scale_fill_distiller(palette = "Reds", trans = "reverse", limits = c(1, 0)) +
  geom_text(aes(label = paste(n, paste("(",round(prop, digits=3),")",sep=""),
                              sep = "\n"))) + 
  ggtitle("number of orthogroups that are exactly the same, mixed - ARS") +
  xlab("method 1") +
  ylab("method 2") +
  theme_bw()

metrics_mix_A2$program1 <- factor(metrics_mix_A2$program1,
                                  levels = c("BR", "OFb", "OFd", "OFm", "SPd", "SPm", "ON"))
metrics_mix_A2$program2 <- factor(metrics_mix_A2$program2,
                                  levels = c("BR", "OFb", "OFd", "OFm", "SPd", "SPm", "ON"))
m10 <- ggplot(metrics_mix_A2, aes(x=program2, y=program1, fill= mean)) + 
  geom_tile() +
  scale_fill_distiller(palette = "Greys", trans = "reverse", limits = c(1, 0)) +
  geom_text(aes(label = round(mean, digits=3)), color = "white") + 
  ggtitle("mean ARS value, mixed") +
  xlab("method 1") +
  ylab("method 2") +
  theme_bw()

metrics_mix_J1$program2 <- factor(metrics_mix_J1$program2,
                                  levels = c("BR", "OFb", "OFd", "OFm", "SPd", "SPm", "ON"))
m11 <- ggplot(metrics_mix_J1, aes(x=program1, y=program2, fill= prop)) + 
  geom_tile() +
  scale_fill_distiller(palette = "Reds", trans = "reverse", limits = c(1, 0)) +
  geom_text(aes(label = paste(n, paste("(",round(prop, digits=3),")",sep=""),
                          sep = "\n"))) +
  ggtitle("number of orthogroups that are exactly the same, mixed - JI") +
  xlab("method 1") +
  ylab("method 2") +
  theme_bw()

metrics_mix_J2$program1 <- factor(metrics_mix_J2$program1,
                                  levels = c("BR", "OFb", "OFd", "OFm", "SPd", "SPm", "ON"))
metrics_mix_J2$program2 <- factor(metrics_mix_J2$program2,
                                  levels = c("BR", "OFb", "OFd", "OFm", "SPd", "SPm", "ON"))
m12 <- ggplot(metrics_mix_J2, aes(x=program2, y=program1, fill= mean)) + 
  geom_tile() +
  scale_fill_distiller(palette = "Greys", trans = "reverse", limits = c(1, 0)) +
  geom_text(aes(label = round(mean, digits=3)), color = "white") + 
  ggtitle("mean JI value, mixed") +
  xlab("method 1") +
  ylab("method 2") +
  theme_bw()

## diploid histogram ##
metrics_dip %>%
  group_by(program1, program2) %>% 
  filter(!grepl('na', Jaccard)) %>%
  mutate_at("Jaccard", as.numeric) -> metrics_dip_J3

metrics_dip_J3$program1 <- factor(metrics_dip_J3$program1,
                               levels = c("BR", "OFb", "OFd", "OFm", "SPd", "SPm", "ON"))
metrics_dip_J3$program2 <- factor(metrics_dip_J3$program2,
                               levels = c("BR", "OFb", "OFd", "OFm", "SPd", "SPm", "ON"))
h1 <- metrics_dip_J3 %>%
  ggplot( aes(x=Jaccard)) +
  geom_bar() +
  scale_x_binned() +
  facet_grid(program1~program2) +
  theme_bw()
print(h1)

metrics_dip %>%
  group_by(program1, program2) %>% 
  filter(!grepl('na', RS)) %>%
  mutate_at("RS", as.numeric) -> metrics_dip_R3

metrics_dip_R3$program1 <- factor(metrics_dip_R3$program1,
                                  levels = c("BR", "OFb", "OFd", "OFm", "SPd", "SPm", "ON"))
metrics_dip_R3$program2 <- factor(metrics_dip_R3$program2,
                                  levels = c("BR", "OFb", "OFd", "OFm", "SPd", "SPm", "ON"))

h2 <- metrics_dip_R3 %>%
  ggplot( aes(x=RS)) +
  geom_bar() +
  scale_x_binned() +
  facet_grid(program1~program2) +
  theme_bw()
print(h2)

## mix histogram ##
metrics_mix %>%
  group_by(program1, program2) %>% 
  filter(!grepl('na', Jaccard)) %>%
  mutate_at("Jaccard", as.numeric) -> metrics_mix_J3

metrics_mix_J3$program1 <- factor(metrics_mix_J3$program1,
                               levels = c("BR", "OFb", "OFd", "OFm", "SPd", "SPm", "ON"))
metrics_mix_J3$program2 <- factor(metrics_mix_J3$program2,
                               levels = c("BR", "OFb", "OFd", "OFm", "SPd", "SPm", "ON"))
h3 <- metrics_mix_J3 %>%
  ggplot( aes(x=Jaccard)) +
  geom_bar() +
  scale_x_binned() +
  facet_grid(program1~program2) +
  theme_bw()
print(h3)

metrics_mix %>%
  group_by(program1, program2) %>% 
  filter(!grepl('na', RS)) %>%
  mutate_at("RS", as.numeric) -> metrics_mix_R3

metrics_mix_R3$program1 <- factor(metrics_mix_R3$program1,
                                  levels = c("BR", "OFb", "OFd", "OFm", "SPd", "SPm", "ON"))
metrics_mix_R3$program2 <- factor(metrics_mix_R3$program2,
                                  levels = c("BR", "OFb", "OFd", "OFm", "SPd", "SPm", "ON"))

h4 <- metrics_mix_R3 %>%
  ggplot( aes(x=RS)) +
  geom_bar() +
  scale_x_binned() +
  facet_grid(program1~program2) +
  theme_bw()
print(h4)

## print plots ##
plots_metrics_dip <- paste0(work_dir, "comparison/230928_combine/metricsSummaryDiploid_231004.pdf")
pdf(plots_metrics_dip)
print(m1)
print(m2)
print(m3)
print(m4)
print(m5)
print(m6)
print(h1)
print(h2)
dev.off()

plots_metrics_mix <- paste0(work_dir, "comparison/230928_combine/metricsSummaryMixed_231004.pdf")
pdf(plots_metrics_mix)
print(m7)
print(m8)
print(m9)
print(m10)
print(m11)
print(m12)
print(h3)
print(h4)
dev.off()


####################################################################
# comparison metrics between species pairs against baseline method #
####################################################################
spPair_dip_path <- paste0(work_dir, "comparison/230928_diploid/diploid_speciesPairsJI_230929.csv")
spPair_mix_path <- paste0(work_dir, "comparison/230928_full/mixed_speciesPairsJI_230929.csv")

spPair_dip <- read.csv(spPair_dip_path)
spPair_mix <- read.csv(spPair_mix_path)

## heatmap - diploid ##
spPair_dip$method <- factor(spPair_dip$method,
                            levels = c("BR", "OFb", "OFd", "OFm", "SPd", "SPm", "ON"))
sp1 <- ggplot(spPair_dip, aes(x=method, y=SpeciesPair, fill= proportionOfExactOG)) + 
  geom_tile() +
  scale_fill_distiller(palette = "Blues", trans = "reverse", limits = c(1, 0)) +
  geom_text(aes(label = numberOfExactOG)) + 
  ggtitle("number of orthogroups that are exactly the same between method and baseline, diploid - JI") +
  xlab("method") +
  ylab("species pair") +
  theme_bw()

sp1_1 <- ggplot(spPair_dip, aes(x=method, y=SpeciesPair, fill= proportionOfExactOG)) + 
  geom_tile() +
  scale_fill_distiller(palette = "Blues", trans = "reverse", limits = c(1, 0)) +
  geom_text(aes(label = round(proportionOfExactOG, digits=3))) + 
  ggtitle("number of orthogroups that are exactly the same between method and baseline, diploid - JI") +
  xlab("method") +
  ylab("species pair") +
  theme_bw()

sp2 <- ggplot(spPair_dip, aes(x=method, y=SpeciesPair, fill= meanJI)) + 
  geom_tile() +
  scale_fill_distiller(palette = "Greens", trans = "reverse", limits = c(1, 0)) +
  geom_text(aes(label = round(meanJI,digits=3))) + 
  ggtitle("mean value, diploid - JI - between method and baseline") +
  xlab("method") +
  ylab("species pair") +
  theme_bw()

## heatmap - mixed ##
spPair_mix$method <- factor(spPair_mix$method,
                            levels = c("BR", "OFb", "OFd", "OFm", "SPd", "SPm", "ON"))
sp3 <- ggplot(spPair_mix, aes(x=method, y=SpeciesPair, fill= proportionOfExactOG)) + 
  geom_tile() +
  scale_fill_distiller(palette = "Blues", trans = "reverse", limits = c(1, 0)) +
  geom_text(aes(label = numberOfExactOG)) + 
  ggtitle("number of orthogroups that are exactly the same between method and baseline, mixed - JI") +
  xlab("method") +
  ylab("species pair") +
  theme_bw()

sp4 <- ggplot(spPair_mix, aes(x=method, y=SpeciesPair, fill= proportionOfExactOG)) + 
  geom_tile() +
  scale_fill_distiller(palette = "Blues", trans = "reverse", limits = c(1, 0)) +
  geom_text(aes(label = round(proportionOfExactOG, digits=3))) + 
  ggtitle("proportion of orthogroups that are exactly the same between method and baseline, mixed - JI") +
  xlab("method") +
  ylab("species pair") +
  theme_bw()

sp5 <- ggplot(spPair_mix, aes(x=method, y=SpeciesPair, fill= meanJI)) + 
  geom_tile() +
  scale_fill_distiller(palette = "Greens", trans = "reverse", limits = c(1, 0)) +
  geom_text(aes(label = round(meanJI,digits=3))) + 
  ggtitle("mean value, mixed - JI - between method and baseline") +
  xlab("method") +
  ylab("species pair") +
  theme_bw()

plots_speciesPair_dip <- paste0(work_dir, "comparison/230928_combine/heatmap_SpeciesPairs_Diploid_231003.pdf")
pdf(plots_speciesPair_dip)
print(sp1)
print(sp1_1)
print(sp2)
dev.off()
plots_speciesPair_mix <- paste0(work_dir, "comparison/230928_combine/heatmap_SpeciesPairs_Mixed_231003.pdf")
pdf(plots_speciesPair_mix)
print(sp3)
print(sp4)
print(sp5)
dev.off()