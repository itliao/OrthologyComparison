work_dir <-"/path_to_directory/"
setwd (work_dir)

library(ggplot2)
library(ComplexUpset)
library(dplyr)
library(tidyr)
library(tibble)

##########################
# brassicaceae - diploid #
##########################

# function for automatically making figures and saving figures in separate files

fileList_to_graph_d <- function(mainDir) {
  for (i in listDirectory){
    #paths for tables
    category <- i[1]
    program <- i[2]
    
    allCountsPath <- paste0(mainDir, category, "/", program, "_allCounts.tsv")
    binaryPath <- paste0(mainDir, category, "/", program, "_binary.tsv")
    tenMaxPath <- paste0(mainDir, category, "/", program, "_tenMax.tsv")
    
    #read in files
    allCounts <- read.table(allCountsPath, header=T)
    binary <- read.table(binaryPath, header=T)
    tenMax <- read.table(tenMaxPath, header=T)
    
    #UpsetPlot - use allCounts
    
    # first, create a new column called one2one for OG with exactly 1 gene copy per species
    if (program == "orthnet"){
      allCounts <- allCounts %>% 
        add_column(one2one = 
                     if_else(.$Ath == 1 & .$Chi == 1 & .$Cru == 1 & .$Tar == 1 & .$Aar == 1, "1", "0"),
                   .after="Aar")
      
      group <- colnames(allCounts)[c(2:6)]
      allCounts[group] <- allCounts[group] !=0 
      
      print(table(allCounts$one2one))
    } else {
      allCounts <- allCounts %>% 
        add_column(one2one = 
                     if_else(.$Ath == 1 & .$Chi == 1 & .$Cru == 1 & .$Tar == 1 & .$Aar == 1, "1", "0"),
                   .after="Tar")
      
      group <- colnames(allCounts)[c(2:6)]
      allCounts[group] <- allCounts[group] !=0 
      
      print(table(allCounts$one2one))}
    
    upset_out = paste0(mainDir,"/",category, "/figure/", program, "_complexUpset_231017.pdf")
    pdf(upset_out)
    
    # upset plot with ALL interactions
    print(
      upset(allCounts, group, name='OG species composition', width_ratio=0.1)
    )
    
    # upset plot with 10 intersections
    print(
      upset(allCounts, 
            group, 
            name='OG species composition', 
            width_ratio=0.1,
            n_intersections=10,
            min_degree=2,
            max_degree=5)
    )
    
    # upset plot for OG with just 2 species
    print(
      upset(
        allCounts, 
        group, 
        name="OG species composition",
        width_ratio=0.1,
        min_degree=2,
        max_degree=2
      )
    )
    
    # upset plot for OG with just 3 species
    print(
      upset(
        allCounts, 
        group, 
        name="OG species composition",
        width_ratio=0.1,
        min_degree=3,
        max_degree=3
      )
    )
    
    # upset plot for OG with just 3 species, interaction size > 20
    print(
      upset(
        allCounts, 
        group, 
        name="OG species composition",
        width_ratio=0.1,
        min_degree=3,
        max_degree=3,
        min_size=20
      )
    )
    
    # upset plot for OG with just 4 species
    print(
      upset(
        allCounts, 
        group, 
        name="OG species composition",
        width_ratio=0.1,
        min_degree=4,
        max_degree=4
      )
    )
    
    # upset plot for OG with just 5 species
    print(
      upset(
        allCounts,
        group, 
        name="OG species composition",
        width_ratio=0.1,
        min_degree=4,
        max_degree=4,
        min_size = 20
      )
    )
    
    # upset plot for OG with all species with bar colored for first category, interaction size > 20  
    print(
      upset(
        allCounts,
        group,
        name = "OG species composition",
        width_ratio=0.1,
        min_degree=5,
        max_degree=5,
        annotations = list(
          'GeneCopy'=(list(
            aes=aes(x=intersection, fill=one2one),
            geom=list(
              geom_bar(stat='count', position='fill', na.rm=TRUE),
              geom_text(aes(label=!!aes_percentage(relative_to='intersection'),
                            group=one2one,),
                        stat='count',
                        position=position_fill(vjust = .5)),
              scale_y_continuous(labels=scales::percent_format()),
              scale_fill_manual(values=c("1"='#c51b7d', "0"='#4d9221')))
          )
          )
        )
      )
    )
    dev.off()
    
    # regular histogram summary - use binary
    
    binary$numPerOG = apply(binary[,c(2:6)], 1, sum)
    
    if (program == "orthofinder" || program == "orthnet"){
      onePlus = binary[binary$numPerOG != 1, ]
      
      spNum <- 
        onePlus %>% count(numPerOG)
      
      histogram_out = paste0(mainDir, "/",category, "/figure/", program, "_histogram.pdf")
      pdf(histogram_out)
      
      print(
        ggplot(spNum, aes(x=numPerOG, y=n)) + 
          geom_bar(stat = "identity") + 
          geom_text(aes(label = n), vjust = -0.2) +
          ggtitle("Distribution of number of species per OG") +
          xlab("number of species per OG") +
          ylab("Frequency") +
          theme_bw()
      )
      dev.off()
      
    } else {
      spNum <- 
        binary %>% count(numPerOG)
      
      histogram_out = paste0(mainDir, "/",category, "/figure/", program, "_histogram.pdf")
      pdf(histogram_out)
      
      print(
        ggplot(spNum, aes(x=numPerOG, y=n)) + 
          geom_bar(stat = "identity") + 
          geom_text(aes(label = n), vjust = -0.2) +
          ggtitle("Distribution of number of species per OG") +
          xlab("number of species per OG") +
          ylab("Frequency") +
          theme_bw()
      )
      dev.off()
    }
    
    
    # Heatmaps & Stacked bars - use tenMax
    
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
      OGtenMaxOnePlus %>% 
        select("Ath", "Chi", "Cru", "Tar", "Aar") %>%
        gather(key = species , value = numGenes) %>%
        group_by(species) %>%
        count(numGenes) %>%
        ungroup %>%
        spread(species, n) -> sumMatrix1
      
      # graphical format
      OGtenMaxOnePlus %>% 
        select("Ath", "Chi", "Cru", "Tar", "Aar") %>%
        gather(key = species , value = numGenes) %>%
        group_by(species) %>%
        count(numGenes) -> sumMatrix2
      
      heatmap_stackBar_out = paste0(mainDir,"/",category, "/figure/", program, "_heatmapStackBar.pdf")
      pdf(heatmap_stackBar_out)
      
      #stacked barplot
      print(
        ggplot(sumMatrix2, aes(fill=factor(numGenes), y=n, x=species)) + 
          geom_bar(position="stack", stat="identity") +
          scale_fill_brewer(palette="PiYG") +
          ggtitle("number of genes per species in OG - stacked bars") +
          xlab("Species") +
          ylab("number of genes per species per OG") +
          theme_bw()
      )
      
      #heatmap
      print(
        ggplot(sumMatrix2, aes(x=species, y=numGenes, fill= n)) + 
          geom_tile() +
          scale_fill_distiller(palette = "PiYG") +
          geom_text(aes(label = n)) +
          ggtitle("number of genes per species in OG - heatmap") +
          xlab("Species") +
          ylab("number of genes per species per OG") +
          theme_bw()
      )
      
      dev.off()
      
    } else {
      # table format
      tenMax %>% 
        select("Ath", "Chi", "Cru", "Tar", "Aar") %>%
        gather(key = species , value = numGenes) %>%
        group_by(species) %>%
        count(numGenes) %>%
        ungroup %>%
        spread(species, n) -> sumMatrix1
      
      # graphical format
      tenMax %>% 
        select("Ath", "Chi", "Cru", "Tar", "Aar") %>%
        gather(key = species , value = numGenes) %>%
        group_by(species) %>%
        count(numGenes) -> sumMatrix2
      
      heatmap_stackBar_out = paste0(mainDir, "/",category, "/figure/", program, "_heatmapStackBar.pdf")
      pdf(heatmap_stackBar_out)
      #stacked barplot
      print(
        ggplot(sumMatrix2, aes(fill=factor(numGenes), y=n, x=species)) + 
          geom_bar(position="stack", stat="identity") +
          scale_fill_brewer(palette="PiYG") +
          ggtitle("number of genes per species in OG - stacked bars") +
          xlab("Species") +
          ylab("number of genes per species per OG") +
          theme_bw()
      )
      
      #heatmap
      print(
        ggplot(sumMatrix2, aes(x=species, y=numGenes, fill= n)) + 
          geom_tile() +
          scale_fill_distiller(palette = "PiYG") +
          geom_text(aes(label = n)) +
          ggtitle("number of genes per species in OG - heatmap") +
          xlab("Species") +
          ylab("number of genes per species per OG") +
          theme_bw()
      )
      dev.off()
    }
  }
}

# list of the directory names and the programs
# for broccoli - prior to reading in file, need to remove the "#" from "#OG_name"

listDirectory = list(c("BR","broccoli"),
                     c("OF_blast","orthofinder"),
                     c("OF_diamond","orthofinder"),
                     c("OF_mmseqs","orthofinder"),
                     c("ON", "orthnet"),
                     c("SP_diamond","sonicparanoid"),
                     c("SP_mmseqs","sonicparanoid"))

# call function
diploidDir = paste0(work_dir, "diploidBrass/")
fileList_to_graph_d(diploidDir)

####################################
# brassicaceae - higher ploidy set #
####################################

# function for automatically making figures and saving figures in separate files

fileList_to_graph_f <- function(mainDir) {
  for (i in listDirectory){
    #paths for tables
    category <- i[1]
    program <- i[2]
    
    allCountsPath <- paste0(mainDir, category, "/", program, "_allCounts.tsv")
    binaryPath <- paste0(mainDir, category, "/", program, "_binary.tsv")
    tenMaxPath <- paste0(mainDir, category, "/", program, "_tenMax.tsv")
    
    #read in files
    allCounts <- read.table(allCountsPath, header=T)
    binary <- read.table(binaryPath, header=T)
    tenMax <- read.table(tenMaxPath, header=T)
    
    #UpsetPlot - use allCounts
    
    # first, create a new column called one2one for OG with exactly 1 gene copy per species
    allCounts <- allCounts %>% 
      add_column(one2one = if_else(.$Ath == 1 & .$Chi == 1 & .$Cru == 1 & .$Tar == 1 & .$Bra == 1 & .$Csa == 1 & .$Sal == 1 & .$Aar == 1, "1", "0"),
                 .after="Tar")
    
    group <- colnames(allCounts)[c(2:9)]
    allCounts[group] <- allCounts[group] !=0 
    
    print(table(allCounts$one2one))
    
    
    upset_out = paste0(mainDir,"/",category, "/figure/", program, "_complexUpset.pdf")
    pdf(upset_out)
    
    # upset plot with ALL interactions
    print(
      upset(allCounts, group, name='OG species composition', width_ratio=0.1)
    )
    
    # upset plot with 10 intersections
    print(
      upset(allCounts, 
            group, 
            name='OG species composition', 
            width_ratio=0.1, 
            n_intersections=10,
            min_degree=2,
            max_degree=8
            )
    )
    
    # upset plot for OG with just 2 species
    print(
      upset(
        allCounts, 
        group, 
        name="OG species composition",
        width_ratio=0.1,
        min_degree=2,
        max_degree=2
      )
    )
    
    # upset plot for OG with just 3 species
    print(
      upset(
        allCounts, 
        group, 
        name="OG species composition",
        width_ratio=0.1,
        min_degree=3,
        max_degree=3
      )
    )
    
    
    # upset plot for OG with just 4 species
    print(
      upset(
        allCounts, 
        group, 
        name="OG species composition",
        width_ratio=0.1,
        min_degree=4,
        max_degree=4
      )
    )
    
    # upset plot for OG with just 5 species
    print(
      upset(
        allCounts, 
        group, 
        name="OG species composition",
        width_ratio=0.1,
        min_degree=5,
        max_degree=5
      )
    )
   
    # upset plot for OG with just 6 species
    print(
      upset(
        allCounts, 
        group, 
        name="OG species composition",
        width_ratio=0.1,
        min_degree=6,
        max_degree=6
      )
    )
    
    # upset plot for OG with just 7 species
    print(
      upset(
        allCounts, 
        group, 
        name="OG species composition",
        width_ratio=0.1,
        min_degree=7,
        max_degree=7
      )
    )
    
    # upset plot for OG with all 8 species
    print(
      upset(
        allCounts, 
        group, 
        name="OG species composition",
        width_ratio=0.1,
        min_degree=8,
        max_degree=8
      )
    )
    
    # upset plot for OG with all species with bar colored for first category, interaction size > 20  
    print(
      upset(
        allCounts,
        group,
        name = "OG species composition",
        width_ratio=0.1,
        min_degree=8,
        max_degree=8,
        min_size = 20,
        annotations = list(
          'GeneCopy'=(list(
            aes=aes(x=intersection, fill=one2one),
            geom=list(
              geom_bar(stat='count', position='fill', na.rm=TRUE),
              geom_text(aes(label=!!aes_percentage(relative_to='intersection'),
                            group=one2one,),
                        stat='count',
                        position=position_fill(vjust = .5)),
              scale_y_continuous(labels=scales::percent_format()),
              scale_fill_manual(values=c("1"='#c51b7d', "0"='#4d9221')))
          )
          )
        )
      )
    )
    dev.off()
    
    # regular histogram summary - use binary
    
    binary$numPerOG = apply(binary[,c(2:9)], 1, sum)
    
    if (program == "orthofinder" || program == "orthnet"){
      onePlus = binary[binary$numPerOG != 1, ]
      
      spNum <- 
        onePlus %>% count(numPerOG)
      
      histogram_out = paste0(mainDir, "/",category, "/figure/", program, "_histogram.pdf")
      pdf(histogram_out)
      
      print(
        ggplot(spNum, aes(x=numPerOG, y=n)) + 
          geom_bar(stat = "identity") + 
          geom_text(aes(label = n), vjust = -0.2) +
          ggtitle("Distribution of number of species per OG") +
          xlab("number of species per OG") +
          ylab("Frequency") +
          theme_bw()
      )
      dev.off()
      
    } else {
      spNum <- 
        binary %>% count(numPerOG)
      
      histogram_out = paste0(mainDir, "/",category, "/figure/", program, "_histogram.pdf")
      pdf(histogram_out)
      
      print(
        ggplot(spNum, aes(x=numPerOG, y=n)) + 
          geom_bar(stat = "identity") + 
          geom_text(aes(label = n), vjust = -0.2) +
          ggtitle("Distribution of number of species per OG") +
          xlab("number of species per OG") +
          ylab("Frequency") +
          theme_bw()
      )
      dev.off()
    }
    
    
    # Heatmaps & Stacked bars - use tenMax
    
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
      OGtenMaxOnePlus %>% 
        select("Ath", "Chi", "Cru", "Tar", "Bra", "Csa", "Sal", "Aar") %>%
        gather(key = species , value = numGenes) %>%
        group_by(species) %>%
        count(numGenes) %>%
        ungroup %>%
        spread(species, n) -> sumMatrix1
      
      # graphical format
      OGtenMaxOnePlus %>% 
        select("Ath", "Chi", "Cru", "Tar", "Bra", "Csa", "Sal", "Aar") %>%
        gather(key = species , value = numGenes) %>%
        group_by(species) %>%
        count(numGenes) -> sumMatrix2
      
      heatmap_stackBar_out = paste0(mainDir,"/",category, "/figure/", program, "_heatmapStackBar.pdf")
      pdf(heatmap_stackBar_out)
      
      #stacked barplot
      print(
        ggplot(sumMatrix2, aes(fill=factor(numGenes), y=n, x=species)) + 
          geom_bar(position="stack", stat="identity") +
          scale_fill_brewer(palette="PiYG") +
          ggtitle("number of genes per species in OG - stacked bars") +
          xlab("Species") +
          ylab("number of genes per species per OG") +
          theme_bw()
      )
      
      #heatmap
      print(
        ggplot(sumMatrix2, aes(x=species, y=numGenes, fill= n)) + 
          geom_tile() +
          scale_fill_distiller(palette = "PiYG") +
          geom_text(aes(label = n)) +
          ggtitle("number of genes per species in OG - heatmap") +
          xlab("Species") +
          ylab("number of genes per species per OG") +
          theme_bw()
      )
      
      dev.off()
      
    } else {
      # table format
      tenMax %>% 
        select("Ath", "Chi", "Cru", "Tar", "Bra", "Csa", "Sal", "Aar") %>%
        gather(key = species , value = numGenes) %>%
        group_by(species) %>%
        count(numGenes) %>%
        ungroup %>%
        spread(species, n) -> sumMatrix1
      
      # graphical format
      tenMax %>% 
        select("Ath", "Chi", "Cru", "Tar", "Bra", "Csa", "Sal", "Aar") %>%
        gather(key = species , value = numGenes) %>%
        group_by(species) %>%
        count(numGenes) -> sumMatrix2
      
      heatmap_stackBar_out = paste0(mainDir, "/",category, "/figure/", program, "_heatmapStackBar.pdf")
      pdf(heatmap_stackBar_out)
      #stacked barplot
      print(
        ggplot(sumMatrix2, aes(fill=factor(numGenes), y=n, x=species)) + 
          geom_bar(position="stack", stat="identity") +
          scale_fill_brewer(palette="PiYG") +
          ggtitle("number of genes per species in OG - stacked bars") +
          xlab("Species") +
          ylab("number of genes per species per OG") +
          theme_bw()
      )
      
      #heatmap
      print(
        ggplot(sumMatrix2, aes(x=species, y=numGenes, fill= n)) + 
          geom_tile() +
          scale_fill_distiller(palette = "PiYG") +
          geom_text(aes(label = n)) +
          ggtitle("number of genes per species in OG - heatmap") +
          xlab("Species") +
          ylab("number of genes per species per OG") +
          theme_bw()
      )
      dev.off()
    }
  }
}

# list of the directory names and the programs
# for broccoli - prior to reading in file, need to remove the "#" from "#OG_name"

listDirectory = list(c("BR","broccoli"),
                     c("OF_blast","orthofinder"),
                     c("OF_diamond","orthofinder"),
                     c("OF_mmseqs","orthofinder"),
                     c("ON", "orthnet"),
                     c("SP_diamond","sonicparanoid"),
                     c("SP_mmseqs","sonicparanoid"))

# call function
fullDir = paste0(work_dir, "fullBrass/")
fileList_to_graph_f(fullDir)
