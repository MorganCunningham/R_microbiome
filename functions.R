#Script with functions


#Calling Packages

library(tidyverse)
library(vegan)
library(tidyverse)
library(reshape2)
library(knitr)
library(corrplot)
library(Hmisc)
library(glue)
library(indicspecies)
library(xtable)

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.18")
library(BiocManager)
if (!require("ANCOMBC", quietly = TRUE))
  BiocManager::install("ANCOMBC")
library(ANCOMBC)


#install.packages("Polychrome")
library(Polychrome)


set.seed(123)



#makes data ready for analysis
#Inputs(Raw data, taxaLVL, year)

processData = function(rawData, taxaLvl, path){

  
  
  data = rawData %>%
    #filter by bacteria
    mutate(bacteria=ifelse(str_detect(taxLineage, 'Bacteria'), 'yes', 'no' )) %>%
    #filter by taxa rank
    filter(taxRank == taxaLvl, bacteria == "yes") %>%
    #remove columns
    select(-c(taxRank, taxID,depth, taxLineage, bacteria)) %>%
    #rotate dataframe
    column_to_rownames("name")%>% 
    t() %>%
    as.data.frame() %>%
    rownames_to_column("sampleID")%>%
    #create column names from sampleID
    separate(sampleID, c("readType", "sampleID"), sep = "_")%>%
    filter(
        str_detect(readType, "cladeReads")) %>%
    mutate(
        Crop = str_extract(sampleID, 'APP|CAC|CAS|COR|CRA|HBB|LBB|SOY'),
        Environment=str_extract(sampleID, 'e|u'),
        Timepoint=str_extract(sampleID, 't1|t2|t3|t4'),
        Replicate=str_extract(sampleID, '01|02|03|04|05|06|07|08|09|10'),
        Year = getYear(path)
        ) %>%
    #rename factors
    mutate(
        Environment = recode_factor(Environment, "e" = "Exposed", "u" = "Unexposed"),
         Timepoint = recode_factor(Timepoint, "t1" = "T1", "t2" = "T2", "t3" = "T3", "t4" = "T4")) %>%
    #bring everything to front
    select(sampleID, readType, Crop, Environment, Timepoint, Replicate, Year, everything())
  
return(data)
}

#Gets sample year from path
#Input: file fath
#output: Sampling year
getYear = function(path){
  year = str_extract(path, "2020|2021")
  return(year)
  
}

#Input: dataframe from processData
#Output: relative abundance
relAbund = function(processedData){

metaData = processedData %>%
    select(sampleID, readType, Crop, Environment, Timepoint, Replicate, Year)
  
relData = processedData %>%
  select(-c(sampleID, readType, Crop, Environment, Timepoint, Replicate, Year)) %>%
  apply(1,function(x) x/sum(x)) %>% t() %>% as.data.frame() %>%
  cbind(metaData) %>%
  select(sampleID, readType, Crop, Environment, Timepoint, Replicate, Year, everything())

return(relData)
}
    
#Plots relative abundance for crops with 1 timepoint
#Input:processedData
#Uses: relAbund
#output: rel abund plot
relAbundPlot = function(processedData){
  
  
  relData = relAbund(processedData)
  
#define colour pallette of varying colours
  # c25 <- c(
  #   "dodgerblue2", "#E31A1C", # red
  #   "green4",
  #   "#6A3D9A", # purple
  #   "#FF7F00", # orange
  #   "white", "gold1",
  #   "skyblue2", "#FB9A99", # lt pink
  #   "palegreen2",
  #   "#CAB2D6", # lt purple
  #   "#FDBF6F", # lt orange
  #   "gray70", "khaki2",
  #   "maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
  #   "darkturquoise", "green1", "yellow4", "yellow3",
  #   "darkorange4", "brown"
  # )
  
P50 = as.data.frame(createPalette(50,  c("#ff0000", "#00ff00", "#0000ff")))


  
plot =  relData %>%
    pivot_longer(8:ncol(relData), names_to = "Taxa", values_to = "Rel_abund",) %>%
    group_by(Environment, Timepoint, Taxa)%>%
    summarise(mean_abund = mean(Rel_abund))%>%
    ggplot(aes(y = mean_abund*100, x = Environment, fill = Taxa)) +
    geom_bar(stat = 'identity', position = 'stack', colour = "black") +
    #ggtitle("Relative abundance by site and timepoint (2021)") +
    ylab("Relative abundance (%)") +
    xlab("Site")+ 
    #guides(fill=guide_legend(title="Site")) +
    #facet_wrap(~Timepoint, nrow = 1)+
  theme_light()+  
  theme(
      
      plot.title = element_text(size = 20),
      axis.title = element_text(size = 15))+
    scale_fill_manual(values = P50$`createPalette(50, c("#ff0000", "#00ff00", "#0000ff"))`)
  
  return(plot)
}

#Plots relative abundance for crops with 1 timepoint
relAbundPlotTPS = function(dataframe){
  
  #define colour pallette of varying colours
  c25 <- c(
    "dodgerblue2", "#E31A1C", # red
    "green4",
    "#6A3D9A", # purple
    "#FF7F00", # orange
    "white", "gold1",
    "skyblue2", "#FB9A99", # lt pink
    "palegreen2",
    "#CAB2D6", # lt purple
    "#FDBF6F", # lt orange
    "gray70", "khaki2",
    "maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
    "darkturquoise", "green1", "yellow4", "yellow3",
    "darkorange4", "brown"
  )
  
  plot =  dataframe %>%
    pivot_longer(8:ncol(dataframe), names_to = "Taxa", values_to = "Rel_abund",) %>%
    group_by(Environment, Timepoint, Taxa)%>%
    summarise(mean_abund = mean(Rel_abund))%>%
    ggplot(aes(y = mean_abund*100, x = Environment, fill = Taxa)) +
    geom_bar(stat = 'identity', position = 'stack', colour = "black") +
    #ggtitle("Relative abundance by site and timepoint (2021)") +
    ylab("Relative abundance (%)") +
    xlab("Site")+ 
    guides(fill=guide_legend(title="Site")) +
    facet_wrap(~Timepoint, nrow = 1)+
    theme(
      
      plot.title = element_text(size = 20),
      axis.title = element_text(size = 15))+
    scale_fill_manual(values = c25)
  
  return(plot)
}


#Input: Raw data
#output: NMDS scores
getNMDS = function(processedData){

  taxaData = processedData[8:ncol(processedData)]

  #nmds
  nmds = metaMDS(taxaData, distance = "bray")
  
  #dataframe
  nmdsScores = as.data.frame(scores(nmds)$sites) %>%
    mutate(SampleID = processedData$SampleID,
           Timepoint = processedData$Timepoint, 
           Environment = processedData$Environment,
           Replicate = processedData$Replicate,
           Year = processedData$Year)
  
  return(nmdsScores)
  
}

#input: processed Data
#Uses getNMDS function
#output: NMDS plot
plotNMDS = function(processedData){

  nmdsScores = getNMDS(processedData)


nmdsPlot = ggplot(nmdsScores, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(aes(colour = Environment), size = 3.5, alpha = 0.75) + 
  theme(panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        legend.key=element_blank(), axis.text = element_text(colour = "black"),
        axis.title = element_text(colour = "black")) + labs(colour = "Environment") + scale_radius(range = c(1.5,6)) 

return(nmdsPlot)
}


permanovaHelper = function(processedData){
  
  taxaData = processedData %>%
    select(-c(sampleID, readType, Crop, Environment, Timepoint, Replicate, Year))
  
  
  
  permanova = adonis2(taxaData ~ Environment+Replicate, processedData, method = "bray")
  
  return(permanova) 
  
}


runIndicSp = function(processedData){

  
  taxaData = processedData %>%
    select(-c(sampleID, readType, Crop, Environment, Timepoint, Replicate, Year))
  
  indicSpecies = multipatt(taxaData, processedData$Environment, func = "r.g", control = how(nperm=9999))
  
  return(indicSpecies)
  
}

#In: processed data
#Uses run IndicSp
#Out: c(plot,significant sp.)
processResultsIndicSp = function(processedData){
  
  
  relData = relAbund(processedData)
  #Gets results from indicator species analysis
  indicResults = runIndicSp(relData)
  
  #Gets significant taxa in indicator species analysis for later filtering
  significantTaxa = indicResults[["sign"]] %>%
  filter(p.value < 0.05) %>% as.data.frame() %>%
    rownames_to_column("Taxa") %>% 
    rename(Unexposed = s.Unexposed,
           Exposed = s.Exposed)
  
  
  indicPlot = relData %>%
  pivot_longer(8:ncol(relData), names_to = "Taxa", values_to = "Rel_Abund") %>%
  filter(Taxa %in% significantTaxa$Taxa) %>%
  ggplot(aes(y = Rel_Abund*100, x = Environment, fill = Environment)) +
  geom_boxplot() +
  #ggtitle("Environment indicator species analysis - significant species (2021)") +
  ylab("Relative abundance (%)") +
  xlab("Site")+ 
  guides(fill=guide_legend(title="Site")) +
  facet_wrap(~Taxa, scales = "free")+ 
    theme_light() +
    scale_fill_manual(values = c("deepskyblue", "#E7AB4B","#7F636E" , "black"))+
    theme(legend.position= 'none',)
  
  
  return(list(significantTaxa, indicPlot))

  
}

#Input Processed data
#output phyloseq obj for runANCOMBC
createPhyloSeqObj = function(dataProcessed){
  
  taxaData = dataProcessed %>%
    select(-c(sampleID, readType, Crop, Environment, Timepoint, Replicate, Year)) %>%
    otu_table(taxa_are_rows = F)
    
  metaData = dataProcessed %>%
    select(sampleID, readType, Crop, Environment, Timepoint, Replicate, Year) %>%
    sample_data()
  
  phyloSeqObj = phyloseq(taxaData, metaData)
  
  return(phyloSeqObj)
  
}

helperANCOMBC2 = function(phyloSeqObj){
  
  resultsANCOMBC = ancombc2(phyloSeqObj, 
                            fix_formula = "Environment+Replicate",
                            p_adj_method = "BH")
  
  
  return(resultsANCOMBC)
  
}


#Runs ancombc2 by calling ancombc related functions
#input: processed Data
#Uses: createPhyloSeqObj
#output: ancombc2 results
runANCOMBC2 = function(dataProcessed){
  
  
  #create phyloseqobj
  phyloSeqObj = createPhyloSeqObj(dataProcessed)
  
  resultsANCOMBC2 = ancombc2(phyloSeqObj, 
                            fix_formula = "Environment+Replicate",
                            p_adj_method = "BH")
  
  
  return(resultsANCOMBC2)
  
}

#Gets significant results and plots ANCOMBC2 results
processANCOMBC2 = function(processedData){
  
  resultsANCOMBC2 = runANCOMBC2(processedData)
  
  
  resultSheet = resultsANCOMBC2$res %>%
    select(taxon, ends_with("Unexposed"))
  
  
  ANCOMBC2plot = resultsANCOMBC2$res %>%
    select(taxon, ends_with("Unexposed")) %>%
    #filter(p_EnvironmentUnexposed < 0.05) %>%
    mutate(`Direction LFC` = as.factor(ifelse(lfc_EnvironmentUnexposed >0, "Positive LFC", "Negative LFC")),
           Significance = ifelse(q_EnvironmentUnexposed <0.05, "*", "")) %>%
    ggplot(aes(x = reorder(taxon, lfc_EnvironmentUnexposed), y = lfc_EnvironmentUnexposed, fill = `Direction LFC`)) +
    geom_bar(stat = "identity", color = "black", 
             position = position_dodge(width = 0.4)) +
    geom_errorbar(aes(ymin = lfc_EnvironmentUnexposed - se_EnvironmentUnexposed, ymax = lfc_EnvironmentUnexposed + se_EnvironmentUnexposed), 
                  width = 0.2, position = position_dodge(0.05), color = "black")+
    theme_light()+
    theme(plot.title = element_text(hjust = 0.5),
          panel.grid.minor.y = element_blank(),
          axis.text.x = element_text(angle = 60, hjust = 1))+
    xlab("Taxa")+ylab("Log Fold Change")+
    geom_text(aes(label = Significance, y = lfc_EnvironmentUnexposed), colour = "black", size = 15)
  
  
  return(list(resultSheet, ANCOMBC2plot))
}

#Analysis: String from ANCOMBC, IndicSP, AlphaDiversity, BetaDiversity, RelAbund
saveFiles = function(processedData, sheet,plot, analysis, taxaLvl){
  
    #save files
    ##specify main directory
    main_dir = getwd()
    
    #specify the experiment type (crop, eccc)
    sub_dir = paste(unique(processedData$Crop), unique(processedData$Year), taxaLvl, sep = "_")
    
    output_dir <- file.path(main_dir, "Outputs", sub_dir, analysis)
    
    ##check if the directory already exists and create if it doesn't
    if (!dir.exists(output_dir)){
      dir.create(output_dir, recursive = TRUE)
    } else {
      print("Dir already exists!")
    }
    
    #create file names
    plot_name = paste(unique(processedData$Crop), unique(processedData$Year), 
                      taxaLvl, analysis, "plot.png",sep = "_")
    sheet_name = paste(unique(processedData$Crop), unique(processedData$Year), 
                       taxaLvl, analysis, "sheet.csv",sep = "_")
    
    #save files to the directory
  
    
    plot = plot
    filename = plot_name
    path = output_dir
    height = 10
    width = 10
    ggsave(filename = filename, plot = plot, path = path, width = width, height = height, dpi=700)
    
    
    
    write.csv(sheet, paste(output_dir, sheet_name, sep = '/'))
  
  }

#Runs everything
callSave = function(path, taxaLvl){
  
  rawData = read.csv(path, header = TRUE)
  processedData = processData(rawData, taxaLvl, path)
  
  #Save relative abundance
  print("Saving Relative abundance")
  saveFiles(processedData, 
            sheet = relAbund(processedData), 
            plot = relAbundPlot(processedData), 
            "relativeAbundance", 
            taxaLvl)
  
  #Save alpha diversity
  
  
  #Save betadiversity 
  print("Saving Beta Diversity")
  saveFiles(processedData, 
            sheet = permanovaHelper(processedData), 
            plot = plotNMDS(processedData), 
            "BetaDiversity", 
            taxaLvl)
  
  #Save indicator species
  print("Saving indicator species")
  saveFiles(processedData, 
            sheet = processResultsIndicSp(processedData)[[1]], 
            plot = processResultsIndicSp(processedData)[[2]], 
            "IndicatorSpecies", 
            taxaLvl)
  
  #Save ANCOMBC2
  print("Saving ANCOMBC2")
  saveFiles(processedData, 
            sheet = processANCOMBC2(processedData)[[1]], 
            plot = processANCOMBC2(processedData)[[2]], 
            "ANCOMBC2", 
            taxaLvl)
  
  
  #Write Markdown File
  
  
}

path = "~/AAFC WORK/act2_taxonomic_tables/act2_taxonomic_tables/app_2021_aggregated_counts.csv"
taxaLvl = "S"
callSave(path, taxaLvl)



rawData = read.csv(path, header = TRUE)
processedData = processData(rawData, taxaLvl, path)



relData = relAbund(processedData)
relPlot = relAbundPlot(data)

NMDSplot = plotNMDS(data)

permanovaResults = permanovaHelper(data)

indic = processResultsIndicSp(relData)

saveFiles(processedData, 
          sheet = indic[[1]], 
          plot = indic[[2]], 
          "IndicSp", 
          taxaLvl)

saveFiles(processedData, 
          sheet = relAbund(processedData), 
          plot = relAbundPlot(processedData), 
          "relAbund", 
          taxaLvl)








summary(indic2)

resultsANCOMBC2 = runANCOMBC2(processedData)


ANCOMBC2Results$res

ANCOMBC2Results$zero_ind








