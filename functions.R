#Script with functions


#Calling Packages

library(tidyverse)





#makes data ready for analysis
#Inputs(Raw data, taxaLVL, year)
processData = function(dataframe, taxaLvl, year){

  data = dataframe %>%
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
        Year = year) %>%
    #rename factors
    mutate(
        Environment = recode_factor(Environment, "e" = "Exposed", "u" = "Unexposed"),
         Timepoint = recode_factor(Timepoint, "t1" = "T1", "t2" = "T2", "t3" = "T3", "t4" = "T4")) %>%
    #bring everything to front
    select(sampleID, readType, Crop, Environment, Timepoint, Replicate, Year, everything())
  
return(data)
}


#Input: dataframe from processData
#Output: relative abundance
relAbund = function(dataframe){

metaData = data %>%
    select(sampleID, readType, Crop, Environment, Timepoint, Replicate, Year)
  
relData = data %>%
  select(-c(sampleID, readType, Crop, Environment, Timepoint, Replicate, Year)) %>%
  apply(1,function(x) x/sum(x)) %>% t() %>% as.data.frame() %>%
  merge(metaData) %>%
  select(sampleID, readType, Crop, Environment, Timepoint, Replicate, Year, everything())

return(relData)
}
    
  
preNMDS = function(dataframe){
  
}

plotNMDS = function{
  
}

permanovaHelper = function{
  
  
}



path = "~/AAFC WORK/act2_taxonomic_tables/act2_taxonomic_tables/app_2021_aggregated_counts.csv"
taxaLvl = "S"
dataframe = read.csv(path, header = TRUE)
year = str_extract(path, "2020|2021")

data = processData(dataframe, taxaLvl, year)
relData = relAbund(data)

