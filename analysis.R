library(gt23)
library(RMySQL)
library(dplyr)
library(ggplot2)
source('lib.R')

dbConn  <- dbConnect(MySQL(), group='specimen_management')
samples <- dbGetQuery(dbConn, 'select * from gtsp where Trial="Ghassemi_CART"')
#dbGetQuery(dbConn, 'select * from gtsp where SpecimenAccNum="GTSP0734"')

if(! file.exists('intSites.rds')){
  intSites <- getDBgenomicFragments(samples$SpecimenAccNum, 'specimen_management', 'intsites_miseq') %>%
            stdIntSiteFragments() %>%
            collapseReplicatesCalcAbunds() %>%
            annotateIntSites()
  saveRDS(intSites, file ='intSites.rds')
} else {
  intSites <- readRDS('intSites.rds')
}

# The heatmap makers sepearte samples by patient.
# Here we convert the GRange object to a data frame and use the patient column to 
# define which groups of samples we want to compare.
d <- data.frame(intSites)
d$patient <- ifelse(grepl('D7', d$timePoint), 'D7', 'D9')

# The epigenetic heatmap maker looks for a file named epiCellTypes which defines which epi-markers 
# to include. Available markers are located here, microb191:/data/internal/epigeneticData.
# Both the epigenetic and genomic heatmap makers look for an INSPIIRED control file to define
# the location of the epimakers and database settings which are now no longer needed for 
# the genomic heatmap maker. Software paths and config file paths are function paramters.

dir.create('epiGenHeatMap')
set.seed(2356)
createEpiGenomicHeatMapData(d, outputDir='epiGenHeatMap',Rscript_path = '/home/adrian/anaconda3/envs/r4/bin/Rscript',controlFile="INSPIIRED.yml",softwarePath = 'epi_heatmap_from_file.R')
epiGenomicHeatMap <- createGenomicHeatMapPlot('epiGenHeatMap', sampleOrder = c('D7', 'D9'))

dir.create('genHeatMap')
createGenomicHeatMapData(d, outputDir='genHeatMap',Rscript_path = '/home/adrian/anaconda3/envs/r4/bin/Rscript',softwarePath = 'genomic_heatmap_from_file.R',controlFile = 'INSPIIRED.yml')
genomicHeatMap <- createGenomicHeatMapPlot('genHeatMap', sampleOrder = c('D7', 'D9'))

# pdf(file='test.pdf')
# genomicHeatMap
# dev.off()
