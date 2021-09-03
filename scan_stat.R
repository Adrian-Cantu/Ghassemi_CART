library(gt23)
library(RMySQL)
library(dplyr)
library(ggplot2)
source('lib.R')

library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(Homo.sapiens)
all_genes <- genes(TxDb.Hsapiens.UCSC.hg38.knownGene)


list_genes <- function(chr,g_start,g_end) {
  temp <- intSites %>% GenomicRanges::as.data.frame(.) %>% filter(seqnames==chr & start>=g_start & start<=g_end)
  return(paste(unique(temp$nearestFeature),collapse = ' '))
}

list_all_genes <- function(seqnames,start,end,genes) {
  tmp_R <-GRanges(paste0(seqnames,':',start,'-',end,':*'))
  tmp <- subsetByOverlaps(genes, tmp_R)
  tmp_g <-  as.data.frame(org.Hs.egSYMBOL) %>% filter(gene_id %in% tmp$gene_id)
  #return(list(tmp$gene_id))
  #return(paste(tmp$gene_id,collapse = ' '))
  return(paste(tmp_g$symbol,collapse = " "))
}



dbConn  <- dbConnect(MySQL(), group='specimen_management')
samples <- dbGetQuery(dbConn, 'select * from gtsp where Trial="Ghassemi_CART"')

if(! file.exists('intSites.rds')){
  intSites <- getDBgenomicFragments(samples$SpecimenAccNum, 'specimen_management', 'intsites_miseq') %>%
    stdIntSiteFragments() %>%
    collapseReplicatesCalcAbunds() %>%
    annotateIntSites()
  saveRDS(intSites, file ='intSites.rds')
} else {
  intSites <- readRDS('intSites.rds')
}

d <- data.frame(intSites)
d$patient <- ifelse(grepl('D7', d$timePoint), 'D7', 'D9')

g1 <- d %>% filter(patient == 'D7' )  %>% makeGRangesFromDataFrame(.,keep.extra.columns=TRUE)
g2 <- d %>% filter(patient == 'D9' )  %>% makeGRangesFromDataFrame(.,keep.extra.columns=TRUE)

g1 %>% GenomicRanges::as.data.frame() %>%nrow()
g2 %>% GenomicRanges::as.data.frame() %>%nrow()


d7vsd9 <- scanStats(g1,g2,gr1.label = "D7",gr2.label = "D9", kvals = "4L:20L") %>%
  GenomicRanges::as.data.frame(.) %>%
  arrange(clusterSource) %>%
#  filter(width < 1000000)  %>%
  rowwise() %>%
  mutate(genesIntsites=list_genes(as.character(seqnames),start,end),
         genesEntrez=list_all_genes(seqnames = seqnames,start = start,end = end,all_genes)) %>%
  ungroup()
d7vsd9
