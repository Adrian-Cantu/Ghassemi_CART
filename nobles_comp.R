library(gt23)
library(RMySQL)
#library(dplyr)
#library(ggplot2)
library(tidyverse)
source('lib.R')

dbConn  <- dbConnect(MySQL(), group='specimen_management')


ss<-c("GTSP0560","GTSP0561","GTSP0562","GTSP0563","GTSP0564","GTSP0565","GTSP0566","GTSP0567","GTSP0568","GTSP0569","GTSP0570","GTSP0571","GTSP0572","GTSP0573","GTSP0574","GTSP0575",
      "GTSP0576","GTSP0577","GTSP0578","GTSP0579","GTSP0580","GTSP0581","GTSP0582","GTSP0583","GTSP0584","GTSP0585","GTSP0586","GTSP0587","GTSP0588","GTSP0589","GTSP0590","GTSP0591",
      "GTSP0592","GTSP0593","GTSP0594","GTSP0595","GTSP0596","GTSP0597","GTSP0598","GTSP0599","GTSP0600","GTSP0601","GTSP0602","GTSP0603","GTSP0604","GTSP0605","GTSP0606","GTSP0607",
      "GTSP0608","GTSP0609","GTSP0610","GTSP0611","GTSP0613","GTSP0614","GTSP0615","GTSP0616","GTSP0617","GTSP0619","GTSP0620","GTSP0621","GTSP0624","GTSP0625","GTSP0626","GTSP0627",
      "GTSP0628","GTSP0629","GTSP0630","GTSP0631","GTSP0632","GTSP0633","GTSP0634","GTSP0635","GTSP0636","GTSP0638","GTSP0639","GTSP0640","GTSP0641","GTSP0642","GTSP0643","GTSP0644",
      "GTSP0645","GTSP0646","GTSP0647","GTSP0648","GTSP0649","GTSP0650","GTSP0651","GTSP0652","GTSP0653","GTSP0654","GTSP0655","GTSP0734","GTSP0735","GTSP0736","GTSP0737","GTSP0738",
      "GTSP0739","GTSP0741","GTSP0742","GTSP0744","GTSP0746","GTSP0747","GTSP0748","GTSP0749","GTSP0750","GTSP0752","GTSP0754","GTSP0755","GTSP0757","GTSP0758","GTSP1166","GTSP1167",
      "GTSP1173","GTSP1175","GTSP1177","GTSP1178","GTSP1180","GTSP1183","GTSP1187","GTSP1188","GTSP1190","GTSP1191","GTSP1192","GTSP1193","GTSP1196","GTSP1197","GTSP1198","GTSP1203",
      "GTSP1204","GTSP1206","GTSP1207","GTSP1209","GTSP1210","GTSP1211","GTSP1212","GTSP1213","GTSP1215","GTSP1216","GTSP1219","GTSP1220","GTSP1222","GTSP1223","GTSP1224","GTSP1225",
      "GTSP1226","GTSP1228","GTSP1229","GTSP1232","GTSP1233","GTSP1234","GTSP1235","GTSP1236","GTSP1238","GTSP1240","GTSP1241","GTSP1405","GTSP1406","GTSP1407","GTSP1408","GTSP1409",
      "GTSP1410","GTSP1411","GTSP1413","GTSP1414","GTSP1415","GTSP1416","GTSP1603","GTSP1604","GTSP1605","GTSP1606","GTSP1607","GTSP1608","GTSP1609","GTSP2275","GTSP2276","GTSP2277",
      "GTSP2278","GTSP2279","GTSP2280","GTSP2281","GTSP2282","GTSP2283","GTSP2284","GTSP2285","GTSP2648","GTSP2649","GTSP2650","GTSP2651","GTSP2652","GTSP2653","GTSP2654","GTSP2655",
      "GTSP2656","GTSP2657","GTSP2658","GTSP2659","GTSP2660","GTSP2661","GTSP2662","GTSP2664","GTSP2665","GTSP2666")



arr <-paste0("SELECT * FROM gtsp WHERE SpecimenAccNum IN ('",paste(ss,collapse = "','"),"')")
nob_sample_tab <- dbGetQuery(dbConn,arr)
openxlsx::write.xlsx(nob_sample_tab, file = 'Nobles2020_samples.xlsx')
t0_samples <- nob_sample_tab %>% filter(Timepoint=='d0',SpecimenInfo=='Pre-Infusion_Product',VCN>0.17)

if(! file.exists('intSites_nob.rds')){
  intSites_nob <- getDBgenomicFragments(t0_samples$SpecimenAccNum, 'specimen_management', 'intsites_miseq') %>%
            stdIntSiteFragments() %>%
            collapseReplicatesCalcAbunds() %>%
            annotateIntSites()
  saveRDS(intSites_nob, file ='intSites_nob.rds')
} else {
  intSites_nob <- readRDS('intSites_nob.rds')
}
intSites <- readRDS('intSites.rds')
# The heatmap makers sepearte samples by patient.
# Here we convert the GRange object to a data frame and use the patient column to 
# define which groups of samples we want to compare.
d <- data.frame(intSites)
d$patient <- ifelse(grepl('D7', d$timePoint), 'D7', 'D9')

d_nob <- data.frame(intSites_nob)
dd <- rbind(d,d_nob)

# The epigenetic heatmap maker looks for a file named epiCellTypes which defines which epi-markers 
# to include. Available markers are located here, microb191:/data/internal/epigeneticData.
# Both the epigenetic and genomic heatmap makers look for an INSPIIRED control file to define
# the location of the epimakers and database settings which are now no longer needed for 
# the genomic heatmap maker. Software paths and config file paths are function paramters.

dir.create('epiGenHeatMap_dddd')
set.seed(2356)
createEpiGenomicHeatMapData(dd, outputDir='epiGenHeatMap_dddd',Rscript_path = '/home/adrian/anaconda3/envs/r4/bin/Rscript',controlFile="INSPIIRED.yml",softwarePath = 'epi_heatmap_from_file.R')
epiGenomicHeatMap <- createGenomicHeatMapPlot('epiGenHeatMap_dddd', sampleOrder = unique(dd$patient))

dir.create('genHeatMap_dddd')
createGenomicHeatMapData(dd, outputDir='genHeatMap_dddd',Rscript_path = '/home/adrian/anaconda3/envs/r4/bin/Rscript',softwarePath = 'genomic_heatmap_from_file.R',controlFile = 'INSPIIRED.yml')
genomicHeatMap <- createGenomicHeatMapPlot('genHeatMap_dddd', sampleOrder = unique(dd$patient))

# pdf(file='test.pdf')
# genomicHeatMap
# dev.off()

gene_dist <- 10000

onco_yes_nob <- (d_nob %>% filter(nearestOncoFeatureDist <  gene_dist) %>% summarise(t=sum(estAbund)))$t
onco_no_nob <- (d_nob %>% filter(nearestOncoFeatureDist >=  gene_dist) %>% summarise(t=sum(estAbund)))$t
onco_yes_nob + onco_no_nob

onco_yes_D7 <- (d %>% filter(patient=='D7' & nearestOncoFeatureDist <  gene_dist) %>% summarise(t=sum(estAbund)))$t
onco_no_D7  <- (d %>% filter(patient=='D7' & nearestOncoFeatureDist >= gene_dist) %>% summarise(t=sum(estAbund)))$t

onco_yes_D9 <- (d %>% filter(patient=='D9' & nearestOncoFeatureDist <  gene_dist) %>% summarise(t=sum(estAbund)))$t
onco_no_D9  <- (d %>% filter(patient=='D9' & nearestOncoFeatureDist >= gene_dist) %>% summarise(t=sum(estAbund)))$t


dat <- data.frame(
  "nob" = c(onco_yes_nob, onco_no_nob),
  "D7" = c(onco_yes_D7,onco_no_D7),
  'D9' = c(onco_yes_D9,onco_no_D9),
  row.names = c("onco_yes", "onco_no"),
  stringsAsFactors = FALSE
)

as.data.frame(t(dat)) %>% mutate(prop_onco=onco_yes/(onco_yes+onco_no))
chisq.test(dat) 

fisher.test(dat %>% select(nob,D7))$p.value

f_dat <- data.frame(
  "nob" = c(0, fisher.test(dat %>% select(nob,D7))$p.value,fisher.test(dat %>% select(nob,D9))$p.value),
  "D7" = c(fisher.test(dat %>% select(nob,D7))$p.value,0,fisher.test(dat %>% select(D9,D7))$p.value),
  'D9' = c(fisher.test(dat %>% select(nob,D9))$p.value,fisher.test(dat %>% select(D7,D9))$p.value,0),
  row.names = c("nob", "D7",'D9'),
  stringsAsFactors = FALSE
)
f_dat

## geting the order of the epifeatures from the svg file
## grep -A 3 -P 'xlink:href="row.\d+.svg"' main.svg | grep -v -e '^--$' | grep desc1 | perl -lne '/desc1>(.*)<\/desc1/; print $1' | tac | tr '\n' ',' | sed 's/,$/")\n/' | sed 's/,/","/' | sed 's/^/c("/'

epi_f_order <- c("Act_CD4_Tip60.10Kb","Rest_CD4_Tip60.10Kb","Rest_CD4_MOF.10Kb","Rest_CD4_p300.10Kb","Rest_CD4_PCAF.10Kb","Rest_CD4_CBP.10Kb","Act_CD4_HDAC6.10Kb","Rest_CD4_HDAC6.10Kb","Rest_CD4_HDAC3.10Kb","Rest_CD4_HDAC2.10Kb","Rest_CD4_HDAC1.10Kb","RestingNucleosomes.10Kb","ActivatedNucleosomes.10Kb","NRSF.10Kb","HP1b_promoters.10Kb","HP1a_promoters.10Kb","Brd4_promoters.10Kb","Brd3_promoters.10Kb","Brd2_promoters.10Kb",
                 "H4K91ac.10Kb","H4K8ac.10Kb","H4K5ac.10Kb","H4K16ac.10Kb","H4K12ac.10Kb","H3K9ac.10Kb","H3K4ac.10Kb","H3K36ac.10Kb","H3K27ac.10Kb","H3K23ac.10Kb","H3K18ac.10Kb","H3K14ac.10Kb","H2BK5ac.10Kb","H2BK20ac.10Kb","H2BK12ac.10Kb","H2BK120ac.10Kb","H2AK9ac.10Kb","H2AK5ac.10Kb","CTCF.10Kb","H2AZ.10Kb","PolII.10Kb","H4K20me3.10Kb","H4K20me1.10Kb","H4R3me2.10Kb","H3K79me3.10Kb","H3K79me2.10Kb","H3K79me1.10Kb","H3K36me3.10Kb","H3K36me1.10Kb","H3K27me3.10Kb","H3K27me2.10Kb","H3K27me1.10Kb","H3K9me3.10Kb","H3K9me2.10Kb","H3K9me1.10Kb","H3K4me3.10Kb","H3K4me2.10Kb","H3K4me1.10Kb","H3R2me2.10Kb","H3R2me1.10Kb","H2BK5me1.10Kb")
b_prot <- c("Act_CD4_Tip60.10Kb","Rest_CD4_Tip60.10Kb","Rest_CD4_MOF.10Kb","Rest_CD4_p300.10Kb","Rest_CD4_PCAF.10Kb","Rest_CD4_CBP.10Kb","Act_CD4_HDAC6.10Kb",
            "Rest_CD4_HDAC6.10Kb","Rest_CD4_HDAC3.10Kb","Rest_CD4_HDAC2.10Kb","Rest_CD4_HDAC1.10Kb","RestingNucleosomes.10Kb","ActivatedNucleosomes.10Kb","NRSF.10Kb",
            "HP1b_promoters.10Kb","HP1a_promoters.10Kb","Brd4_promoters.10Kb","Brd3_promoters.10Kb","Brd2_promoters.10Kb")
epi_g_order <-c('Bound Proteins','Histone Post-translational Modification')
kk <- readRDS('epiGenHeatMap_dddd/roc.res.rds')
epiROC <- as.data.frame(kk$ROC) %>% rownames_to_column(var = "feature") %>% pivot_longer(!feature,values_to='val', names_to='sample') %>%
  mutate(group=ifelse(feature %in% b_prot,'Bound Proteins','Histone Post-translational Modification')) %>%
  mutate(group=factor(group,levels=epi_g_order)) %>%
  mutate(feature=factor(feature,levels=rownames(kk$ROC))) %>% mutate(sample=factor(sample,levels=colnames(kk$ROC)))


png(height = 10, width = 10,units = 'in', res=300, file = 'epi_roc.png')
epiROC %>%
  ggplot( aes(sample,feature, fill= val)) + 
  geom_tile() +
  scale_fill_gradientn(colours=c('blue','white','red'),na.value = "transparent",breaks=c(0,0.5,1),labels=c(0,0.5,1), limits=c(0,1)) +
  facet_grid(rows=vars(epiROC$group),scales = "free_y", space = "free_y",switch='y') + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90),
        axis.title.x=element_blank(),
        axis.title.y=element_blank()) +
  labs(fill="ROC area") 
dev.off()

###

geno_kk <- readRDS('genHeatMap_dddd/roc.res.rds')
geno_f_order <- rownames(geno_kk$ROC)[3:nrow(geno_kk$ROC)]

geno_g_boundary <- geno_f_order[2:4]
geno_g_density <- geno_f_order[5:9]
geno_cg <- geno_f_order[10:14]
geno_cpg <- geno_f_order[15:19]
geno_dnase <- geno_f_order[20:23]
geno_g_order=c('cancer\nasso-\nciated\ngenes','gene boundary','gene density','CG content','CpG islands','DNaseI sites')

genoROC <- as.data.frame(geno_kk$ROC) %>% rownames_to_column(var = "feature") %>% filter(feature %in% geno_f_order ) %>%
  pivot_longer(!feature,values_to='val', names_to='sample') %>%
  mutate(group=ifelse(feature %in% geno_g_boundary,geno_g_order[2],
               ifelse(feature %in% geno_g_density,geno_g_order[3],
               ifelse(feature %in% geno_cg,geno_g_order[4],
               ifelse(feature %in% geno_cpg,geno_g_order[5],
               ifelse(feature %in% geno_dnase,geno_g_order[6],
                      geno_g_order[1])))))) %>%
  mutate(group=factor(group,levels=geno_g_order)) %>%
  mutate(feature=factor(feature,levels=rev(geno_f_order),ordered=TRUE)) %>% mutate(sample=factor(sample,levels=colnames(geno_kk$ROC)))

png(height = 10, width = 10,units = 'in', res=300, file = 'geno_roc.png')
genoROC %>%
  ggplot( aes(sample,feature, fill= val)) + 
  geom_tile() +
  scale_fill_gradientn(colours=c('blue','white','red'),na.value = "transparent",breaks=c(0,0.5,1),labels=c(0,0.5,1), limits=c(0,1)) +
  facet_grid(rows=vars(genoROC$group),scales = "free_y", space = "free_y",switch='y') + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90),
        axis.title.x=element_blank(),
        axis.title.y=element_blank()) +
  labs(fill="ROC area") 
dev.off()

