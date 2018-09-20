#!/usr/bin/env Rscript

# load libs
source("libs.R")
# load processed data
source("processing.R")


## starting conc averages
tab %>% #
    filter(role=="Unknown" & hours==0) %>% #
    group_by(assay, season) %>% #
    summarise(mean=mean(copiesLitre))
    
## max cq
tab %>% #
    filter(role=="Unknown") %>% #
    select(cq) %>% #
    summary()
    
    
## summary table of env variables
treats %>% filter(!grepl('neg-', treatment)) %>% #
    group_by(season, treatment) %>% #
    mutate(bDNA=if_else(treatment=="e1-twothirds" | treatment=="harbour-twothirds", NaN, bDNA)) %>% # remove background DNA from mixes
    summarise(bDNA=round(unique(bDNA),digits=1), range_pH=paste(round(min(pH),digits=2), round(max(pH),digits=2), sep="--"), avg_ph=paste0(round(mean(pH),digits=2), " [", round(sd(pH),digits=2), "]"), #
    range_EC=paste(round(min(EC),digits=1), round(max(EC),digits=1),sep="--"), avg_EC=paste0(round(mean(EC),digits=1), " [", round(sd(EC),digits=2), "]")) %>% 
    select(-c(range_EC,range_pH)) %>% 
    xtable(caption="blahhh", digits=c(1,1,1,1,1,1)) %>%  print(include.rownames=FALSE)
# avgs
treats %>% filter(!grepl('neg-', treatment) & treatment!="asu") %>% #
    group_by(season) %>% 
    summarise(mean.ph=mean(pH), sd.ph=sd(pH), mean.ec=mean(EC), sd.ec=sd(EC))
#vs
treats %>% filter(!grepl('neg-', treatment) & treatment=="e1-full" | treatment=="harbour-full") %>% #
    group_by(treatment) %>% 
    summarise(mean.ph=mean(pH), sd.ph=sd(pH), mean.ec=mean(EC), sd.ec=sd(EC))

    
## summarise inhibition qPCR tests
treats %>% group_by(season, treatment) %>% filter(!is.na(inhib_effic)) %>% summarise(meanEff=mean(inhib_effic), maxEff=max(inhib_effic))
treats %>% filter(!is.na(inhib_effic)) %>% summarise(meanEff=mean(inhib_effic), maxEff=max(inhib_effic))


## calculate negative control amplications
# no NTC control amplifications 
tab %>% filter(role=="NTC") %>% select(cq) %>% summary()


## calculate bionegative control amplications
# total bionegs = 96
tab %>% filter(role=="BioNegative") %>% group_by(seasontank, assay, hours, treatment) %>% summarise(cnts=n()) %>% dim()

# negs with positive amplifications
tab %>% filter(role=="BioNegative" & !is.na(cq)) %>% #
    group_by(tank, assay, hours, treatment) %>% #
    summarise(cnts=n(), mpl=mean(copiesLitre)) %>% #
    arrange(treatment,hours,assay) %>% #
    #filter(assay=="crab-hex" & treatment=="neg-harbour") %>% #
    #arrange(seasontank, hours) %>% #
    print(n=Inf, width=Inf)

    
tab %>% filter(role=="BioNegative" & !is.na(cq) & treatment!="neg-harbour") %>% #
    group_by(assay) %>% #
    summarise(cnts=n(), mpl=mean(copiesLitre))

    
### Calulate probability of amplification - LOD/LOQ
# grab the assay efficiency from the qPCR output files with:
# { for i in *.csv; do echo "$(sed -n '9p' "$i")"; } done > output

# correct missing data in copies
tab <- tab %>% #
    mutate(copies=if_else(well=="A8" | well=="B8" | well=="C8", true=1, false=copies)) %>% #
    mutate(copies=if_else(well=="D7" | well=="E7" | well=="F7", true=10, false=copies)) %>% #
    mutate(copies=if_else(well=="A7" | well=="B7" | well=="C7", true=100, false=copies)) %>% #
    mutate(copies=if_else(well=="D6" | well=="E6" | well=="F6", true=1000, false=copies)) %>% #
    mutate(copies=if_else(well=="A6" | well=="B6" | well=="C6", true=10000, false=copies)) %>% #
    mutate(copies=if_else(well=="D5" | well=="E5" | well=="F5", true=100000, false=copies)) %>% #
    mutate(copies=if_else(well=="A5" | well=="B5" | well=="C5", true=1000000, false=copies))
    
# read in the assay efficiency
eff <- read_csv(file="../data/degradation-paper/efficiency.csv")
# summarise
eff %>% group_by(assay) %>% summarise(meanEff=mean(efficiency), sdEff=sd(efficiency), meanR2=mean(R2), sdR2=sd(R2))

# summarize the proportion of successful amplifications 
tab %>% filter(role=="Standard") %>% #
    group_by(assay, copies) %>% #
    summarise(totz=n(), countz=sum(!is.na(cq))) %>% #
    mutate(propz=countz/totz)


## zero amplifications 
pz <- function(x){round(length(which(x==0))/length(x), digits=2)}# proportion
pc <- function(x){length(which(!is.na(x)))}# number

tab %>% 
    filter(role=="Unknown") %>%
    group_by(treatment, hours, assay, season, tank) %>%
    summarise(numz = pc(copies)) %>% 
    #summarise(propz = pz(copies)) %>% 
    filter(numz < 3) %>%
    #arrange(hours, treatment, desc(propz), season, assay, tank) %>% 
    arrange(season, hours, tank, desc(numz), treatment, assay) %>% 
    #write_csv(path="data_missing.csv")    
    print(n=Inf)

# count the zero amps
tab %>% 
    filter(role=="Unknown") %>%
    group_by(season,hours) %>% 
    summarise(nas=sum(is.na(copiesLitre)), notNas=sum(!is.na(copiesLitre))) %>% 
    mutate(propNas=nas/notNas, tot=nas+notNas)
    