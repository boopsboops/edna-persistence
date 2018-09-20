#!/usr/bin/env Rscript

### load the data frames using readr

treats <- read_csv(file="../data/treatments_tanks.csv")#treatments per tank
tab <- read_csv(file="../data/qpcr_data.csv")# qpcr data
lit.rates <- read_csv(file="../data/literature_rates.csv")# rates of degrad from literature


## prep the data frame

# seperate the tank and hour labels and make new columns
tab <- tab %>% mutate(tank=str_split_fixed(string=sample, pattern="-", n=2)[,2], hours=str_split_fixed(string=sample, pattern="-", n=2)[,1])
#fix NAs
tab <- tab %>% mutate(tank=if_else(tank=="", NA_character_, tank), hours=if_else(hours=="", NA_character_, hours))

# create dummies for season and tank
tab <- tab %>% mutate(seasontank=if_else(role=="Unknown", true=paste0(season,tank), false=NA_character_))
treats <- treats %>% mutate(seasontank=paste0(season,random_tank_assignment))

# match and insert into data frame
tab <- tab %>% mutate(treatment=treats$treatment[match(seasontank, treats$seasontank)])

# add the environmental variables
tab <- tab %>% mutate(bDNA=treats$bDNA[match(seasontank, treats$seasontank)], #
    pH=treats$pH[match(seasontank, treats$seasontank)], #
    conductivity=treats$EC[match(seasontank, treats$seasontank)]) #%>% print(n=Inf, width=Inf)

# convert the number of hours into numeric hours
tab <- tab %>% mutate(hours=as.numeric(str_replace_all(hours ,"h", "")))

# add a biological negative control role
tab <- tab %>% mutate(role=if_else(grepl('neg-', treatment), true="BioNegative", false=role))

# add a new column for copies/L of seawater given 1 uL in each PCR reaction, 50 uL extraction volume, and 600 mL filtration volume
tab <- tab %>% mutate(copiesLitre = ((copies*50)*(1000/600)))

# make new seasontank
tab <- tab %>% mutate(tank=paste(season,tank,sep="."))

# check
#tab %>% head(n=100) %>% print(n=Inf, width=Inf)#


## count NAs 
tab %>% 
    filter(role=="Unknown") %>% #
    group_by(hours) %>% 
    summarise(nas=sum(is.na(cq)))

    
## prep data for model 
# get lowest copies/L
mcop <- min(filter(tab, role=="Unknown")$copiesLitre,na.rm=TRUE)# set to lowest conc
# replace missing data with zeros
# normalise
tab.means <- tab %>% 
    filter(role=="Unknown") %>% #
    filter(!(is.na(copiesLitre) & hours<96)) %>% #
    mutate(copiesLitre=if_else(is.na(copiesLitre), true=mcop, false=copiesLitre)) %>% # replac non-amplification NAs at 96 and 192 h to low concentrations (min copies in all reactions)#
    group_by(season, treatment, assay, hours, tank) %>% 
    summarise(meanCopiesLitre=mean(copiesLitre, na.rm=TRUE)) %>%
    filter(!(is.na(meanCopiesLitre))) %>% 
    ungroup() %>% 
    group_by(season, treatment, assay, tank) %>% 
    mutate(prop0=(meanCopiesLitre/meanCopiesLitre[hours==0]), startConc=meanCopiesLitre[hours==0]) %>% # normalise
    ungroup()

## minus ASU data
tab.means.asu <- tab.means %>% filter(treatment!="asu")

## minus 0h data
tab.means.rz <- tab.means %>% filter(hours!=0)

## minus 0h data and ASU
tab.means.asu.rz <- tab.means.rz %>% filter(treatment!="asu")

## add the env variables back in 
treats <- treats %>% mutate(season_tank=paste(season,random_tank_assignment,sep="."))
tab.means.asu.rz <- tab.means.asu.rz %>% mutate(bDNA=treats$bDNA[match(tank, treats$season_tank)], #
    pH=treats$pH[match(tank, treats$season_tank)], #
    conductivity=treats$EC[match(tank, treats$season_tank)]) #

#tab.means %>% head(n=100) %>% print(n=Inf, width=Inf)#glimpse()
