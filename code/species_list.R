#!/usr/bin/env Rscript
# script to download all UK decapod and fish sequences using a list from GBIF sequences from GenBank

# load libs
library("tidyverse")
library("rentrez")
library("rgbif")
library("rfishbase")
library("ggmap")
library("bold")
library("parallel")
library("ape")
library("taxize")
library("traits")
library("data.table")
#rm(list=ls())
sessionInfo()

## Gbif
# get classes for each group
classes <- c("Malacostraca","Actinopterygii","Elasmobranchii","Cephalaspidomorphi","Holocephali")

# get keys
keys <- sapply(classes, function(x) name_backbone(name=x)$classKey)

# set gbif account with your details
me <- ""
pass <- ""
mymail <- ""

# get a download ticket
# see https://github.com/ropensci/rgbif/blob/master/vignettes/downloads.Rmd
searches <- lapply(keys, function(x) occ_download(paste0("taxonKey=",x), "hasCoordinate=TRUE", "country=GB", user=me, pwd=pass, email=mymail))
# might have to run one by one and cat
crab.res <- occ_download(paste0("taxonKey=",keys[1]), "hasCoordinate=TRUE", "country=GB", user=me, pwd=pass, email=mymail)
fish.res <- occ_download(paste0("taxonKey=",keys[2]), "hasCoordinate=TRUE", "country=GB", user=me, pwd=pass, email=mymail)
shark.res <- occ_download(paste0("taxonKey=",keys[3]), "hasCoordinate=TRUE", "country=GB", user=me, pwd=pass, email=mymail)
lamp.res <- occ_download(paste0("taxonKey=",keys[4]), "hasCoordinate=TRUE", "country=GB", user=me, pwd=pass, email=mymail)
chim.res <- occ_download(paste0("taxonKey=",keys[5]), "hasCoordinate=TRUE", "country=GB", user=me, pwd=pass, email=mymail)
searches <- list(crab.res, fish.res, shark.res, lamp.res, chim.res)

# check the searches
meta.res <- lapply(searches, occ_download_meta)

# get citations
write.table(lapply(lapply(meta.res, gbif_citation), function(x) x$download), file="../temp/gbif_citations.txt", sep="\n")

# download
lapply(meta.res, function(x) occ_download_get(x$key, overwrite=TRUE))

# function reads in and cleans the data
zip_read_clean <- function(meta){
    unzip(zipfile=paste0(meta$key, ".zip"), files="occurrence.txt", overwrite=TRUE)
    tib <- as.tibble(fread("occurrence.txt", data.table=FALSE, sep="auto", select=c("basisOfRecord", "taxonRank", "order", "family", "genus","specificEpithet","decimalLatitude","decimalLongitude")))
    tib <- tib %>% filter(taxonRank=="SPECIES" & basisOfRecord!="FOSSIL_SPECIMEN") #%>% distinct(genus, specificEpithet)
    tib
}

# clean up
all.list <- lapply(meta.res, zip_read_clean)

# join
all.tib <- bind_rows(all.list)

# clean longlat
all.tib <- all.tib %>% filter(decimalLatitude > 40 & decimalLatitude < 70) %>% filter(decimalLongitude > -20 & decimalLongitude < 20)

# name new 
all.tib <- all.tib %>% mutate(species=paste(genus, specificEpithet))

# make unique and clean bad species names
all.clean <- unique(all.tib$species)[grep("[A-Z][a-z]* [a-z]*", unique(all.tib$species))]

# delete the files when not needed
#file.remove(sapply(meta.res, function(x) paste0(x$key,".zip")))

# get data off Genbank
range <- "300:1000"
query <- paste0("(", all.clean, "[ORGN] AND COI[ALL] AND ", range, "[SLEN]) OR (", all.clean, "[ORGN] AND cox1[ALL] AND ", range, "[SLEN]) OR (", all.clean, "[ORGN] AND CO1[ALL] AND ", range, "[SLEN])")

# search genbank
search.res <- mcmapply(FUN=function(x) entrez_search(db="nuccore", term=x, retmax=100), query, SIMPLIFY=FALSE, USE.NAMES=FALSE, mc.cores=8)

# removed the zero counts
search.res.nz <- search.res[which(lapply(search.res, function(x) x$count) > 0)]

# get IDs and remove dups
search.ids <- unique(unlist(lapply(search.res.nz, function(x) x$ids)))

# chunk up into 500s to stop server from rejecting request 
chunk <- 300
id.split <- unname(split(search.ids, ceiling(seq_along(search.ids)/chunk)))

# download using traits
ncbi.tabs <- mcmapply(FUN=ncbi_byid, id.split, SIMPLIFY=FALSE, USE.NAMES=FALSE, mc.cores=8)

# join all the data sets
ncbi.df <- as.tibble(bind_rows(ncbi.tabs))

# add the order and species
ncbi.df <- ncbi.df %>% mutate(order=all.tib$order[match(taxon, all.tib$species)], family=all.tib$family[match(taxon, all.tib$species)])

#write out a fasta
dtmp <- strsplit(ncbi.df$sequence, "")
names(dtmp) <- paste(ncbi.df$order, ncbi.df$family, str_replace_all(ncbi.df$taxon, " ", "_"), ncbi.df$acc_no, sep="_")
dat <- as.DNAbin(dtmp)
print(dat)
#write.dna(dat, file="../data/qpcr_specificity.fas", format="fasta", colw=99999)
#add the SEADNA seqs by hand

# check
dat <- read.dna(file="../data/qpcr_specificity.fas", format="fasta")
print(dat)#18,675 seqs
nam <- sapply(str_split(names(dat), pattern="_", n = Inf, simplify = FALSE), function(x) paste(x[3],x[4], sep="_"))
length(unique(nam))# 759 spp.
