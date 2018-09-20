#!/usr/bin/env Rscript
# clear memory
rm(list=ls())
# libs
library("tidyverse")
library("stringr")
library("forcats")
library("lubridate")
library("nlme")
library("MuMIn")
library("scales")
library("ggthemes")
library("xtable")
library("segmented")
library("gridExtra")
library("emmeans")
library("broom")
#library("broom.mixed")

#sink("../../temp/degradation-paper/RsessionInfo.txt")
#print(sessionInfo())
#sink()# stop sink

# disable printing sci notation
options(scipen = 999)