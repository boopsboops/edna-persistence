#!/usr/bin/env Rscript

## load libs
source("libs.R")
## load processed data
source("processing.R")

## Literature review of half lives
#half life converter 
convert_to_halflife <- function(rate,period){
    if(period=="days"){
    half <- round(log(0.5)/(-abs(rate/24)), digits=1)
    }
    if(period=="hours"){
    half <- round(log(0.5)/(-abs(rate)), digits=1)
    }
    return(half)
}

# latex table
lit.rates %>% 
    mutate(halfLifeLit=if_else(period=="hours", true=convert_to_halflife(rate, period="hours"), false=convert_to_halflife(rate, period="days"))) %>% 
    mutate(halfLife=if_else(is.na(halfLife), true=halfLifeLit, false=halfLife)) %>% 
    mutate(study=paste0(str_replace_all(string=study, pattern=" \\([0-9][0-9][0-9][0-9]\\)", replacement=""), str_replace_all(string=studyTex, pattern="citet\\{", replacement="cite{"))) %>% 
    select(-c(studyTex,halfLifeLit,rate,period)) %>% 
    arrange(halfLife) %>% 
    xtable(caption="blah", digits=c(0,0,0,0,0,0,0,1,1)) %>%  print.xtable(include.rownames=FALSE, sanitize.text.function=identity, digits=1)

# grob table with plot
lit.rates.tab <- lit.rates %>% 
    mutate(halfLifeLit=if_else(period=="hours", true=convert_to_halflife(rate, period="hours"), false=convert_to_halflife(rate, period="days"))) %>% 
    mutate(halfLife=if_else(is.na(halfLife), true=halfLifeLit, false=halfLife)) %>% 
    mutate(halfLife=round(halfLife, digits=1)) %>% 
    select(-c(studyTex,halfLifeLit,rate,period)) %>% 
    arrange(halfLife) %>% 
    mutate(key=parse_factor(paste(study,organism,waterSource,environment,fragment,temperature,pH,halfLife), levels=paste(study,organism,waterSource,environment,fragment,temperature,pH,halfLife)))


# set cols -- http://tools.medialab.sciences-po.fr/iwanthue/index.php
mcols <- c("#4fa500","#00a0ee")

# plot the lit rates fig
pl5 <- lit.rates.tab %>% ggplot(aes(x=key, y=halfLife)) +
    geom_linerange(aes(ymin=0, ymax=halfLife, colour=environment), size=2) + 
    scale_y_continuous(position="left", trans="log10") + #
        theme(axis.title.y=element_blank(), 
        #axis.ticks.y=element_blank(),
        axis.ticks.x=element_blank(),
        #axis.text.y=element_blank(),
        axis.text.x=element_blank(),
        axis.title.x=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
    scale_color_manual(values=mcols) 
plot(pl5)
#ggsave(filename="../temp/lit_plot.svg", plot=pl5, width=297, height=80, units="mm")

# make the table
pl6 <- lit.rates.tab %>% select(-c("environment","key")) %>% 
    mutate(temperature=replace_na(temperature, ""), pH=replace_na(pH, "")) %>% 
    mutate(temperature=str_replace_all(temperature, "--", "-"), pH=str_replace_all(pH, "--", "-"), fragment=str_replace_all(fragment, "--", "-")) %>%
    rename(Study=study, "Water\nsource"=waterSource, Organism=organism, "Fragment\nlength (bp)"=fragment, "Temperature\n(degC)"=temperature, "Half life\n(hours)"=halfLife) %>% 
    tableGrob(rows=NULL, theme=ttheme_default(base_size=8, core=list(fg_params=list(hjust=0, x=0.07)), colhead=list(fg_params=list(hjust=0, x=0.07))))
plot(pl6)
#ggsave(filename="../temp/figureslit_tab.svg", plot=pl6, height=500, width=200, units="mm")
