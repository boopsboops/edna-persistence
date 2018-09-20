#!/usr/bin/env Rscript

## load libs
source("libs.R")# this also removes all R objects in memory
## load processed data
source("processing.R")# this loads up the data


## choose a dataset 
#tab.tmp <- tab # full processed data
#tab.tmp <- tab.means # full processed data with per tank averages
#tab.tmp <- tab.means.asu # with ASU data removed
#tab.tmp <- tab.means.rz# with zero hours data removed
tab.tmp <- tab.means.asu.rz # with zero hours AND asu data removed (with the env variables added back)


## Start with a simple generalised least squares model, including all possible main and interactive effects of your experimental factors
m0 <- gls(log(prop0) ~ hours * assay * treatment * season * log(startConc), data=tab.tmp, method="REML", control=lmeControl(opt="optim",msMaxIter=500)) # model with no random effects
#summary(m0)

## model diagnostics
    par(mfrow=c(2,3),mar=c(4,4,3,1))
    qqnorm(resid(m0)); qqline(resid(m0))
    plot(resid(m0) ~ fitted(m0), main="Fitted v Residuals", xlab="Fitted values", ylab="Residuals"); abline(0,0)
    plot(resid(m0) ~ factor(tab.tmp$hours), main="Hours v Residuals", xlab="Hours", ylab="Residuals"); abline(0,0)
    plot(resid(m0) ~ factor(tab.tmp$assay), main="Assay v Residuals", xlab="Assay", ylab="Residuals"); abline(0,0)
    plot(resid(m0) ~ factor(tab.tmp$treatment), main="Treatment v Residuals", xlab="Treatment", ylab="Residuals", xaxt="n"); abline(0,0)
    axis(1, at=1:4, labels=c("e1_f","e1_2/3","h_f","h_2/3"))
    plot(resid(m0) ~ factor(tab.tmp$season), main="Season v Residuals", xlab="Season", ylab="Residuals"); abline(0,0)


## Find the optimal variance weighting for your model (do this with the full fixed effects structure and method='REML'). Use varIdent for categorical variables.
m0.w0 <- gls(log(prop0) ~ hours * assay * treatment * season * log(startConc), data=tab.tmp, method="REML", control=lmeControl(opt="optim",msMaxIter=500))   # model with no random effects
m0.w1 <- gls(log(prop0) ~ hours * assay * treatment * season * log(startConc), weights=varIdent(form=~1|assay), data=tab.tmp, method="REML", control=lmeControl(opt="optim",msMaxIter=500))   # model with no random effects
m0.w2 <- gls(log(prop0) ~ hours * assay * treatment * season * log(startConc), weights=varIdent(form=~1|treatment), data=tab.tmp, method="REML", control=lmeControl(opt="optim",msMaxIter=500))   # model with no random effects
m0.w3 <- gls(log(prop0) ~ hours * assay * treatment * season * log(startConc), weights=varIdent(form=~1|season), data=tab.tmp, method="REML", control=lmeControl(opt="optim",msMaxIter=500))   # model with no random effects
m0.w4 <- gls(log(prop0) ~ hours * assay * treatment * season * log(startConc), weights=varIdent(form=~1|assay*treatment), data=tab.tmp, method="REML", control=lmeControl(opt="optim",msMaxIter=500))   # model with no random effects
m0.w5 <- gls(log(prop0) ~ hours * assay * treatment * season * log(startConc), weights=varIdent(form=~1|assay*season), data=tab.tmp, method="REML", control=lmeControl(opt="optim",msMaxIter=500))   # model with no random effects
m0.w6 <- gls(log(prop0) ~ hours * assay * treatment * season * log(startConc), weights=varIdent(form=~1|treatment*season), data=tab.tmp, method="REML", control=lmeControl(opt="optim",msMaxIter=500))   # model with no random effects
m0.w7 <- gls(log(prop0) ~ hours * assay * treatment * season * log(startConc), weights=varIdent(form=~1|treatment*season*assay), data=tab.tmp, method="REML", control=lmeControl(opt="optim",msMaxIter=500))   # model with no random effects
model.sel(m0.w0,m0.w1,m0.w2,m0.w3,m0.w4,m0.w5,m0.w6,m0.w7) %>% as_tibble(rownames="m") %>%  select(m, delta)

# Allowing different variance weighting for each treatment*season level significantly improves model fitting
m0.w8 <- gls(log(prop0) ~ hours * assay * treatment * season * log(startConc), weights=varIdent(form=~1|treatment*season), cor=corAR1(), data=tab.tmp, method="REML", control=lmeControl(opt="optim",msMaxIter=500))   # model with no random effects
model.sel(m0.w6,m0.w8) %>% as_tibble(rownames="m") %>%  select(m, delta)
# time autocorelation does not

## Find the optimal random effects structure (do this with the full fixed effects structure and method='REML')
m.r1 <- lme(log(prop0) ~ hours * assay * treatment * season * log(startConc), random=~1|tank, weights=varIdent(form=~1|treatment*season), data=tab.tmp, method="REML", control=lmeControl(opt="optim",msMaxIter=500))   # same degradation rate for each tank, but with different starting values
m.r2 <- lme(log(prop0) ~ hours * assay * treatment * season * log(startConc), random=~1+hours|tank, weights=varIdent(form=~1|treatment*season), data=tab.tmp, method="REML", control=lmeControl(opt="optim",msMaxIter=500))   # different degradation rates for each tank
model.sel(m0.w6,m.r1,m.r2) %>% as_tibble(rownames="m") %>%  select(m, delta)
# Random intercept and slope model is best


## new model selection
# create all conbinations
x1 <- paste(combn(c("hours","assay","season","treatment","log(startConc)"), 1, FUN=paste, collapse=':'), collapse=" + ")
x2 <- paste(combn(c("hours","assay","season","treatment","log(startConc)"), 2, FUN=paste, collapse=':'), collapse=" + ")
x3 <- paste(combn(c("hours","assay","season","treatment","log(startConc)"), 3, FUN=paste, collapse=':'), collapse=" + ")
x4 <- paste(combn(c("hours","assay","season","treatment","log(startConc)"), 4, FUN=paste, collapse=':'), collapse=" + ")
x5 <- paste(combn(c("hours","assay","season","treatment","log(startConc)"), 5, FUN=paste, collapse=':'), collapse=" + ")
f0 <- "log(prop0) ~ "
form <- as.formula(paste0(f0,paste(x1,x2,x3,x4,x5,sep=" + ")))
# fit model
mm <- lme(fixed=form, random=~1+hours|tank, weights=varIdent(form=~1|treatment*season), data=tab.tmp, method="ML", control=lmeControl(opt="optim",msMaxIter=500))
mm %>% anova

# stepwise deletion of interactions
mm <- update(mm, . ~ . - hours:assay:season:treatment:log(startConc)); anova(mm)
mm <- update(mm, . ~ . - assay:season:treatment:log(startConc)); anova(mm)
mm <- update(mm, . ~ . - hours:season:treatment:log(startConc)); anova(mm)
mm <- update(mm, . ~ . - hours:assay:season:log(startConc)); anova(mm)
mm <- update(mm, . ~ . - hours:assay:season:treatment); anova(mm)
mm <- update(mm, . ~ . - hours:assay:treatment:log(startConc)); anova(mm)
mm <- update(mm, . ~ . - season:treatment:log(startConc)); anova(mm)
mm <- update(mm, . ~ . - assay:treatment:log(startConc)); anova(mm)
mm <- update(mm, . ~ . - assay:season:log(startConc)); anova(mm)
mm <- update(mm, . ~ . - assay:season:treatment); anova(mm)
mm <- update(mm, . ~ . - hours:treatment:log(startConc)); anova(mm)
mm <- update(mm, . ~ . - hours:season:log(startConc)); anova(mm)
mm <- update(mm, . ~ . - hours:season:treatment);anova(mm)
mm <- update(mm, . ~ . - hours:assay:log(startConc)); anova(mm)
mm <- update(mm, . ~ . - hours:assay:treatment); anova(mm)
mm <- update(mm, . ~ . - hours:assay:season); anova(mm)
mm <- update(mm, . ~ . - treatment:log(startConc)); anova(mm)
mm <- update(mm, . ~ . - season:log(startConc)); anova(mm)
mm <- update(mm, . ~ . - assay:log(startConc));anova(mm)
#mm <- update(mm, . ~ . - hours:log(startConc));anova(mm)
# all interactions now significant!

# to print the final model 
mm$call

## Switch to method='REML' for final model output 
# run new model 
m0 <- lme(log(prop0) ~ hours + assay + season + treatment + 
    log(startConc) + hours:assay + hours:season + hours:treatment + 
    hours:log(startConc) + assay:season + assay:treatment + season:treatment,
    random=~1+hours|tank, weights=varIdent(form=~1|treatment*season), 
    data=tab.tmp, method="REML", control=lmeControl(opt="optim",msMaxIter=500))


## model diagnostics 
#pdf("../temp/model_diagnostics.pdf", width=12, height=9, useDingbats=FALSE)
    par(mfrow=c(2,3),mar=c(4,4,3,1))
    qqnorm(resid(m0)); qqline(resid(m0))
    plot(resid(m0, type="normalized") ~ fitted(m0), main="Fitted v Residuals", xlab="Fitted values", ylab="Residuals"); abline(0,0)
    plot(resid(m0, type="normalized") ~ factor(tab.tmp$hours), main="Hours v Residuals", xlab="Hours", ylab="Residuals"); abline(0,0)
    plot(resid(m0, type="normalized") ~ factor(tab.tmp$assay), main="Assay v Residuals", xlab="Assay", ylab="Residuals"); abline(0,0)
    plot(resid(m0, type="normalized") ~ factor(tab.tmp$treatment), main="Treatment v Residuals", xlab="Treatment", ylab="Residuals", xaxt="n"); abline(0,0)
    axis(1, at=1:4, labels=c("e1_f","e1_2/3","h_f","h_2/3"))
    plot(resid(m0) ~ factor(tab.tmp$season), main="Season v Residuals", xlab="Season", ylab="Residuals"); abline(0,0)
#dev.off()

# Looks good, apart from potential heavy-tailed distribution of residuals? Save model output to file
#sink(file="../temp/model_diagnostics.txt")
#summary(m0);anova(m0)
#sink()


## posthoc tests

emt.hours <- emtrends(model=m0, specs="hours", var="hours")
emt.assay <- emtrends(model=m0, specs="assay", var="hours")
emt.season <- emtrends(model=m0, specs="season", var="hours")
emt.treatment <- emtrends(model=m0, specs="treatment", var="hours")
emt.assay.season <- emtrends(model=m0, specs=c("assay","season"), var="hours")
emt.sat <- emtrends(model=m0, specs=c("season","assay","treatment"), var="hours")


# FUN to calc half lives
half_life <- function(slope){
    hl <- log(0.5)/slope
    return(hl)
    }#

# FUN to add the half lives to the table
add_half_life <- function(df){
    df <- df %>% as_tibble %>% mutate(est.halflife=half_life(hours.trend), u.conf.halflife=half_life(upper.CL), l.conf.halflife=half_life(lower.CL))
    return(df)
    }#

# list them
emt.list <- list(emt.hours, emt.assay, emt.season, emt.treatment, emt.assay.season, emt.sat)

# add and bind them
emt.all <- bind_rows(lapply(emt.list, add_half_life))

print_clean <- function(df){
    df <- df %>% mutate(coefsCI=paste0(round(hours.trend, digits=3), " [",round(upper.CL, digits=3),",",round(lower.CL, digits=3),"]"), #
        halfsCI=paste0(round(est.halflife, digits=1), " [",round(l.conf.halflife, digits=1),",",round(u.conf.halflife, digits=1),"]"))
        df <- df %>% mutate(treatment=if_else(is.na(treatment), "All", as.character(treatment)), #
            season=if_else(is.na(season), "All", as.character(season)), #
            assay=if_else(is.na(assay), "All", as.character(assay)))
    df <- df %>% select(treatment, season, assay, coefsCI, halfsCI)
    df <- df %>% arrange(treatment,season,assay)
    return(df)
    }#

# print for xtable
emt.all %>% print_clean() %>% xtable(caption="blah") %>% print(include.rownames=FALSE)


# plot
pl4 <- emt.all %>% 
    filter(!is.na(treatment) & !is.na(assay)) %>% 
    mutate(assayseason=paste0(season,assay)) %>% 
    ggplot(aes(x=parse_factor(assayseason,levels=c("summercrab-hex","wintercrab-hex","summershanny-fam","wintershanny-fam")), y=est.halflife, colour=parse_factor(treatment,levels=c("e1-full","e1-twothirds","harbour-twothirds","harbour-full")))) + 
    geom_pointrange(aes(ymin=l.conf.halflife, ymax=u.conf.halflife), position=position_dodge(0.55), size=1.3) + 
    #geom_point(aes(x=assayseason, y=est),data=lme.rates, position=position_dodge(0.55), shape=3, size=7) +
    theme_bw() + scale_color_manual(values=c(ptol_pal()(5)[c(2,3,4,5)])) + #
    scale_y_continuous(limits=c(0,70), breaks=seq(0,70,10)) +
        theme(axis.text.x=element_text(size=12), #
        axis.text.y=element_text(size=12), #
        axis.title.x=element_text(size=14), #
        axis.title.y=element_text(size=14), #
        strip.text.x=element_text(size = 14), #
        strip.text.y=element_text(size = 14), #
        strip.background=element_blank(), #
        legend.text=element_text(size=12), #
        legend.title=element_blank()) +
        labs(x="Species assay by season", y="eDNA half life (hours)")
plot(pl4)
#ggsave(filename="../temp/half-lives.svg", height=6, width=10, plot=pl4)


## statistical tests on slopes

# make pairwise comparisons, and bind
# remove all
emt.pairs <- bind_rows(lapply(lapply(emt.list[c(2,3,4)], pairs), as_tibble))

# format
emt.pairs <- emt.pairs %>% separate(col=contrast, into=c("contrast1","contrast2"), sep = " - ") 

# format for xtable
emt.pairs %>% 
    mutate(p=as.character(round(p.value, digits=4)), 
    estimate=paste0(as.character(round(estimate, digits=3)), " [", as.character(round(SE, digits=3)), "]"),
    t.ratio=as.character(round(t.ratio, digits=3))) %>% 
    select(-c(df,p.value,SE)) %>% 
    add_column(Variable="x", .before="contrast1") %>% 
    xtable(caption="blahhh", digits=c(3,3,3,3,3,3,3)) %>%  print(include.rownames=FALSE)


## Plot

# add in the fitted/fixed values
m0.aug <- m0 %>% augment(data=tab.tmp, effects="fixed") %>% as.tibble()
# convert to factors
m0.aug <- m0.aug %>% mutate(treatment=parse_factor(treatment, levels=c("e1-full","e1-twothirds","harbour-twothirds","harbour-full")), 
    assay=parse_factor(assay,levels=c("shanny-fam","crab-hex")), 
    season=parse_factor(season,levels=c("winter","summer")))
m0.aug <- m0.aug %>% group_by(season,treatment,assay,hours) %>% mutate(.fitted.avg=mean(.fitted), .fixed.avg=mean(.fixed))
m0.aug %>% print(n=40)


# Plot Assay*Season (see treatments)
pl1 <- m0.aug %>% ggplot(aes(x=hours, y=log(prop0))) +
    geom_point(aes(colour=treatment), position=position_dodge(12), size=3, alpha=0.9, pch=16) +
    geom_line(aes(x=hours,y=.fitted.avg,colour=treatment), size=1,alpha=1) +
    #geom_smooth(aes(colour=treatment), method=lm, se=FALSE, size=1, alpha=0.09) + #linetype="dashed", 
    facet_grid(rows=vars(assay), cols=vars(season)) +
    scale_x_continuous(breaks=c(0,12,24,48,96,192),labels=c(0,12,24,48,96,192)) + #
    theme_bw() + scale_color_manual(values=c(ptol_pal()(5)[c(2,3,4,5)])) + # just the four of the ptol colours
    theme(axis.text.x=element_text(size=12), #
        axis.text.y=element_text(size=12), #
        axis.title.x=element_text(size=14), #
        axis.title.y=element_text(size=14), #
        strip.text.x=element_text(size = 12), #
        strip.text.y=element_text(size = 12), #
        strip.background=element_blank(), #
        #panel.grid.major = element_blank(),
        #panel.grid.minor = element_blank(),
        legend.text=element_text(size=12), #
        legend.title=element_blank(),
        plot.background=element_blank()) +
        labs(x="Hours", y="Proportion of time-zero eDNA concentration\n (natural log transformed)")
plot(pl1)
#ggsave(filename="../temp/rates.svg", width=12, height=9, plot=pl1)


# plot on a NON-LOG SCALE
tab.tmp <- tab.means
pl3 <- tab.tmp %>% ggplot(aes(x=hours, y=meanCopiesLitre)) +#
    geom_point(aes(colour=parse_factor(treatment, levels=c("asu","e1-full","e1-twothirds","harbour-twothirds","harbour-full"))), position=position_dodge(12), size=3, alpha=0.9, pch=16) + #
    geom_smooth(aes(colour=parse_factor(treatment, levels=c("asu","e1-full","e1-twothirds","harbour-twothirds","harbour-full"))), method="glm", method.args=list(family=gaussian(link="log")), se=FALSE) +
    scale_x_continuous(breaks=c(0,12,24,48,96,192),labels=c(0,12,24,48,96,192)) + #
    scale_y_continuous(labels=comma) +
    facet_grid(rows=vars(parse_factor(assay,levels=c("shanny-fam","crab-hex"))), cols=vars(parse_factor(season, levels=c("winter","summer"))), scales="free") +
    theme_bw() + scale_color_ptol() + scale_fill_ptol() + #
    theme(axis.text.x=element_text(size=12), #
        axis.text.y=element_text(size=12), #
        axis.title.x=element_text(size=14), #
        axis.title.y=element_text(size=14), #
        strip.text.x=element_text(size = 14), #, face="bold"
        strip.text.y=element_text(size = 14), #, face="bold"
        strip.background=element_blank(), #
        legend.text=element_text(size=12), #
        legend.title=element_blank(), #
        plot.background=element_blank()) +
        labs(x="Hours", y="eDNA copies/L")
plot(pl3)
#ggsave(filename="../temp/non-log.svg", width=12, height=9, plot=pl3)


## test sig of other variables and ASU synthetic
# reset data
tab.tmp <- tab.means.asu.rz
# simple model for covariates # weights=varIdent(form=~1|treatment*season)
m1 <- lme(log(prop0) ~ hours + assay + bDNA + pH + conductivity + log(startConc), random=~1+hours|tank, data=tab.tmp, method="REML", control=lmeControl(opt="optim",msMaxIter=500))
m1 <- lme(log(prop0) ~ hours + season + treatment + assay + bDNA + pH + conductivity + log(startConc), random=~1+hours|tank, data=tab.tmp, method="REML", control=lmeControl(opt="optim",msMaxIter=500))
anova(m1); summary(m1); fixef(m1)


# make a quick model of the synthetic degrad
tab.tmp <- tab.means.rz %>% filter(treatment=="asu")
m2 <- lme(fixed=log(prop0) ~ hours + assay + season + log(startConc), random=~1+hours|tank, weights=varIdent(form=~1|season), data=tab.tmp, method="REML", control=lmeControl(opt="optim",msMaxIter=500))
emt.hours <- emtrends(model=m2, specs="hours", var="hours")
emt.assay.season <- emtrends(model=m2, specs=c("assay","season"), var="hours")
half_life(as_tibble(emt.hours)$hours.trend); half_life(as_tibble(emt.hours)$lower.CL); half_life(as_tibble(emt.hours)$upper.CL)
summary(m2); fixef(m2); anova(m2)

# plot asu
m2.aug <- m2 %>% augment(data=tab.tmp, effects="fixed") %>% as.tibble()
m2.aug <- m2.aug %>% mutate(assay=parse_factor(assay,levels=c("shanny-fam","crab-hex")), 
    season=parse_factor(season,levels=c("winter","summer")))
m2.aug <- m2.aug %>% group_by(season,assay,hours) %>% mutate(.fitted.avg=mean(.fitted), .fixed.avg=mean(.fixed))

pl.m2 <- m2.aug %>% ggplot(aes(x=hours, y=log(prop0))) +
    geom_point(aes(colour=season), position=position_dodge(12), size=3, alpha=0.9, pch=16) +
    geom_line(aes(x=hours,y=.fitted.avg,colour=season, group=tank), size=1,alpha=1) +
    #geom_smooth(aes(colour=season), method=lm, se=FALSE, size=1, alpha=0.09, linetype="dashed") + #linetype="dashed", 
    facet_grid(rows=vars(assay)) +
    scale_x_continuous(breaks=c(0,12,24,48,96,192),labels=c(0,12,24,48,96,192)) + #
    theme_bw() + scale_color_ptol() + # just the four of the ptol colours
    theme(axis.text.x=element_text(size=12), #
        axis.text.y=element_text(size=12), #
        axis.title.x=element_text(size=14), #
        axis.title.y=element_text(size=14), #
        strip.text.x=element_text(size = 12), #
        strip.text.y=element_text(size = 12), #
        strip.background=element_blank(), #
        #panel.grid.major = element_blank(),
        #panel.grid.minor = element_blank(),
        legend.text=element_text(size=12), #
        legend.title=element_blank(),
        plot.background=element_blank()) +
        labs(x="Hours", y="Proportion of time-zero eDNA concentration\n (natural log transformed)")
plot(pl.m2)
#ggsave(filename="../temp/synthetic.svg", width=7, height=9, plot=pl.m2)


# quick plot of startConcs
# crab lower startconc than shanny, winter higher startconc than summer
# crab degrades faster than shanny, winter degrades faster than summer
tab.tmp %>% mutate(seasonAssay=paste0(season,assay)) %>% 
ggplot(aes(x=tank, y=log(startConc))) +
    geom_point(aes(colour=seasonAssay), size=5)
    