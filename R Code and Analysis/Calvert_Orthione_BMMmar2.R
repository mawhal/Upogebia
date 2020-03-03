###################################################################
# Hakai Institute - Calvert Island Observatory
# Biodiversity and range shifts
# Orthione griffensis infestataion of Upogebia pugettensis
# documenting the presence and severity of infestation
# first observation in July 2018
# additional data collected in August 2018 to document prevalence
# All collections made in Pruth Bay

# load libraries
library(tidyverse)
library(ggplot2)
library(cowplot)
library(lubridate)
library(lme4)
library(plyr)
library(psych) # for logistic function
## read data sets --------------------------------------------------------------------------------------------------------
#@ first calvert sampling 2018
calvert1 <- read.csv( "C:/Users/FABS/Documents/Ubpgebia/data/PugettensisCalvertAugV1D1.csv", stringsAsFactors = FALSE )
# remove empty rows
calvert1 <- calvert1[ calvert1$Date != "", ]
# select and rename columns
d1 <- calvert1 %>% select( site=Site, date=Date, len.total=Size.mm., parasite.sex=ParasiteType )
d1$location = "Calvert Island"
d1$site.type = "Reference Site"
d1$method <- "excavation"
d1$date <- dmy( d1$date )
d1$len.cara <- 0.3 * d1$len.total # this conversion comes from John Chapman

# convert parasite sex column to two presence-absence columns
fmmat <- matrix(NA,ncol=2,nrow=nrow(d1))
for( i in 1:nrow(d1) ){
  if( d1$parasite.sex[i]=="Both" )  fmmat[i,] <- c(1,1)
  if( d1$parasite.sex[i]=="Female" )  fmmat[i,] <- c(1,0)
  if( d1$parasite.sex[i]=="None" )  fmmat[i,] <- c(0,0)
}
d1$pf <- fmmat[,1]
d1$pm <- fmmat[,2]

## second calvert sampling
calvert2 <- read.csv("C:/Users/FABS/Documents/Ubpgebia/data/Upogebia+Orthione_20180814_Ben+Matt.csv")
# select and rename columns
d2 <- calvert2 %>% select( date=date.collected, method, sex=sex.shrimp, len.cara=length.carapace, len.total=length.total, 
                           pf=parasite.female, pm=parasite.male )
d2$site = "Pruth Bay"
d2$location = "Calvert Island"
d2$site.type = "Reference Site"
d2$date <- ymd( d2$date )

## third Calvert sampling
calvert3 <- read.csv("C:/Users/FABS/Documents/Ubpgebia/data/201901_Upogebia_Ben+Matt.csv", stringsAsFactors = FALSE )
# select and rename columns
d3 <- calvert3 %>% select( date=date.collected, method, sex=sex.shrimp, len.cara=length.carapace, len.total=length.total, 
                           pf=parasite.female, pm=parasite.male )
d3$site = "Pruth Bay"
d3$location = "Calvert Island"
d3$site.type = "Reference Site"
d3$date <- ymd( d3$date )



# additional data on parasitism from shellfish farms and reference sites on Vancouver Island
vancouverisland<-read.csv("C:/Users/FABS/Documents/Ubpgebia/data/Additional Pugettensis Data Updated Sept 3.csv")
# select and rename columns
d4 <- vancouverisland %>% select( site=Site, site.type=Site.Type, date=Date, location=Location, len.total=Size.mm., 
                        parasite.sex=ParasiteType   )
d4$date <- mdy( d4$date )
d4$len.cara <- 0.3 * d4$len.total
# convert parasite sex column to two presence-absence columns
fmmat <- matrix(NA,ncol=2,nrow=nrow(d4))
for( i in 1:nrow(d4) ){
  if( d4$parasite.sex[i]=="Both" )  fmmat[i,] <- c(1,1)
  if( d4$parasite.sex[i]=="Female" )  fmmat[i,] <- c(1,0)
  if( d4$parasite.sex[i]=="" )  fmmat[i,] <- c(0,0)
}
d4$pf <- fmmat[,1]
d4$pm <- fmmat[,2]


## ----------------------------------------------------------------------------------------------------------------------
## Combine these datasets
# join the first Calvert Sites
d12  <- full_join( d1, d2 )
d123 <- full_join( d12, d3 ) 
# join the Vancouver Island sites
d    <- full_join( d123, d4 )

# make date a factor 
d$date <- factor(d$date)
# make a factor to show whether parasites measured in field or lab
d$workspace <- ifelse( is.na(d$sex), "field", "lab" )
# change order of locations
d$location <- factor( d$location, levels=c("Calvert Island","Baynes Sound") )
# add a single parasite category (0=no parasite, 1=female only, 2=male+female
d$parasite.category <- d$pf + d$pm
#  separate calvert data
calvert <- d[d$location=="Calvert Island",]
# rename Shellfish farms and Reference sites
d$site.type[ d$site.type=="Shellfish Farm"] <- "Shellfish farm"
d$site.type[ d$site.type=="Reference Site"] <- "No aquaculture"
## ----------------------------------------------------------------------------------------------------------------------



# are all infestations represented by female parasites alone?
with( d, pf >= pm )
mistake <- which(with( d, pf < pm ))
# this is likely a mistake. Immature females look like males, so code this as a female, NOT a male
d$pf[ mistake ] <- 1
d$pm[ mistake ] <- 0
# 



### General Rates
ggplot(d, aes(as.factor(pf))) + geom_bar() + 
  facet_grid( location~site.type ) +
  xlab("") + theme_classic(base_size =15 )
# express as frequency
d.stats <- d %>%
  group_by(location,site.type, date) %>%
  mutate( n = n() ) %>%
  group_by(location,site.type, date, pf) %>%
  summarise(p = n()/n[1], n=n() ) %>%
  group_by(location,site.type,date) %>%
  mutate( total=sum(n)) %>%
  filter(pf==1)
# average over all data in each location and site.type
d.stats2 <- d %>%
  group_by(location,site.type) %>%
  mutate( n = n() ) %>%
  group_by(location,site.type,  pf) %>%
  summarise(p = n()/n[1], n=n() ) %>%
  group_by(location,site.type) %>%
  mutate( total=sum(n)) %>%
  filter(pf==1)

windows(3.5,3.5)
ggplot(d.stats, aes(x=location,y=p,col=site.type)) + geom_point(pch=1,size=3.5) +
  geom_text( label=as.character(d.stats$total), hjust=c(1.5,1.5,-0.75,1.5,2,1.5,1.5,1.5,2), show.legend=FALSE ) +
  geom_point( data=d.stats2,  size=3 ) +
  scale_color_manual(values=c('black', 'gray48')) +
  theme( legend.position = "top", legend.title = element_blank(), legend.text = element_text(size=10), legend.spacing.x = unit(0.2, 'cm') ) +
  guides( text = guide_legend(override.aes = list(text = NULL)) ) +
  xlab("Region") + ylab("Parasite prevalence")

with( d, sum(pf)/length(pf) )
with( d123, sum(pf)/length(pf) )
with( d123[d123$len.cara>=11.4,], sum(pf)/length(pf) )
by( d, d$site.type, function(z) sum(z$pf)/length(z$pf) )
by( d[d$location=="Calvert Island",], d$date[d$location=="Calvert Island"], function(z) c( sum(z$pf)/length(z$pf), length(z$pf) ) )
range( d$len.cara[d$location=="Calvert Island" & d$pf==1], na.rm=TRUE )
by( d[d$location=="Calvert Island" & d$len.cara>=11.4,], d$date[d$location=="Calvert Island" & d$len.cara>=11.4], function(z) c( sum(z$pf)/length(z$pf), length(z$pf) ) )
by( calvert, calvert$workspace, function(z) c( sum(z$pf)/length(z$pf), length(z$pf) ) )


#### Rate by size
# boxplots
ggplot(data= d, aes(x=pf, y=len.cara, fill=factor(pf) ))  +geom_boxplot()+ theme_classic(base_size =15 )+xlab("")
ggplot(data= d, aes(x=pf, y=len.cara, fill=factor(pf) ))  +geom_violin()+ theme_classic(base_size =15 )+xlab("")
# from Lemay
# stacked bar plots
calvert$bin <- cut( calvert$len.cara, seq(4,27,by=1) )
ptab <- with(calvert, table(pf,bin))
barplot(ptab)
pdf <- data.frame(ptab)
pdf$lower <- gsub("[(]","",pdf$bin)
pdf$lower <- as.numeric( gsub(",.*","",pdf$lower) )
stack <- ggplot( pdf, aes(x=lower,y=Freq,fill=pf)) + geom_bar(stat="identity",col='black') +
  xlab( "Host carapace length (mm)" ) + ylab("\nFrequency") +
  scale_fill_manual(name="Parasite",values=c("gray75","gray25"), 
                    breaks=c(0,1), labels=c("absent","present") ) +
  theme_classic() +
  theme(legend.justification=c(1,1), legend.position=c(1,1)) 
stack
# anova style of parasitized category with test
size.cat <- ggplot( calvert, aes(x=as.factor(parasite.category),y=len.cara)) + #geom_point(col='slateblue',alpha=0.5) +
  stat_summary( fun.data = "mean_cl_boot", colour = "black", size = 1.5)  + 
  annotate( "text", label = "A", x = 1, y = 15.25, size=5) +
  annotate( "text", label = "B", x = 2, y = 18.75, size=5) +
  annotate( "text", label = "B", x = 3, y = 21, size=5) +
  theme_classic() + 
  ylab( "Host carapace\nlength (mm)") + xlab( "Parasite load" )
calvert$parasite.category <- factor(calvert$parasite.category)
a1 <- aov( len.cara ~ (parasite.category), data=calvert )
anova(lm( len.cara ~ (parasite.category), data=calvert ))
summary(lm( len.cara ~ (parasite.category), data=calvert ))
TukeyHSD( a1 )
windows( 4.5, 4.5 )
plot_grid( stack, size.cat, ncol=1, labels = "AUTO", align="b" )


# prevalence by size
ggplot( data=d, aes(x=len.cara,y=pf,fill=location) ) +
  geom_vline( linetype=2,xintercept=12 ) +
  geom_point( col='slateblue', size = 3, alpha=0.3 ) +
  geom_smooth( method='glm', method.args = list(family=binomial),
               col = 'black', size=0.5 ) +
  ylab( "Parasite prevalence" ) + xlab( "Shrimp carapace length (mm)" ) +
  theme_classic() + xlim(c(0,max(d$len.cara)))


# prevalence by size and date at Calvert only
labels <- format( sort(unique(d123$date)), "%d %b %Y"  )
names(labels) <- sort(unique(d123$date))
windows(4,4)
ggplot( data=d123, aes(x=len.cara,y=pf) ) + facet_wrap(~date, labeller=labeller(date=labels)) +
  geom_vline( linetype=2,xintercept=12 ) +
  geom_smooth( method='glm', method.args = list(family=quasibinomial),
               col = 'black', size=0.5 ) +
  geom_point( col='slateblue', size = 3, alpha=0.3 ) +
  ylab( "Parasite prevalence" ) + xlab( "Shrimp carapace length (mm)" ) +
  theme_classic() + xlim(c(0,max(d$len.cara)))

# 2018-07-29 sampling was on eastern side of Pruth Bay (watershed 626, soft sed site PBE)
# sampling done by rolling boulders and hand-retrieving shrimp

# prevalence by size and sex at Calvert only
ggplot( data=d123[ !is.na(d123$sex), ], aes(x=len.cara,y=pf,fill=sex, col=sex, group=sex) ) + facet_wrap(~sex) +
  geom_vline( linetype=2,xintercept=12 ) +
  geom_point( size = 3, alpha=0.3 ) +
  geom_smooth( method='glm', method.args = list(family=quasibinomial), size=0.5, col='black' ) +
  ylab( "Parasite prevalence" ) + xlab( "Shrimp carapace length (mm)" ) 
by( d2, d2$sex, function(z) sum(z$pf)/length(z$pf) )
by( d3, d3$sex, function(z) sum(z$pf)/length(z$pf) )
by( d123, d123$sex, function(z) sum(z$pf)/length(z$pf) )



dsex <- d123 %>% 
  filter( !is.na(sex) ) %>% 
  mutate( sex=factor(sex) )

glmer.sex <- glmer( pf ~ sex + (1|date), data=dsex, family=binomial )
summary(glmer.sex)
glm.sex <- glm( pf ~ sex + date, data=dsex, family=quasibinomial )
summary(glm.sex)


##### MODELS

# model comparing parasite presence on Calvert and Vancouver Island
d$site.type <- factor(d$site.type)
d$baynes.com <- factor(with(d, paste(location,site.type) ))
d$baynes.com <- relevel( d$baynes.com, ref = "Calvert Island No aquaculture" )
glmer.site <- glmer( pf ~ baynes.com + (1|date), data=d, family='binomial' )
summary(glmer.site)

# a model of prevalence by size
glm1 <- glm( pf ~ len.cara, data=d123, family=binomial )
summary( glm1 )
# predict the extremes of body size in the dataset
psych::logistic( predict( glm1, newdata = list(len.cara=range(d$len.cara,na.rm=TRUE) ) ) )

# does prevalnce differ by collection date, location, or site type?
glm2 <- glm( pf ~ len.cara+date, data=d123, family=binomial )
summary(glm2)
glm2.5 <- glm( pf ~ len.cara+factor(date), data=d[d$location=="Calvert Island",], family=binomial )
summary(glm2.5)
glm3 <- glm( pf ~ len.cara+location, data=d, family=binomial )
summary(glm3)
glm4 <- glm( pf ~ len.cara+site.type, data=d, family=binomial )
summary(glm4)
glm4 <- glm( pf ~ len.cara+site.type, data=d, family=binomial )
summary(glm4)

# date matters, how does it influence parameter estimates when treated as a random effect?
glmer1 <- glmer( pf ~ len.cara + (1|date), data=calvert, family=binomial )
summary(glmer1)
glmer0 <- glmer( pf ~ 1 +(1|date), data=calvert, family=binomial )
anova( glmer0, glmer1 )
AICctab( glmer0, glmer1, nobs=nrow(calvert) )
MuMIn::r.squaredGLMM( glmer1 )
# glmer1 <- glmer( pf ~ as.factor(sex) + (1|date), data=d123, family=binomial )
summary(glmer1)
nd <- with( d123, expand.grid(len.cara=range(d123$len.cara,na.rm=TRUE) , date=unique(date) ) ) 
psych::logistic( predict( glmer1, newdata = nd ))
glmer2 <- glmer( pf ~ len.cara + (1|date), data=d123[d123$len.cara>=11.4,], family=binomial )
summary(glmer2)

# no "females" smaller that 12mm carapace length. Maybe filter to larger sizes (>15)
# we determined sex based on pleopods


### ----------------------------------------------------------------------------------------------------------------------
# MODEL PREDICTIONS
# Set up prediction frame
pred.frame <- with(calvert, expand.grid(len.cara=unique(len.cara))) # needs to be same as model fitted
X <- model.matrix(~len.cara,data=pred.frame) # calculate model matrix (formula needs to be same as the model fitted)
pred <- data.frame(pred.frame,mean=(X%*%fixef(glmer1)))  # these are the point predictions
# Calculate variance between observations within each site (for each combination of fixed effects)
V <- vcov(glmer1)
pred.var1 <- diag(X %*% V %*% t(X)) # XVX^T
# diag(X %*% tcrossprod(V,X)) # just another way

# Attach to dataframe, calculate standard errors and confidence intervals (using 1.96*sigma may be anti-conservative...)
predictions <- data.frame(pred,pred.se1=sqrt(pred.var1))
predictions <- with(predictions, data.frame(predictions,
                                            mean=logistic(mean),
                                            p1.lo = logistic(mean-1.96*pred.se1),
                                            p1.hi = logistic(mean+1.96*pred.se1 )) )

# plot 
len.pred <- ggplot( data=calvert, aes(x=len.cara,y=pf) ) +
  geom_vline( linetype=3,xintercept=12 ) +
  geom_point( col='slategrey', size = 3, alpha=0.3 ) +
  geom_line(  data=predictions, aes(x=len.cara,y=mean.1) ) +
  geom_line(  data=predictions, aes(x=len.cara,y=p1.lo), lty=2 ) +
  geom_line(  data=predictions, aes(x=len.cara,y=p1.hi), lty=2 ) +
  ylab( "Parasite\npresence" ) + xlab( "Host carapace length (mm)" ) +
  theme_classic() +
  xlim(c(0,max(d$len.cara))) +
  scale_y_continuous( breaks = c(0,1))

# base plot
# reorder prediction
# predictions <- predictions %>% arrange( len.cara )
# plot( pf~len.cara, data=calvert )
# lines( mean.1~len.cara, data=predictions )
# lines( p1.lo~len.cara, data=predictions, lty=2 )
# lines( p1.hi~len.cara, data=predictions, lty=2 )
plot_grid( stack, size.cat, len.pred, ncol=1, align="hv", labels="auto")

windows(5,5)
bottomrow <- plot_grid( size.cat, len.pred, ncol=2, align="hv", labels=c("b","c"))
plot_grid( stack, bottomrow, ncol=1,  labels = c("a",""))


### ----------------------------------------------------------------------------------------------------------------------
#ADDITIONAL PLOTS

#plot of prevalence by sample event on Calvert Island

d.calver.month <- d %>%
  filter(location == "Calvert Island") %>%
  group_by(date) %>%
  summarize(n = length(date),
            inf = sum(pf)) %>% 
  mutate(inf_rate = inf/n)

ggplot( data=d.calver.month, aes(x=date,y=inf_rate) ) +
  geom_col() +
  ylab( "Parasite prevalence" ) + xlab( "Date" ) +
  theme_classic() 



