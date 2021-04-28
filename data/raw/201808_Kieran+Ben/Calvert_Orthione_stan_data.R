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
library(psych) # for logistic function
## read data sets --------------------------------------------------------------------------------------------------------
# first calvert sampling 2018
calvert1<-read.csv("PugettensisCalvertAugV1D1.csv")
# select and rename columns
d1 <- calvert1 %>% select( site=Site, date=Date, len.total=Size.mm., parasite.sex=ParasiteType )
d1$location = "Calvert Island"
d1$site.type = "Reference Site"
d1$date <- dmy( d1$date )
d1$len.cara <- 0.3 * d1$len.total
# convert parasite sex column to two presence-absence columns
fmmat <- matrix(NA,ncol=2,nrow=nrow(d1))
for( i in 1:nrow(d1) ){
  if( d1$parasite.sex[i]=="Both" )  fmmat[i,] <- c(1,1)
  if( d1$parasite.sex[i]=="Female" )  fmmat[i,] <- c(1,0)
  if( d1$parasite.sex[i]=="None" )  fmmat[i,] <- c(0,0)
}
d1$pf <- fmmat[,1]
d1$pm <- fmmat[,2]

# second calvert sampling
calvert2 <- read.csv("Upogebia+Orthione_20180814_Ben+Matt.csv")
# select and rename columns
d2 <- calvert2 %>% select( date=date.collected, method, sex=sex.shrimp, len.cara=length.carapace, len.total=length.total, 
                           pf=parasite.female, pm=parasite.male )
d2$site = "Pruth Bay"
d2$location = "Calvert Island"
d2$site.type = "Reference Site"
d2$date <- ymd( d2$date )

# additional data on parasitism from shellfish farms and reference sites on Vancouver Island
vancouverisland<-read.csv("Additional Pugettensis Data Updated Sept 3.csv")
# select and rename columns
d3 <- vancouverisland %>% select( site=Site, site.type=Site.Type, date=Date, location=Location, len.total=Size.mm., 
                        parasite.sex=ParasiteType   )
d3$date <- mdy( d3$date )
d3$len.cara <- 0.3 * d3$len.total
# convert parasite sex column to two presence-absence columns
fmmat <- matrix(NA,ncol=2,nrow=nrow(d3))
for( i in 1:nrow(d3) ){
  if( d3$parasite.sex[i]=="Both" )  fmmat[i,] <- c(1,1)
  if( d3$parasite.sex[i]=="Female" )  fmmat[i,] <- c(1,0)
  if( d3$parasite.sex[i]=="" )  fmmat[i,] <- c(0,0)
}
d3$pf <- fmmat[,1]
d3$pm <- fmmat[,2]


## ----------------------------------------------------------------------------------------------------------------------
## Combine these datasets
# join the Calvert Sites
d12 <- full_join( d1, d2 )
# join the Vancouver Island sites
d   <- full_join( d12, d3 ) 

# make date a factor 
d$date <- factor(d$date)
# make a factor to show whether parasites measured in field or lab
d$workspace <- ifelse( is.na(d$sex), "field", "lab" )

# are all infestations represented by female parasites alone?
with( d, pf >= pm )
mistake <- which(with( d, pf < pm ))
# this is likely a mistake. Immature females look like males
d$pf[ mistake ] <- 1
d$pm[ mistake ] <- 0
# 

#  separate calvert data
calvert <- d[d$location=="Calvert Island",]
## ----------------------------------------------------------------------------------------------------------------------

# write relevant data to disk for downstream analysis in STAN
calvert.stan <- calvert %>%
  select( infected=pf, size=len.cara )

N <- 167
M <- 1
X <- calvert.stan$size 
y <- calvert.stan$infected


