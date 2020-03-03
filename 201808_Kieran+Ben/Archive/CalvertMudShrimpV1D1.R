
library(ggplot2)

library(cowplot)

theme_black = function(base_size = 12, base_family = "") {
  
  theme_grey(base_size = base_size, base_family = base_family) %+replace%
    
    theme(
      # Specify axis options
      axis.line = element_blank(),  
      axis.text.x = element_text(size = base_size*0.8, color = "white", lineheight = 0.9),  
      axis.text.y = element_text(size = base_size*0.8, color = "white", lineheight = 0.9),  
      axis.ticks = element_line(color = "white", size  =  0.2),  
      axis.title.x = element_text(size = base_size, color = "white", margin = margin(0, 10, 0, 0)),  
      axis.title.y = element_text(size = base_size, color = "white", angle = 90, margin = margin(0, 10, 0, 0)),  
      axis.ticks.length = unit(0.3, "lines"),   
      # Specify legend options
      legend.background = element_rect(color = NA, fill = "black"),  
      legend.key = element_rect(color = "white",  fill = "black"),  
      legend.key.size = unit(1.2, "lines"),  
      legend.key.height = NULL,  
      legend.key.width = NULL,      
      legend.text = element_text(size = base_size*0.8, color = "white"),  
      legend.title = element_text(size = base_size*0.8, face = "bold", hjust = 0, color = "white"),  
      legend.position = "right",  
      legend.text.align = NULL,  
      legend.title.align = NULL,  
      legend.direction = "vertical",  
      legend.box = NULL, 
      # Specify panel options
      panel.background = element_rect(fill = "black", color  =  NA),  
      panel.border = element_rect(fill = NA, color = "white"),  
      panel.grid.major = element_line(color = "grey35"),  
      panel.grid.minor = element_line(color = "grey20"),  
      panel.margin = unit(0.5, "lines"),   
      # Specify facetting options
      strip.background = element_rect(fill = "grey30", color = "grey10"),  
      strip.text.x = element_text(size = base_size*0.8, color = "white"),  
      strip.text.y = element_text(size = base_size*0.8, color = "white",angle = -90),  
      # Specify plot options
      plot.background = element_rect(color = "black", fill = "black"),  
      plot.title = element_text(size = base_size*1.2, color = "white"),  
      plot.margin = unit(rep(1, 4), "lines")
      
    )
  
}

Pugettensis<-read.csv("PugettensisCalvertAugV1D1.csv")
Pugettensis
names(Pugettensis)
levels(Pugettensis$Date)
levels(Pugettensis$ParasiteType)
levels(Pugettensis$Parasitized)
levels(Pugettensis$IndividualID)
summary(Pugettensis$ParasiteType)
hist(Pugettensis$Size.mm.)

Additional<-read.csv("Additional Pugettensis Data.csv")
names(Additional)

# add presence absence for parasites
Pugettensis$presence <- as.numeric(Pugettensis$Parasitized)-1
Additional$presence <- as.numeric(Additional$Parasitized)-1

library(ggplot2)
library(tidyverse)

### General Rates
ggplot(Pugettensis, aes(Parasitized, fill=Parasitized))+geom_bar()+theme_classic(base_size =15 )+xlab("")
ggplot(Additional, aes(Parasitized, fill=Parasitized))+geom_bar()+theme_classic(base_size =15 )+xlab("")

ggplot(Pugettensis, aes(ParasiteType, fill=ParasiteType))+geom_bar()+theme_classic(base_size =15 )+xlab("")

with( Pugettensis, sum(presence)/length(presence) )
with( Additional, sum(presence)/length(presence) )


#### Rate by size
ggplot(data= Pugettensis, aes(x=Parasitized, y=Size.mm., fill=Parasitized))  +geom_boxplot()+ theme_classic(base_size =15 )+xlab("")#+ ggtitle("") + ylab("")+ xlab("") +scale_colour_manual(values=c("#666666","#CC0000","#0000CC"))  +theme(plot.title = element_text(size = rel(1)))+theme(axis.text.x = element_text(angle = 5, hjust = 0.5,colour = "black",size=15))+guides(colour=FALSE, fill=FALSE)#+ facet_grid(SoundSourceStimulis~.)
ggplot(data= Pugettensis, aes(x=Parasitized, y=Size.mm., fill=Parasitized))  +geom_violin()+ theme_classic(base_size =15 )+xlab("")#+ ggtitle("") + ylab("")+ xlab("") +scale_colour_manual(values=c("#666666","#CC0000","#0000CC"))  +theme(plot.title = element_text(size = rel(1)))+theme(axis.text.x = element_text(angle = 5, hjust = 0.5,colour = "black",size=15))+guides(colour=FALSE, fill=FALSE)#+ facet_grid(SoundSourceStimulis~.)
ggplot(data= Additional,  aes(x=Parasitized, y=Size.mm., fill=Parasitized))  +geom_violin()+ theme_classic(base_size =15 )+xlab("")#+ ggtitle("") + ylab("")+ xlab("") +scale_colour_manual(values=c("#666666","#CC0000","#0000CC"))  +theme(plot.title = element_text(size = rel(1)))+theme(axis.text.x = element_text(angle = 5, hjust = 0.5,colour = "black",size=15))+guides(colour=FALSE, fill=FALSE)#+ facet_grid(SoundSourceStimulis~.)

ggplot(data= Pugettensis, aes(x=ParasiteType, y=Size.mm., fill=ParasiteType))  +geom_boxplot()+ theme_classic(base_size =15 )+xlab("")#+ ggtitle("") + ylab("")+ xlab("") +scale_colour_manual(values=c("#666666","#CC0000","#0000CC"))  +theme(plot.title = element_text(size = rel(1)))+theme(axis.text.x = element_text(angle = 5, hjust = 0.5,colour = "black",size=15))+guides(colour=FALSE, fill=FALSE)#+ facet_grid(SoundSourceStimulis~.)


# prevalence by size

# windows(9,4)
# classic
size1 <- ggplot( data=Pugettensis, aes(x=Size.mm.,y=presence) ) +
  geom_vline( linetype=2,xintercept=40 ) +
  geom_point( col='slateblue', size = 3, alpha=0.3 ) +
  geom_smooth( method='glm', method.args = list(family=binomial),
               col = 'black', size=0.5 ) +
  ylab( "Parasite prevalence" ) + xlab( "Total shrimp length (mm)" ) +
  theme_classic() + xlim(c(0,80))
size2 <- ggplot( data=Additional, aes(x=Size.mm.,y=presence) ) +
  geom_vline( linetype=2,xintercept=40 ) +
  geom_jitter( aes(col=Site), size = 3, alpha=0.8, width=0.1, height=0.01 ) +
  geom_smooth( method='glm', method.args = list(family=binomial),
               col = 'black', size=0.5 ) +
  ylab( "Parasite prevalence" ) + xlab( "Total shrimp length (mm)" ) +
  theme_classic() + xlim(c(0,80))
plot_grid( size1, size2, rel_widths = c(1,1.5) )

# black out
ggplot( data=Pugettensis, aes(x=Size.mm.,y=presence) ) +
  geom_vline( linetype=2,xintercept=40 ) +
  geom_point( col='chartreuse', size = 3, alpha=0.3 ) +
  geom_smooth( method='glm', method.args = list(family=binomial),
               col = 'white', size=0.5 ) +
  ylab( "Parasite prevalence" ) + xlab( "Total shrimp length (mm)" ) 
  theme_black()
ggplot( data=Additional, aes(x=Size.mm.,y=presence) ) +
  geom_vline( linetype=2,xintercept=40 ) +
  geom_point( col='chartreuse', size = 3, alpha=0.3 ) +
  geom_smooth( method='glm', method.args = list(family=binomial),
               col = 'white', size=0.5 ) +
  ylab( "Parasite prevalence" ) + xlab( "Total shrimp length (mm)" ) +
  theme_black()

# the model
glm1 <- glm( presence ~ Size.mm., data=Pugettensis, family=binomial )
summary( glm1 )
# predict the extremes of body size in the dataset
psych::logistic( predict( glm1, newdata = list(Size.mm.=range(Pugettensis$Size.mm.) ) ) )


### Example Figure

A<-ggplot(Pugettensis, aes(Parasitized, fill=Parasitized))+geom_bar()+theme_classic(base_size =15 )+xlab("")
B<-ggplot(Pugettensis, aes(ParasiteType, fill=ParasiteType))+geom_bar()+theme_classic(base_size =15 )+xlab("")
C<-ggplot(data= Pugettensis, aes(x=Parasitized, y=Size.mm., fill=Parasitized))  +geom_boxplot()+ theme_classic(base_size =15 )+xlab("")#+ ggtitle("") + ylab("")+ xlab("") +scale_colour_manual(values=c("#666666","#CC0000","#0000CC"))  +theme(plot.title = element_text(size = rel(1)))+theme(axis.text.x = element_text(angle = 5, hjust = 0.5,colour = "black",size=15))+guides(colour=FALSE, fill=FALSE)#+ facet_grid(SoundSourceStimulis~.)
D<-ggplot(data= Pugettensis, aes(x=ParasiteType, y=Size.mm., fill=ParasiteType))  +geom_boxplot()+ theme_classic(base_size =15 )+xlab("")#+ ggtitle("") + ylab("")+ xlab("") +scale_colour_manual(values=c("#666666","#CC0000","#0000CC"))  +theme(plot.title = element_text(size = rel(1)))+theme(axis.text.x = element_text(angle = 5, hjust = 0.5,colour = "black",size=15))+guides(colour=FALSE, fill=FALSE)#+ facet_grid(SoundSourceStimulis~.)

plot_grid(A,B,C,D, labels=c("","","",""), ncol = 2, nrow =2,label_size=18)


