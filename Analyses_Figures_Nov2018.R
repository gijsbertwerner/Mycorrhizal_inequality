#Mycorrhizal Inequality Paper - Data Analysis
#(C) Gijsbert Werner, Department of Zoology, University of Oxford, gijsbert.werner@zoo.ox.ac.uk
#November 2018

##This script performs all the statistical analyses presented in Whiteside et al. 2018 and generates all its figures.
##Please make sure that data files are in the same relative paths as here, or change paths to where you have saved them.

##Running the associated Markdown-script Analyses_Figures_Report_V2.Rmd produces a markdown report of all our results. 
##This script was executed in R 3.4.4.

# Load Packages  -------------------------------------------------------

library(ggplot2)
library(reshape2)
library(nlme)
library(lme4)
library(car)
library(ggpubr)
library(gridExtra)
library(boot)
library(dplyr)
stderr <- function(x) sqrt(var(x,na.rm=TRUE)/length(na.omit(x)))
dir.create("./Models")
dir.create("./Figures")

# Load Data ---------------------------------------------------------------

ineq<-read.csv("./Data/Ineq_full_fungal_experimental_data.csv",as.is = T)
head(ineq)
lapply(ineq,class)
ineq$ineq_treatment<-factor(ineq$ineq_treatment,levels = c("None","Medium","High"))

ineq_nut_level_control<-read.csv("./Data/Ineq_nut_level_control_experimental_data.csv",as.is=T)
head(ineq_nut_level_control)
lapply(ineq_nut_level_control,class)

spore_growth<-read.csv("./Data/Spore_growth_data.csv")
head(spore_growth)
lapply(spore_growth,class)

in_vitro_dat<-read.csv("./Data/In_vitro_colour_control_data.csv",as.is=T)
head(in_vitro_dat)
lapply(in_vitro_dat,class)

colon<-read.csv("./Data/Colonisation_Data.csv",as.is=T)
head(colon)
lapply(colon,class)

whole_plant_dat<-read.csv("./Data/Whole_plant_data.csv",as.is=T)
head(whole_plant_dat)
lapply(whole_plant_dat,class)

tox_dat<-read.csv("./Data/Toxicity_data.csv",as.is=T)
head(tox_dat)
lapply(tox_dat,class)

# Figure 2A ---------------------------------------------------------------

# How is the overall, plate-level, P-transfer to the roots affected by plate inequality level?

##Data formatting
# This is a question at the root-compartment level.
ineq_root_comp <- ineq %>% filter(comp=="Root")

#We are interested in P tranferred per unit of root. To get this, we just sum the two source compartments
ineq_root_comp$root_QD_per_mg_root_both_origin<-ineq_root_comp$root_QD_per_mg_root_rich_origin+ineq_root_comp$root_QD_per_mg_root_poor_origin
head(ineq_root_comp)

#Descriptives 

#Have a look at the data
ggstripchart(ineq_root_comp,x = "ineq_treatment",y = "root_QD_per_mg_root_both_origin",
             add=c("mean_se","violin"),color="ineq_treatment",add.params = list(color="darkgray"))
#Looks like a (weak) positive effect of inequality on total transfer. 

##Statistical Modelling

#Linear model
lm_ineq_root_QD_per_mg<-lm(root_QD_per_mg_root_both_origin~ineq_treatment,ineq_root_comp)
save(lm_ineq_root_QD_per_mg,file="./Models/lm_ineq_root_QD_per_mg")

#Visually check if assumptions of the model are met
residualPlots(lm_ineq_root_QD_per_mg)
qqPlot(lm_ineq_root_QD_per_mg)
#This isn't great, but not too bad either. Evaluate GLMs to see if we can get a better fit.

#GLMs

#Gaussian error distributions
glm_ineq_root_QD_per_mg_gaus_identity<-glm(root_QD_per_mg_root_both_origin~ineq_treatment,ineq_root_comp,family = gaussian(link = "identity"))
glm_ineq_root_QD_per_mg_gaus_log<-glm(root_QD_per_mg_root_both_origin~ineq_treatment,ineq_root_comp,family = gaussian(link = "log"))
glm_ineq_root_QD_per_mg_gaus_inverse<-glm(root_QD_per_mg_root_both_origin~ineq_treatment,ineq_root_comp,family = gaussian(link = "inverse"))
#Diagnostics
glm.diag.plots(glm_ineq_root_QD_per_mg_gaus_identity)
glm.diag.plots(glm_ineq_root_QD_per_mg_gaus_log)
glm.diag.plots(glm_ineq_root_QD_per_mg_gaus_inverse)
#Unsurprisingly, none of these really look better than the the linear model we considered before (the identity should be the exact same), consider Gamma-distributions

#Gamma error distribution
glm_ineq_root_QD_per_mg_gamma_inverse<-glm(root_QD_per_mg_root_both_origin~ineq_treatment,ineq_root_comp,family = Gamma(link = "inverse"))
glm_ineq_root_QD_per_mg_gamma_log<-glm(root_QD_per_mg_root_both_origin~ineq_treatment,ineq_root_comp,family = Gamma(link = "log"))
glm_ineq_root_QD_per_mg_gamma_identity<-glm(root_QD_per_mg_root_both_origin~ineq_treatment,ineq_root_comp,family = Gamma(link = "identity"))
#Diagnostics
glm.diag.plots(glm_ineq_root_QD_per_mg_gamma_inverse)
glm.diag.plots(glm_ineq_root_QD_per_mg_gamma_log)
glm.diag.plots(glm_ineq_root_QD_per_mg_gamma_identity)
#The gamma distribution with identity link looks very good. Analyse it further. 

#Further model analysis
summary(glm_ineq_root_QD_per_mg_gamma_identity) #High inequality signficantly more transfer than the intercept (i.e. no inequality)
Anova(glm_ineq_root_QD_per_mg_gamma_identity,test.statistic = "F")

#Save the models
save(glm_ineq_root_QD_per_mg_gamma_identity,file="./Models/glm_ineq_root_QD_per_mg_gamma_identity")

###Figure

#Calculate the plot data. Now use the two source compartments separately for visualisation purposes
PlotData_Fig2A<-ineq_root_comp %>% 
  melt(id.vars="ineq_treatment",
       measure.vars=c("root_QD_per_mg_root_poor_origin","root_QD_per_mg_root_rich_origin","root_QD_per_mg_root_both_origin")) %>% 
  group_by(ineq_treatment,variable) %>%
  summarise(mean_root_QD=mean(value,na.rm = T),
            se_root_QD=stderr(value))

PlotData_Fig2A$variable<-factor(PlotData_Fig2A$variable,levels = c("root_QD_per_mg_root_rich_origin",
                                                                   "root_QD_per_mg_root_poor_origin",
                                                                   "root_QD_per_mg_root_both_origin"))

#Plot the figure
Fig2a<-ggplot()+
  geom_col(aes(x=ineq_treatment,y=mean_root_QD,fill=variable),
           data=PlotData_Fig2A %>% filter(variable %in% c("root_QD_per_mg_root_rich_origin","root_QD_per_mg_root_poor_origin")),
           position = position_stack(reverse = F))+
  geom_errorbar(aes(x=ineq_treatment,ymin=mean_root_QD-se_root_QD,ymax=mean_root_QD+se_root_QD),
                data=PlotData_Fig2A %>% filter(variable=="root_QD_per_mg_root_both_origin"),width=0.25)+
  scale_fill_manual(values = c("#00BFC4","#F8766D"),
                    name="Nutrient\nCompartment",
                    breaks=c("root_QD_per_mg_root_rich_origin","root_QD_per_mg_root_poor_origin"),
                    labels=c("Rich","Poor"))+
  labs(x="Level of Inequality",y="P transfer to roots (nmol QD-apatite per mg root)")+
  #scale_y_continuous(limits = c(0,0.08),breaks = seq(from=0.01,to=0.07,by=0.01),labels = c("","0.02","","0.04","","0.06",""))+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),
        axis.line = element_line(colour="black"))+
  geom_text(aes(x=3,y=0.0012,label="*",size=30))
Fig2a

png("./Figures/Figure2a.png",res = 300,width = 5,height = 5,units = "in")
Fig2a
dev.off()

#Create a table with the percentages from the origin compartments (Table S1)
ineq_root_comp<-ineq_root_comp %>% mutate(perc_rich_origin=(100*root_QD_per_mg_root_rich_origin/root_QD_per_mg_root_both_origin),
                                           perc_poor_origin=(100*root_QD_per_mg_root_poor_origin/root_QD_per_mg_root_both_origin))
ineq_transfer_stats_table <- ineq_root_comp %>%
  group_by(ineq_treatment) %>% 
  summarise(mean_perc_rich_origin=mean(perc_rich_origin),
            mean_perc_poor_origin=mean(perc_poor_origin),
            se_perc_origin=stderr(perc_poor_origin))
ineq_transfer_stats_table
save(ineq_transfer_stats_table,file="./Models/ineq_transfer_stats_table")

#Statistical test, to compare if without inequality both compartments transfer equally well. 
t.test(ineq_root_comp %>% filter(ineq_treatment=="None") %>% select(perc_poor_origin),mu=50)

#Save the analysed data frame
save(ineq_root_comp,file="./Models/ineq_root_comp")


# Figure 2B ---------------------------------------------------------------

#How is the P retention within the hyphae affected by inequality level? 

##Data formatting 
#We now need to get separate columns for all the hyphal retained rich vs poor origin QDs, regardless of where they are. 
ineq_hyphal_retained_molten <-ineq  %>% select(plate,comp,hyphal_QD_total,hyphal_QD_total_from_other_hypal,hyphal_QD_in_root_from_poor_total,hyphal_QD_in_root_from_rich_total,hyphal_bm_mg)  %>%
  melt(id.vars=c("plate","comp"),variable.name="measurement",value.name="QD")
#Sum the right compartments
ineq_hyphal_retained_subdiv_sum <- ineq_hyphal_retained_molten %>% group_by(plate) %>%
  summarise(poor_QD_sum = QD[measurement=="hyphal_QD_total"&comp=="Poor"]+
              QD[measurement=="hyphal_QD_total_from_other_hypal"&comp=="Rich"]+
              QD[measurement=="hyphal_QD_in_root_from_poor_total"&comp=="Root"],
            rich_QD_sum = QD[measurement=="hyphal_QD_total"&comp=="Rich"]+
              QD[measurement=="hyphal_QD_total_from_other_hypal"&comp=="Poor"]+
              QD[measurement=="hyphal_QD_in_root_from_rich_total"&comp=="Root"],
            hyphal_bm_mg_total=sum(QD[measurement=="hyphal_bm_mg"]))
head(ineq_hyphal_retained_subdiv_sum)
#Combine with marker for inequality level, and sum to get the total level
ineq_hyphal_retained<-inner_join(x= ineq_hyphal_retained_subdiv_sum,
                                 y= ineq %>% filter(comp=="Poor") %>% select(plate,ineq_treatment))
ineq_hyphal_retained$hyphal_retention<-ineq_hyphal_retained$poor_QD_sum+ineq_hyphal_retained$rich_QD_sum
head(ineq_hyphal_retained)

##Descriptives
#Plot the data per group to evaluate distribution
ggstripchart(ineq_hyphal_retained,x = "ineq_treatment",y = "hyphal_retention",
             add="mean_se",color="ineq_treatment") 
#Looks like a clear effect of inequality on QD retention within the hyphae

##Statistical modelling
lm_ineq_hyphal_retention<-lm(hyphal_retention~ineq_treatment,ineq_hyphal_retained)
save(lm_ineq_hyphal_retention,file="./Models/lm_ineq_hyphal_retention")

#Visually check if assumptions of the model are met
residualPlot(lm_ineq_hyphal_retention)
qqPlot(lm_ineq_hyphal_retention)
#Ok, this is not surpising looking at the data, but the model assumptions are not at all met. 

#GLMs

#Gamma error distribution
#I will now start exploring Gamma error distributions
glm_ineq_hyphal_retention_gamma_identity<-glm(hyphal_retention~ineq_treatment,ineq_hyphal_retained,family = Gamma(link = "identity"))
glm_ineq_hyphal_retention_gamma_log<-glm(hyphal_retention~ineq_treatment,ineq_hyphal_retained,family = Gamma(link = "log"))
glm_ineq_hyphal_retention_gamma_inverse<-glm(hyphal_retention~ineq_treatment,ineq_hyphal_retained,family = Gamma(link = "inverse"))
#Diagnostics
glm.diag.plots(glm_ineq_hyphal_retention_gamma_identity)
glm.diag.plots(glm_ineq_hyphal_retention_gamma_log)
glm.diag.plots(glm_ineq_hyphal_retention_gamma_inverse)
#The gamma distribution with identity link looks reasonably good. 
save(glm_ineq_hyphal_retention_gamma_identity,file="./Model/glm_ineq_hyphal_retention_gamma_identity")
#There is one relatively influential point in the Cook statistics. Depending on the p-values, our inference might be fairly robust against that.

#Analyse the model 
summary(glm_ineq_hyphal_retention_gamma_identity)
Anova(glm_ineq_hyphal_retention_gamma_identity,test.statistic = "F")
#Very strong effect of inequality on retention. 

#Save the models
save(glm_ineq_hyphal_retention_gamma_identity,file="./Models/glm_ineq_hyphal_retention_gamma_identity")

##Secondary analysis, this time on a per mg basis. 
#Plot the data per group to evaluate distribution
ineq_hyphal_retained$hyphal_retention_per_hyphal_mg<-ineq_hyphal_retained$hyphal_retention/ineq_hyphal_retained$hyphal_bm_mg_total
ggstripchart(ineq_hyphal_retained,x = "ineq_treatment",y = "hyphal_retention_per_hyphal_mg",
             add="mean_se",color="ineq_treatment") 
#Looks like a clear effect of inequality on QD retention within the hyphae

lm_ineq_hyphal_retention_per_hyphal_mg<-lm(hyphal_retention_per_hyphal_mg~ineq_treatment,ineq_hyphal_retained)
save(lm_ineq_hyphal_retention_per_hyphal_mg,file="./Models/lm_ineq_hyphal_retention_per_hyphal_mg")

#Visually check if assumptions of the model are met
residualPlot(lm_ineq_hyphal_retention_per_hyphal_mg)
qqPlot(lm_ineq_hyphal_retention_per_hyphal_mg)

#Gamma error distribution
#I will again start exploring Gamma error distributions
glm_ineq_hyphal_retention_per_hyphal_mg_gamma_identity<-glm(hyphal_retention_per_hyphal_mg~ineq_treatment,ineq_hyphal_retained,family = Gamma(link = "identity"))
glm_ineq_hyphal_retention_per_hyphal_mg_gamma_log<-glm(hyphal_retention_per_hyphal_mg~ineq_treatment,ineq_hyphal_retained,family = Gamma(link = "log"))
glm_ineq_hyphal_retention_per_hyphal_mg_gamma_inverse<-glm(hyphal_retention_per_hyphal_mg~ineq_treatment,ineq_hyphal_retained,family = Gamma(link = "inverse"))
#Diagnostics
glm.diag.plots(glm_ineq_hyphal_retention_per_hyphal_mg_gamma_identity)
glm.diag.plots(glm_ineq_hyphal_retention_per_hyphal_mg_gamma_log)
glm.diag.plots(glm_ineq_hyphal_retention_per_hyphal_mg_gamma_inverse)
#The gamma distribution with identity link looks reasonably good. 
#There is one relatively influential point in the Cook statistics. Depending on the p-values, our inference might be fairly robust against that.
save(glm_ineq_hyphal_retention_per_hyphal_mg_gamma_identity,file="./Models/glm_ineq_hyphal_retention_per_hyphal_mg_gamma_identity")

#Analyse the model 
summary(glm_ineq_hyphal_retention_per_hyphal_mg_gamma_identity)
Anova(glm_ineq_hyphal_retention_per_hyphal_mg_gamma_identity,test.statistic = "F")
#Very strong effect of inequality on retention. 

##Visualisation

#Create the data frames for plotting Figure 2B
ineq_hyphal_retained_stats<-ineq_hyphal_retained %>% 
  group_by(ineq_treatment) %>%
  summarise(mean_retention=mean(hyphal_retention,na.rm=T),
            se_retention=stderr(hyphal_retention),
            mean_poor_retention=mean(poor_QD_sum,na.rm=T),
            se_poor_retention=stderr(poor_QD_sum),
            mean_rich_retention=mean(rich_QD_sum,na.rm=T),
            se_rich_retention=stderr(rich_QD_sum))
ineq_hyphal_retained_stats

ineq_hyphal_retained_stats_melt<-melt(as.data.frame(ineq_hyphal_retained_stats),id.vars="ineq_treatment",
                                measure.vars = c("mean_retention","mean_poor_retention",
                                                 "mean_rich_retention","se_retention"),
                                variable.name = "compartment",value.name="mean_total_retention")
ineq_hyphal_retained_stats_melt

ineq_hyphal_retained_stats_melt_fig<-ineq_hyphal_retained_stats_melt %>% 
  filter(compartment %in% c("mean_poor_retention","mean_rich_retention"))
ineq_hyphal_retained_stats_melt_fig$compartment<-c("poor","poor","poor","rich","rich","rich")
ineq_hyphal_retained_stats_melt_fig

Figure2B<-ggplot()+
  geom_col(aes(x=ineq_treatment,y=mean_total_retention,fill=compartment),
           data=ineq_hyphal_retained_stats_melt_fig,position = position_stack(reverse = TRUE))+
  geom_errorbar(aes(x=ineq_treatment,ymin=mean_retention-se_retention,
                    ymax=mean_retention+se_retention),data=ineq_hyphal_retained_stats,
                width=0.25)+
  scale_fill_discrete(name="Nutrient\nCompartment")+
  #theme(axis.text.x=element_text(angle = -45, hjust = 0))+
  ylab("Total P retained by hyphae (nmol QD-apatite)")+
  xlab("Level of Inequality")+
  scale_y_continuous(breaks = seq(from=0.01,to=0.12,by=0.01),labels = c("","0.02","","0.04","","0.06","","0.08","","0.10","","0.12"))+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),
        axis.line = element_line(colour="black"))+
  guides(fill=guide_legend(reverse=TRUE),colour=guide_legend(reverse=TRUE))+
  geom_text(aes(x=2,y=0.1,label="*",size=30))+
  geom_text(aes(x=3,y=0.03,label="***",size=30))
Figure2B

#Figure 2B plot
png("./Figures/Figure2B.png",res = 300,width = 5,height = 5,units = "in")
Figure2B
dev.off()

#Create a table with the percentages from the origin compartments 
ineq_hyphal_retained<-ineq_hyphal_retained %>% mutate(perc_rich_retained=(100*rich_QD_sum/hyphal_retention),
                                          perc_poor_retained=(100*poor_QD_sum/hyphal_retention))
ineq_hyphal_retained_table <- ineq_hyphal_retained %>%
  group_by(ineq_treatment) %>% 
  summarise(mean_perc_rich_retained=mean(perc_rich_retained,na.rm=T),
            mean_perc_poor_retained=mean(perc_poor_retained,na.rm=T),
            se_perc_retained=stderr(perc_poor_retained))
ineq_hyphal_retained_table
save(ineq_hyphal_retained_table,file="./Models/ineq_hyphal_retained_table")

#Statistical test, to compare if without inequality both compartments transfer equally well. 
t.test(as.data.frame(ineq_hyphal_retained %>% filter(ineq_treatment=="None") %>% select(perc_poor_retained)),
       mu=50)

#Save the data frame
save(ineq_hyphal_retained,file="./Models/ineq_hyphal_retained")

# Figure 2C ---------------------------------------------------------------

#How does compartment nutrient level affect vacuole area. 

#Data formatting. 
#We need to select only the hyphal compartments
ineq_vacuole<-ineq %>% dplyr::filter(comp %in% c("Poor","Rich")) 
head(ineq_vacuole)

#Add the numeric information about what the nutrient concentration per compartment is
ineq_vacuole<- ineq_vacuole %>% 
  mutate(comp_nutrients = case_when(
    ineq_treatment == "High" & comp == "Rich" ~ 90,
    ineq_treatment == "High" & comp == "Poor" ~ 10,
    ineq_treatment == "Medium" & comp == "Rich" ~ 70,
    ineq_treatment == "Medium" & comp == "Poor" ~ 30,
    ineq_treatment == "None" ~ 50)
  ) 
ineq_vacuole

##Descriptives
#Plot the data per group to evaluate distribution
ggstripchart(ineq_vacuole,x = "comp_nutrients",y = "percent_vac",
             add="mean_se",color="ineq_treatment") 

##Statistical modelling
perc_vac_poly_lme<-lme(percent_vac~poly(comp_nutrients,2),random = ~ 1|plate,
                       data = ineq_vacuole,na.action="na.omit")
#Save the mdoel 
save(perc_vac_poly_lme,file="./Models/perc_vac_poly_lme")

#Visually check if assumptions of the model are met
plot(perc_vac_poly_lme)
qqnorm(perc_vac_poly_lme$residuals)
qqline(perc_vac_poly_lme$residuals)

#What does the model tell us? 
summary(perc_vac_poly_lme)
Anova(perc_vac_poly_lme)

##Visualisation

#Print prediction line
polynomial_line<-predict(perc_vac_poly_lme,list(comp_nutrients=seq(1,100,1)),level=0)
predicted_poly_line<-rbind(cbind(seq(1,100,1),polynomial_line))
predicted_poly_line<-as.data.frame(predicted_poly_line,row.names = NULL)
names(predicted_poly_line)<-c("compartment_values","predictions")
predicted_poly_line

##The figure
ineq_vacuole_fig_stats<-ineq_vacuole %>% 
  group_by(comp_nutrients) %>%
  summarise(mean_perc_vac=mean(percent_vac,na.rm = T),
            se_perc_vac=stderr(percent_vac),
            n=n())
ineq_vacuole_fig_stats

Figure2c<-ggplot()+
  geom_pointrange(aes(x=comp_nutrients,ymin=mean_perc_vac-se_perc_vac,
                      ymax=mean_perc_vac+se_perc_vac,y=mean_perc_vac),
                  data=ineq_vacuole_fig_stats)+
  geom_line(aes(x=compartment_values,y=predictions),
            data=predicted_poly_line)+
  ylab("Relative surface area vacuoles (%)")+
  xlab("Nutrient Compartment")+
  scale_y_continuous(breaks = seq(from=5,to=35,by=5),labels = c("","10","","20","","30",""))+
  scale_x_continuous(breaks=c(10,30,50,70,90))+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),
        axis.line = element_line(colour="black"))
Figure2c

png("./Figures/Figure2c.png",res = 300,width = 5,height = 5,units = "in")
Figure2c
dev.off()

save(ineq_vacuole,file="./Models/ineq_vacuole")

###############Joint Figure 2A-C

####Plot a joint figure of 2A-C
Fig2a_legendless<-ggplot()+
  geom_col(aes(x=ineq_treatment,y=mean_root_QD,fill=variable),
           data=PlotData_Fig2A %>% filter(variable %in% c("root_QD_per_mg_root_rich_origin","root_QD_per_mg_root_poor_origin")),
           position = position_stack(reverse = F))+
  geom_errorbar(aes(x=ineq_treatment,ymin=mean_root_QD-se_root_QD,ymax=mean_root_QD+se_root_QD),
                data=PlotData_Fig2A %>% filter(variable=="root_QD_per_mg_root_both_origin"),width=0.25)+
  scale_fill_manual(values = c("#00BFC4","#F8766D"),
                    name="Nutrient\nCompartment",
                    breaks=c("root_QD_per_mg_root_rich_origin","root_QD_per_mg_root_poor_origin"),
                    labels=c("Rich","Poor"))+
  labs(x="Level of Inequality",y="P transfer to roots (nmol QD-apatite per mg root)")+
  #scale_y_continuous(limits = c(0,0.08),breaks = seq(from=0.01,to=0.07,by=0.01),labels = c("","0.02","","0.04","","0.06",""))+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),
        axis.line = element_line(colour="black"))+
  theme(legend.position="none")+  
  geom_text(aes(x=3,y=0.0012,label="*",size=30))
Fig2a


Fig2a_legendless

Figure2B_legendless<-ggplot()+
  geom_col(aes(x=ineq_treatment,y=mean_total_retention,fill=compartment),
           data=ineq_hyphal_retained_stats_melt_fig,position = position_stack(reverse = TRUE))+
  geom_errorbar(aes(x=ineq_treatment,ymin=mean_retention-se_retention,
                    ymax=mean_retention+se_retention),data=ineq_hyphal_retained_stats,
                width=0.25)+
  scale_fill_discrete(name="Nutrient\nCompartment")+
  #theme(axis.text.x=element_text(angle = -45, hjust = 0))+
  ylab("Total P retained by hyphae (nmol QD-apatite)")+
  xlab("Level of Inequality")+
  scale_y_continuous(breaks = seq(from=0.01,to=0.12,by=0.01),labels = c("","0.02","","0.04","","0.06","","0.08","","0.10","","0.12"))+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),
        axis.line = element_line(colour="black"))+
  guides(fill=guide_legend(reverse=TRUE),colour=guide_legend(reverse=TRUE))+
  theme(legend.position="none")+
  geom_text(aes(x=2,y=0.1,label="*",size=30))+
  geom_text(aes(x=3,y=0.03,label="***",size=30))
Figure2B_legendless

#Print the joint figure
pdf("./Figures/Figure2_AllPanels_Horizontal.pdf",width = 15,height = 5)
grid.arrange(Fig2a_legendless,Figure2B_legendless,Figure2c,nrow=1,widths=c(1,1,1))
dev.off()
#Using Illustrator, manually add the legend to the open white space in Figure 2B


# pdf("./Figures/ForLegendOnly.pdf")
# Figure2B
# dev.off()
 
# Figure 3 ---------------------------------------------------------------

#Analyse hyphal non-origin and origin P. How is it distributed? 

#Data formatting
ineq_P_origin<- ineq %>% filter(comp!="Root") %>%
  mutate(perc_other_origin=hyphal_QD_total_from_other_hypal/(hyphal_QD_total+hyphal_QD_total_from_other_hypal)) %>%
  select(plate,comp,ineq_treatment,perc_other_origin) %>% dcast(formula = plate+ineq_treatment~comp) %>%
  mutate(net_movement=Poor-Rich)
head(ineq_P_origin)

##Descriptives
#Look at the percentage of P in each compartment with originates from the other compartment (i.e. which has moved between fugnal compartemtns)
ggstripchart(ineq_P_origin,x = "ineq_treatment",y = c("Poor","Rich"),
             add="mean_se",color="ineq_treatment",combine=T,
             ylab="Percentage_non_origin")
#Look at the net movements
ggstripchart(ineq_P_origin,x = "ineq_treatment",y = c("net_movement"),
             add="mean_se",color="ineq_treatment",combine=T,
             ylab="Percentage_non_origin")
#Remove the ones where we can't calculate net movements. 
ineq_P_origin<-ineq_P_origin %>% filter(!is.na(net_movement))

##Statistical Modelling

ineq_P_origin_melted<-
  ineq_P_origin %>% 
  select(plate,ineq_treatment,Poor,Rich) %>%
  melt(id.vars=c("plate","ineq_treatment"),value.name = "P_transfer",variable.name="direction")
head(ineq_P_origin_melted)

#Evaluate a generalised linear mixed model. 
p_movement_glmer<-
  glmer(P_transfer~ineq_treatment*direction+(1|plate),
    data=ineq_P_origin_melted,family = gaussian(link="log"))

#Are the assumptions met? 
plot(p_movement_glmer)
qqnorm(resid(p_movement_glmer))
qqline(resid(p_movement_glmer))

#What does the model tell us? 
summary(p_movement_glmer)
Anova(p_movement_glmer)
save(p_movement_glmer,file="./Models/p_movement_glmer")
save(ineq_P_origin_melted,file="./Models/ineq_P_origin_melted")

##Figure

#Calculate Figure stats
ineq_P_origin_figure_stats_means<-ineq_P_origin %>%
  group_by(ineq_treatment) %>%
  summarise(mean_rich_to_poor=mean(Poor,na.rm=T),
            mean_poor_to_rich=mean(Rich,na.rm = T),
            mean_net_movement=mean(net_movement,na.rm=T)) %>%
            melt(id.vars="ineq_treatment",value.name = "mean_value") %>%
  mutate(movement_direction=case_when(
    variable == "mean_rich_to_poor" ~ "rich_to_poor",
    variable == "mean_poor_to_rich" ~ "poor_to_rich",
    variable == "mean_net_movement" ~ "net_movement"
  )) %>% select(-variable)
ineq_P_origin_figure_stats_ses<-ineq_P_origin %>%
  group_by(ineq_treatment) %>%
  summarise(se_rich_to_poor=stderr(Poor),
            se_poor_to_rich=stderr(Rich),
            se_net_movement=stderr(net_movement)) %>%
  melt(id.vars="ineq_treatment",value.name = "se_value") %>%
  mutate(movement_direction=case_when(
    variable == "se_rich_to_poor" ~ "rich_to_poor",
    variable == "se_poor_to_rich" ~ "poor_to_rich",
    variable == "se_net_movement" ~ "net_movement"
  )) %>% select(-variable)
ineq_P_origin_figure_stats<-inner_join(ineq_P_origin_figure_stats_means,ineq_P_origin_figure_stats_ses)
ineq_P_origin_figure_stats
ineq_P_origin_figure_stats$movement_direction<-factor(ineq_P_origin_figure_stats$movement_direction,
                                                      levels=c("rich_to_poor",
                                                               "poor_to_rich",
                                                               "net_movement"))

#Plot the figure  
Figure3<-ggplot()+
  geom_col(aes(x=ineq_treatment,y=mean_value,fill=movement_direction),
           data=ineq_P_origin_figure_stats,position=position_dodge(width = 0.9))+
  geom_errorbar(aes(x=ineq_treatment,ymin=mean_value-se_value,ymax=mean_value+se_value,fill=movement_direction),
                data=ineq_P_origin_figure_stats,position=position_dodge(width = 0.9),width=0.35)+
  scale_fill_brewer(name="Movement\nDirection",labels=c("Rich to Poor","Poor to Rich","Net Movement"))+
  ylab("% P originating from other nutrient compartment")+
  xlab("Levels of Inequality")+
  scale_y_continuous(breaks = seq(from=-5,to=65,by=5),
                     labels = c("","0","","10","","20","","30","","40","","50","","60",""))+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),
        axis.line = element_line(colour="black"))
Figure3

png("./Figures/Figure3.png",res = 300,width = 5,height = 5,units = "in")
Figure3
dev.off()

#Save the test results and data analysed
save(kruskal_non_origin_P_net_movement,file="./Models/kruskal_non_origin_P_net_movement")
save(ineq_P_origin,file="./Models/ineq_P_origin")

# SOM ---------------------------------------------------------------------

#######S4 Linear 

##Data formatting
head(ineq_nut_level_control)

#Descriptives
ggstripchart(ineq_nut_level_control,x = "comp_nutrients",
             y = c("root_QD_per_mg_root"),
             add="mean_se",ylab="P transfer to roots (nmol/mg)",combine=T)

##Model building
lm_QD_per_mg_nut_level_control<-lm(root_QD_per_mg_root~comp_nutrients,
                        data=ineq_nut_level_control)
residualPlot(lm_QD_per_mg_nut_level_control)
qqPlot(lm_QD_per_mg_nut_level_control)
#This isn't perfect, but looks reasonable. If the p-values are convincing, probably the departure from normality won't be too influential. 

#Save the model
save(lm_QD_per_mg_nut_level_control,file="./Models/lm_QD_per_mg_nut_level_control")

#What does it tell us?
summary(lm_QD_per_mg_nut_level_control)
Anova(lm_QD_per_mg_nut_level_control)

##Visualisation
ineq_nut_level_control_fig<-ineq_nut_level_control %>% group_by(comp_nutrients) %>% 
  summarise(mean_root_QD_per_mg_root=mean(root_QD_per_mg_root,na.rm = T),
            se_root=stderr(root_QD_per_mg_root))

Figure_S4Linear<-ggplot()+
  geom_pointrange(aes(x=comp_nutrients,ymin=mean_root_QD_per_mg_root-se_root,ymax=mean_root_QD_per_mg_root+se_root,
                      y=mean_root_QD_per_mg_root),
                  data=ineq_nut_level_control_fig,size=0.5)+
  geom_abline(intercept = 2.735853e-04,slope=3.447013e-06)+
  scale_y_continuous(limits = c(0,0.001))+
  xlab("Nutrient Compartment")+
  ylab("P transfer to roots (nmol QD-apatite per mg root)")+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),
        axis.line = element_line(colour="black"))+
  scale_x_continuous(breaks=c(10,30,50,70,90))
Figure_S4Linear

png("./Figures/Figure_S4_Linear.png",res = 300,width = 5,height = 5,units = "in")
Figure_S4Linear
dev.off()

save(ineq_nut_level_control,file="./Models/ineq_nut_level_control")

#######S5 Spore Growth
head(spore_growth)

#Descriptives
ggstripchart(spore_growth,x = "Treatment",
             y = c("Total.Growth"),
             add=c("mean_se","violin"),ylab="Spore Growth (cm)",combine=T)

#Print this figure, no formal statistical modelling is needed here. 
png("./Figures/Figure_S5Spore_Growth.png",res = 300,width = 5,height = 5,units = "in")
ggstripchart(spore_growth,x = "Treatment",
             y = c("Total.Growth"),
             add=c("mean_se","violin"),ylab="Hyphal Growth (cm)",combine=T)
dev.off()

#How many have grown at all, or more than 1 cm? 
spore_growth$grown<-ifelse(spore_growth$Total.Growth>0,1,0)
spore_growth$grown_more_1<-ifelse(spore_growth$Total.Growth>1,1,0)

table(spore_growth$Treatment,spore_growth$grown)
table(spore_growth$Treatment,spore_growth$grown_more_1)

#What is the longest growth? 
max(spore_growth$Total.Growth)

save(spore_growth,file="./Models/spore_growth")

########S6 In-Vitro Colour Control

head(in_vitro_dat)

#Diagnostics
in_vitro_dat_melt<-in_vitro_dat %>% select(plate,cyan_488nm_nmol_mg_dry_root,red_666nm_nmol_mg_dry_root) %>%
  melt(id.vars="plate",variable.name="colour_tested",value.name = "QD_per_mg_root")
in_vitro_dat_melt

#Descriptives
ggstripchart(in_vitro_dat_melt,x = "colour_tested",y="QD_per_mg_root",add="mean_se",
             font.label = list(size=2))

#Statistical modelling
in_vitro_t_test<-t.test(in_vitro_dat$cyan_488nm_nmol_mg_dry_root,
                        in_vitro_dat$red_666nm_nmol_mg_dry_root,paired = T)
in_vitro_t_test
save(in_vitro_t_test,file="./Models/in_vitro_t_test")

#Visualisation
in_vitro_dat_melt_fig<-in_vitro_dat_melt %>% group_by(colour_tested) %>%
  summarise(mean_QD_per_mg_root=mean(QD_per_mg_root),
            se_QD_per_mg_root=stderr(QD_per_mg_root))
in_vitro_dat_melt_fig$colour_tested<-c("488nm (cyan)","666nm (red)")
in_vitro_dat_melt_fig$colour_tested<-factor(in_vitro_dat_melt_fig$colour_tested,levels=c("666nm (red)","488nm (cyan)"))
head(in_vitro_dat_melt_fig)

FigS6_in_vitro<-ggplot()+
  geom_col(aes(x=colour_tested,y=mean_QD_per_mg_root,fill=colour_tested),
           data=in_vitro_dat_melt_fig,
           position = position_stack(reverse = F))+
  geom_errorbar(aes(x=colour_tested,ymin=mean_QD_per_mg_root-se_QD_per_mg_root,ymax=mean_QD_per_mg_root+se_QD_per_mg_root),
                data=in_vitro_dat_melt_fig,width=0.25)+
  labs(x="Colour",y="P transfer to roots (nmol QD-apatite per mg root)")+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),
        axis.line = element_line(colour="black"))
FigS6_in_vitro

png("./Figures/Figure_S6_in_vitro_colour.png",res = 300,width = 5,height = 5,units = "in")
FigS6_in_vitro
dev.off()

########S7 Root Biomass

#We reuse the plate-level data frame from above

#Descriptives
ggstripchart(ineq_root_comp,x="ineq_treatment",
             y=c("total_root_bm_mg"),
             add=c("mean_se","violin"),ylab = "Root biomass (mg)",combine=T)

##Statistical Modelling
lm_total_root_bm_mg_ineq<-lm(total_root_bm_mg~ineq_treatment,data=ineq_root_comp)

#Visually check if model assumptions are met
residualPlot(lm_total_root_bm_mg_ineq)
qqPlot(lm_total_root_bm_mg_ineq)

#What does the model tell us? 
summary(lm_total_root_bm_mg_ineq)
Anova(lm_total_root_bm_mg_ineq)

##Figure
FigS7_Root_BM<-
  ineq_root_comp %>% group_by(ineq_treatment) %>%
  summarise(mean_total_root_bm_mg=mean(total_root_bm_mg,na.rm = T),
            se_total_root_bm_mg=stderr(total_root_bm_mg)) %>%
  ggplot()+
  geom_col(aes(x=ineq_treatment,y=mean_total_root_bm_mg))+
  geom_errorbar(aes(x=ineq_treatment,ymin=mean_total_root_bm_mg-se_total_root_bm_mg,ymax=mean_total_root_bm_mg+se_total_root_bm_mg),width=0.25)+
  labs(x="Level of Inequality",y="Root biomass (mg)")+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),
        axis.line = element_line(colour="black"))
FigS7_Root_BM

png("./Figures/Figure_S7_RootBM.png",res = 300,width = 5,height = 5,units = "in")
FigS7_Root_BM
dev.off()

#Save the model
save(lm_total_root_bm_mg_ineq,file="./Models/lm_total_root_bm_mg_ineq")

##########S8 Colonisation

#Data formatting
colon$ineq_treatment<-factor(colon$ineq_treatment,levels = c("None","Medium","High"))

#Descriptives. 
ggstripchart(data = colon,x="ineq_treatment",y="perc_colonised",
             add = "mean_se",color = "ineq_treatment")

#Analyses
lm_colon<-lm(perc_colonised~ineq_treatment,data=colon)

#Visually check if model assumptions are met
residualPlot(lm_colon)
qqPlot(lm_colon)
#Not perfect, but reasonable enough given the small sample size. 

#Wht does the model tell us? 
summary(lm_colon)
Anova(lm_colon)

#Alternatively, with a glm.
glm_colon<-glm((perc_colonised/100)~ineq_treatment,family = "binomial",data=colon)
summary(glm_colon)
Anova(glm_colon)

#Figure
FigS8_colon<-
  colon %>% group_by(ineq_treatment) %>%
  summarise(mean_perc_colon=mean(perc_colonised),
            se_perc_colon=stderr(perc_colonised)) %>%
  ggplot()+
  geom_col(aes(x=ineq_treatment,y=mean_perc_colon))+
  geom_errorbar(aes(x=ineq_treatment,ymin=mean_perc_colon-se_perc_colon,
                    ymax=mean_perc_colon+se_perc_colon),width=0.25)+
  labs(x="Level of Inequality",y="Mycorrhizal Colonisation (%)")+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),
        axis.line = element_line(colour="black"))
FigS8_colon

png("./Figures/Figure_S8_root_colonisation.png",res = 300,width = 5,height = 5,units = "in")
FigS8_colon
dev.off()

#Save the model and data frame
save(lm_colon,file="./Models/lm_colon")
save(colon,file="./Models/colon")

########### S9 Hyphal Biomass Overall

#We are here evaluating the full hyphal biomass per plate, across all compartments. 
#Calculate the plate level fungal biomass
ineq_hypal_bm<-ineq %>% group_by(plate) %>%
  summarise(plate_hyphal_bm=sum(hyphal_bm_mg),
            ineq_treatment=ineq_treatment[comp=="Poor"])
ineq_hypal_bm

#Descriptives. 
ggstripchart(data = ineq_hypal_bm,x="ineq_treatment",y="plate_hyphal_bm",
             add = "mean_se",color = "ineq_treatment")

#Modelling
lm_plate_bm<-lm(plate_hyphal_bm~ineq_treatment,data=ineq_hypal_bm)
save(lm_plate_bm,file="./Models/lm_plate_bm")

#Visually check if model assumptions are met
residualPlot(lm_plate_bm)
qqPlot(lm_plate_bm)
#Not good, evaluate a GLM

#GLM
glm_plate_bm_gamma_identity<-glm(plate_hyphal_bm~ineq_treatment,data=ineq_hypal_bm,family = Gamma(link = "identity"))
save(glm_plate_bm_gamma_identity,file="./Models/glm_plate_bm_gamma_identity")
#Visually evaluate the model assumptions
glm.diag.plots(glm_plate_bm_gamma_identity)

#What does the model tell us? 
summary(glm_plate_bm_gamma_identity)
Anova(glm_plate_bm_gamma_identity,test.statistic ="F")

#Visualise it
ineq_hypal_bm_stats <- ineq_hypal_bm %>% group_by(ineq_treatment) %>%
  summarise(mean_biomass=mean(plate_hyphal_bm,na.rm=T),
            se_biomass=stderr(plate_hyphal_bm))
ineq_hypal_bm_stats

FigureS9_hyphal_bm_overall<-ggplot()+
  geom_col(aes(x=ineq_treatment,y=mean_biomass),
           data=ineq_hypal_bm_stats)+
  geom_errorbar(aes(x=ineq_treatment,ymin=mean_biomass-se_biomass,ymax=mean_biomass+se_biomass),
                data=ineq_hypal_bm_stats,
                width=0.25)+
  ylab("Total fungal biomass (mg)")+
  xlab("Levels of Inequality")+
  theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_x_discrete(labels=c("None","Medium","High"))
FigureS9_hyphal_bm_overall

png("./Figures/Figure_S9_OverallHyphalBM.png",res = 300,width = 5,height = 5,units = "in")
FigureS9_hyphal_bm_overall
dev.off()

save(ineq_hypal_bm,file="./Models/ineq_hyphal_bm")

########S15 Whole Plant Colour Test

#Data formatting
head(whole_plant_dat)
whole_plant_dat$QD_total_shoot<-whole_plant_dat$QD_mg_shoot*whole_plant_dat$total_shoot_mass_g*1000
whole_plant_dat$QD_total_root<-whole_plant_dat$QD_mg_root*whole_plant_dat$total_root_mass_g*1000
whole_plant_dat$QD_total_plant<-whole_plant_dat$QD_total_shoot+whole_plant_dat$QD_total_root
whole_plant_dat$QD_treatment<-factor(whole_plant_dat$QD_treatment,c("666nm","488nm"))
#For this analysis, we only consider the equal treatment and mycorrhizal plants 
whole_plant_equal_dat<-whole_plant_dat %>% filter(ineq_treatment=="None") 
head(whole_plant_equal_dat)

#Descriptives

#For total QDs in the whole plant
ggstripchart(whole_plant_equal_dat,x = "QD_treatment",y=c("QD_total_plant"),
             add=c("mean_se"),color="QD_treatment",combine=T)

#Statistical Modelling

#Analyse shoots
QD_total_plant_488<-whole_plant_equal_dat$QD_total_plant[whole_plant_equal_dat$QD_treatment=="488nm"]
QD_total_plant_666<-whole_plant_equal_dat$QD_total_plant[whole_plant_equal_dat$QD_treatment=="666nm"]
QD_total_plant_488
QD_total_plant_666

#Paired t.test
whole_plant_5050_colour_t_test<-t.test(QD_total_plant_488,QD_total_plant_666,paired = T)
whole_plant_5050_colour_t_test
save(whole_plant_5050_colour_t_test,file="./Models/whole_plant_5050_colour_t_test")

#Visualisation
FigS15_Colours<-whole_plant_equal_dat %>%  
  group_by(QD_treatment) %>%
  summarise(mean_QDs=mean(QD_total_plant),
            se_QDs=stderr(QD_total_plant)) %>%
  ggplot()+
  geom_col(aes(x=QD_treatment,y=mean_QDs,fill=QD_treatment))+
  geom_errorbar(aes(x=QD_treatment,ymin=mean_QDs-se_QDs,ymax=mean_QDs+se_QDs), 
                width=0.25)+
  ylab("Total P transfer to plant (nmol QD-apatite per plant)")+
  xlab("Inoculation")+
  theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.position = "none")+
  scale_x_discrete(labels=c("666nm (Red)","488nm (Cyan)"))
FigS15_Colours

png("./Figures/FigS15_WholePlantEqualColours.png",res = 300,width = 5,height = 5,units = "in")
FigS15_Colours
dev.off()


#########S16 Colours switch

##Data Formatting
#Subset to only the unequal treatment and the mycorrhizal plants. 
whole_plant_colour_switch_uneq_M<-whole_plant_dat %>% filter(ineq_treatment=="High" & mycorrhizal_treatment=="M")
whole_plant_colour_switch_uneq_M$comp<-as.factor(whole_plant_colour_switch_uneq_M$comp)
head(whole_plant_colour_switch_uneq_M)

##Descriptives
#Have a look at the data
ggstripchart(whole_plant_colour_switch_uneq_M,x = "comp",y=c("QD_total_plant"),
             add=c("mean_se"),color="QD_treatment")

##Statistical Modelling
#Evaluate a linear mixed model
lme_col_switch<-lme(QD_total_plant~comp+QD_treatment,random = ~1|pot,data = whole_plant_colour_switch_uneq_M)

#Visually check if assumptions of the model are met
plot(lme_col_switch)
qqnorm(lme_col_switch$residuals)
qqline(lme_col_switch$residuals)
#Assumptions not met, try a Generalised linear mixed model

#GLMMS
glmer_col_switch<-glmer(QD_total_plant~comp+QD_treatment+(1|pot),
                        data = whole_plant_colour_switch_uneq_M,family=inverse.gaussian(link="identity"))

save(glmer_col_switch,file="./Models/glmer_col_switch")
plot(glmer_col_switch)
qqnorm(resid(glmer_col_switch))
qqline(resid(glmer_col_switch))
#This looks quite ok, not perfect. 

#What does the model tell us? 
summary(glmer_col_switch)
Anova(glmer_col_switch)

##Figure
whole_plant_colour_switch_uneq_M_figstats<-whole_plant_colour_switch_uneq_M %>% group_by(comp,QD_treatment) %>%
  summarise(mean_QD_total_whole_plant=mean(QD_total_plant),
            se_QD_total_whole_plant=stderr(QD_total_plant))
whole_plant_colour_switch_uneq_M_figstats$QD_treatment<-factor(whole_plant_colour_switch_uneq_M_figstats$QD_treatment,c("666nm","488nm"))

FigS16_ColourSwitch<-
  whole_plant_colour_switch_uneq_M_figstats %>% ggplot()+
  geom_col(aes(x=comp,y=mean_QD_total_whole_plant,fill=QD_treatment),position="dodge")+
  geom_errorbar(aes(x=comp,
                    ymin=mean_QD_total_whole_plant-se_QD_total_whole_plant,
                    ymax=mean_QD_total_whole_plant+se_QD_total_whole_plant,group=QD_treatment),
                position=position_dodge(width=0.9),
                width=0.25)+
  ylab("Total QD uptake in whole plant (nmol)")+
  xlab("Source Compartment")+
  scale_fill_discrete(name="Colour\nQD-label")+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),
        axis.line = element_line(colour="black"))
FigS16_ColourSwitch

png("./Figures/FigS16_WholePlantColourSwitch.png",res = 300,width = 5,height = 5,units = "in")
FigS16_ColourSwitch
dev.off()

save(whole_plant_colour_switch_uneq_M,file="./Models/whole_plant_colour_switch_uneq_M")

###### S17 Toxicity Data 

head(tox_dat)

## Data Formatting
#calculate full plant mass
tox_dat$full_plant_mass_g<-tox_dat$shoot_g+tox_dat$root_g

##Descriptives

#Have a look
ggstripchart(data = tox_dat,x="p_source",y="full_plant_mass_g",
             add = "mean_se")

## Statistical Modelling

#Evaluate a linear model
lm_tox<-lm(full_plant_mass_g~p_source,data=tox_dat)
save(lm_tox,file="./Models/lm_tox")

##Visually check if assumptions of the model are met
residualPlot(lm_tox)
qqPlot(lm_tox)
#Assumptions seem met

#What does the model tell us? 
summary(lm_tox)
Anova(lm_tox)

## Figure
tox_dat_figdata <- tox_dat %>% 
  group_by(p_source) %>%
  summarise(mean_full_plant_mass_g=mean(full_plant_mass_g),
            se_full_plant_mass_g=stderr(full_plant_mass_g))
tox_dat_figdata

#Plot the figure
FigureS17_tox<-ggplot()+
  geom_col(aes(x=p_source,y=mean_full_plant_mass_g),
           data=tox_dat_figdata)+
  geom_errorbar(aes(x=p_source,ymin=mean_full_plant_mass_g-se_full_plant_mass_g,ymax=mean_full_plant_mass_g+se_full_plant_mass_g),
                data=tox_dat_figdata,
                width=0.25)+
  ylab("Total plant biomass (g)")+
  xlab("Phosphorus source")+
  theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_x_discrete(labels=c("Apatite","QD-Apatite"))
FigureS17_tox

png("./Figures/FigS17_WholePlantToxicity.png",res = 300,width = 5,height = 5,units = "in")
FigureS17_tox
dev.off()

save(tox_dat,file="./Models/tox_dat")


