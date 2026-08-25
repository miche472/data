library(tidyverse)
library(lubridate)
library(ggpubr)
set_theme(theme_pubr())
set_palette("npg")
#set working directory to where I have this file...currently in my downloads folder
raw_dataframe <- read.csv("../Downloads/Zimmermann.csv")
zim<-raw_dataframe %>% filter(experiment_id=="Zimmerman_2020")

#Females during the light period
zim %>% filter(time_of_day=="light" & sex=="female") %>% 
  ggplot(.,aes(x=body_mass_g,y=energy_expenditure_kcal_hr,color=strain,fill=strain))+geom_point(shape=21, color="black",alpha=0.5)+geom_smooth(method="lm",se=F)+facet_wrap(~sex)


#Light period BW vs EE
zim %>% filter(time_of_day=="light") %>% 
  ggplot(.,aes(x=body_mass_g,y=energy_expenditure_kcal_hr,color=strain,fill=strain))+geom_point(shape=21, color="black",alpha=0.5)+geom_smooth(method="lm",se=F,alpha=0.5)+
  labs(x="Body mass (g)",y="EE (kcal/hr)")

#
zim %>% filter(time_of_day=="light") %>% 
  ggplot(.,aes(x=body_mass_g,y=energy_expenditure_kcal_hr))+geom_point(shape=21, color="black",alpha=0.5,aes(fill=strain))+geom_smooth(method="loess",span=1,se=F,alpha=0.5,color="black")+
  labs(x="Body mass (g)",y="EE (kcal/hr)")
ggsave("output plots/Zimmerman_light_3.png",width=10, height=8, units="in",dpi=600)

zim %>% filter(time_of_day!="Full day") %>% 
  ggplot(.,aes(x=body_mass_g,y=energy_expenditure_kcal_hr))+geom_point(shape=21, color="black",alpha=0.5,aes(fill=strain))+geom_smooth(method="loess",span=1,se=F,alpha=0.5,color="black")+
  labs(x="Body mass (g)",y="EE (kcal/hr)")+facet_wrap(~time_of_day)
ggsave("output plots/Zimmerman_both.png",width=8, height=8, units="in",dpi=600)

zim %>% filter(time_of_day=="light") %>% 
  ggplot(.,aes(x=body_mass_g,y=energy_expenditure_kcal_hr))+geom_point(shape=21, color="black",alpha=0.5,aes(fill=strain))+
  geom_smooth(method="loess",span=1,se=F,alpha=0.5,color="black")+
  labs(x="Body mass (g)",y="EE (kcal/hr)")+facet_wrap(~sex)

zim |> count(strain)

#nzo and b6
nzo_df<-zim |> filter(strain %in% c("C57BL/6J","NZO/HlLtJ"))
#nzo_df<-zim |> filter(strain %in% c("NZO/HlLtJ"))

plot1<-nzo_df %>%
filter(time_of_day=="Full day") %>% 
  ggplot(.,aes(x=body_mass_g,y=energy_expenditure_kcal_hr,group=strain))+geom_point(shape=21, color="black",alpha=0.5,aes(fill=strain))+geom_smooth(method="lm",span=1,se=F,alpha=0.5,color="black")+
  labs(x="Body mass (g)",y="EE (kcal/hr)")+facet_wrap(~sex)

set_palette(plot1, "jco")
ggsave("output plots/nzo_1.png",width=8, height=8, units="in",dpi=600)

plot2<-nzo_df %>%
  filter(time_of_day=="Full day") %>% 
  ggplot(.,aes(x=body_mass_g,y=energy_expenditure_kcal_hr,group=strain))+geom_point( alpha=0.5,aes(fill=strain,shape=sex))+
  geom_smooth(method="lm",span=1,se=F,alpha=0.5,color="black")+
  labs(x="Body mass (g)",y="EE (kcal/hr)")
set_palette(plot2, "jco")
ggsave("output plots/nzo_2.png",width=8, height=8, units="in",dpi=600)
nzo_stats<-nzo_df |> filter(time_of_day=="Full day")



stat_result<-glm(energy_expenditure_kcal_hr~body_mass_g+strain,data=nzo_df)
summary(stat_result)
