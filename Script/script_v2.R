setwd("C:/Users/lejla/Dropbox/project_M/screen/raw_data/260723")
install.packages("ggplate")
library(tidyverse)
library(ggplate)

#setting the theme ============
theme <-theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1 ,size=15, vjust=0.5),
        axis.text.y = element_text(angle = 0, hjust = 1,size=15, vjust=0.5),
        axis.title.x = element_text(size=20,
                                    margin = margin(t = 15, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size=20,
                                    margin = margin(t = 15, r = 0, b = 0, l = 0)),
        plot.title = element_text(size = 20, face = "italic"),
        strip.text.x = element_text(size = 10, face = "bold.italic"),
        strip.text.y = element_text(size = 10),
        strip.background = element_rect(fill = "azure3"),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 10))

theme_set(theme)


data <- read.csv("./data_261007.csv")

#second static number (9) is the number of columns, change it to dynamic number
data <- data %>% 
  pivot_longer(6:9, names_to = "strain", values_to = "OD")

#visualization and quality checking ====================
strains <- unique(data$strain)

#printing a pdf
pdf("data_260430.pdf", width = 7, height = 6)

for (s in strains) {
  p <- data %>%
    filter(strain == s, position == 4) %>%
    plate_plot(
      position = well_96,
      value = OD,
      plate_size = 96,
      plate_type = "round",
      title = paste(s, "-position4")
    )
  
  print(p)  
}


p1 <- data %>% 
  filter(well_type == "strain") %>% 
  ggplot(aes(strain, OD)) +
  geom_boxplot() +
  geom_jitter() +
  scale_y_log10()

print(p1)


dev.off()

#selecting one sus strain
 #data %>% 
#  filter(strain == "B_dorie") %>% 
#  plate_plot(
#    position = well_384,
#    value = OD,
#    plate_size = 384
#  )


##########################################################################################
data <- data %>% 
  filter(!is.na(OD))


background_OD_additive <- data %>% 
  filter(well_type == "additive_background") %>% 
 # mutate(additive = ifelse(additive == "BACI", 0, additive)) %>%  what do you mean??!!
  group_by(additive) %>% 
  summarise( background_additive = mean(OD)) %>% 
  ungroup()

data <- data %>% 
  filter(well_type != "frame")

data <- left_join(data, background_OD_additive, by= "additive")

data <- data %>% 
  mutate(OD_wo_background = OD - background_additive)

#this is also done in my own script, so choose which one to use 
background_OD_strain <- data %>% 
  filter(position != "4", additive == "empty") %>% 
  group_by(strain) %>% 
  summarise(untreated_strain = mean(OD_wo_background)) %>% 
  ungroup()


data <- left_join(data, background_OD_strain, by= "strain")


#this datafram is falmiliar,as this is the one I've been receiving
data <- data %>% 
  mutate(relative_OD = OD_wo_background/untreated_strain)


#writing here will be redundant in the finished combined script
write_csv(data, "calculated_data_260723.csv")



