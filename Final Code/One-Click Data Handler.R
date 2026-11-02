#changes  both here ad the amrkdown:
#   renamed first image export name to od_spread_...
#   outliers_xxxx... now belongs to second image set pdf
#   new section of merging and filtering and further stuff
#   t-test, data export, such and such

library(tidyverse)
library(readxl)
library(ggplate)
library(diptest)
# Downloading the custom package
# Line can be removed after the first time running it
remotes::install_github("BeefMaster64/platereaderr")
library(platereaderr)

#importing the data and creating save location ================

#select the data to be worked with, this has the name of 'xxxxxx.txt' where x's are the date of measure
message("Please select the raw input file!")
working_data <- file.choose()
#the message + file import could create its own function for the one-click code

#choose save location folder for exporting the final data and images
message("Please select the file output folder for saving data!\nA separate fodler wil be created inside the selected folder")
output_dir <- choose.dir(caption = "Select a save folder")

#checking to generate images or not
answer <- menu(choices = c("Yes", "No"),
               title = "Do you wish to generate quality checking images?")

generate_images <- answer == 1

#extracting date of measurement
#basename removes all parts of a directory but the last (aka. the file name)
measurement_date <- working_data |> 
  basename() |> 
  str_remove(pattern = "\\.txt$")

#creating a new folder in the save location based on measurement date
save_folder <- file.path(output_dir, measurement_date)

dir.create(save_folder, showWarnings = FALSE)


#parsing the data from the raw output file
od_tables <- od_plate_reader(working_data)

#splitting the data into 4 smaller plates
#optional flattening, but I don't like it all that much, might make further work harder 
#still kept it, might help simplify code later
od_tables_96 <- od_tables |> 
  map(split_384_plate) |> 
  list_flatten(name_spec = "{outer}_{inner}")


#assigning the additives to the correct measurements =================
#turning the od_tables_96 into a long format
od_measurements_long <- imap_dfr(od_tables_96, function(plate, table_name) {
  
  plate |> 
    rownames_to_column(var = "Row") |> 
    pivot_longer(cols = -Row,
                 names_to = "Column",
                 values_to = "OD") |>  #makes it the same layout as the reference above
    mutate(Column = as.integer(Column)) |> 
    
    #Defining the 'plate position'
    mutate(Plate_num = as.integer(str_extract(table_name, "[1-4]$")),
           
           #naming the strain used in the measurement
           Strain = str_remove(table_name, "_plate[1-4]$"))
  
})

#pairing up the measurements with the correct additive based on well position
od_additives <- od_measurements_long |> 
  left_join(y = layout_reference,
            by = c("Row", "Column"))

#pivoting to resemble the final data frame
#Keeping the well positio markers for now
#in case it might be needed (I don't think so)
od_with_well_pos <- od_additives |> 
  pivot_wider(names_from = Strain,
              values_from = OD) |> 
  arrange(Plate_num, Row, Column)

#creating the well type descriptors
#'empty' wells on plate 1,2,3 are empty 'strain' wells'
#'empty' wells on palte 4 are 'additive_background'
#'-' wells are 'frame'
#'the rest, aka. additive well on plate 1,2,3 are 'additive'
#'getting the order right was tricky :p
od_table_wide <- od_with_well_pos |> 
  mutate(well_type = case_when(Additive == "-" ~ "frame",
                               Plate_num == 4 ~ "additive_background",
                               Additive == "empty" & Plate_num %in% 1:3 ~ "strain_control",
                               !is.na(Additive) & Plate_num %in% 1:3 ~ "additive",
                               .default = NA_character_), .after = Additive) |> 
  unite(col = "pos_96", Row, Column, sep = "")


#saving this middle od table in case it's needed
write_csv(x = od_table_wide, file = file.path(save_folder, paste0("plate_", measurement_date, ".csv")))


#merge with the script_v2.R code ================ 
od_table_long <- od_table_wide |> 
  pivot_longer(cols = 5:ncol(od_table_wide), names_to = "strain", values_to = "OD") |> 
  rename(additive = Additive,
         plate = Plate_num)

#visualization 
#only generates if enabled (generate_images = TRUE)
if (generate_images) {
  
  #plotting the control plate composition
  strains <- unique(od_table_long$strain)
  plate_images <- list()
  
  for (i in strains) {
    
    plate_images[[i]] <- od_table_long %>%
      filter(strain == i, plate == 4) %>%
      plate_plot(
        position = pos_96,
        value = OD,
        plate_size = 96,
        plate_type = "round",
        title = paste(i, "-position4"))  
  }
  
  #the theme was predetermined, no idea if it does anything
  plate_box <- ggplot(data = od_table_long |> 
                        filter(well_type == "strain_control"),
                      aes(x = strain, y = OD)) +
    geom_boxplot() +
    geom_point(position = position_jitter(width = 0.1)) +
    scale_y_log10() +
    theme(axis.text.x = element_text(angle = 40, vjust = 1, hjust = 1)) +
    labs(x = "Strains", y = "OD value",
         title = "OD value spread of untreated strain controls") +
    theme(axis.text.x = element_text(angle = 40, vjust = 1, hjust = 1),
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
  
  pdf(file = file.path(save_folder, paste0("plate_check_", measurement_date, ".pdf")), width = 7, height = 6)
  for (plot in plate_images) {
    print(plot)
  }
  plate_box
  dev.off()
  
}

#background subtraction and relative OD calculation ===========
additive_background <- od_table_long |> 
  filter(well_type == "additive_background") |> 
  group_by(additive) |> 
  summarise(background_od = mean(OD))


background_data <- od_table_long |> 
  filter(well_type != "frame") |> 
  left_join(additive_background, by = "additive") |> 
  mutate(od_wo_background = OD - background_od)

untreated_strain <- background_data |> 
  filter(well_type == "strain_control") |> 
  group_by(strain) |> 
  summarise(untreated_od = mean(od_wo_background))

#calculated data contains ALL columns that one would ever need for double checking
calculated_data <-  background_data |> 
  mutate(identifier = paste0("id_", measurement_date, " ", strain), .before = additive) |> 
  left_join(untreated_strain, by = "strain") |>
  mutate(relative_od = od_wo_background/untreated_od)

#clean data contains only the columns used for further analysis
clean_data <- calculated_data |> 
  filter(well_type != "additive_background") |> 
  select(-pos_96, -plate, -well_type, -OD, -background_od)

#exporting the clean results
write_csv(x = clean_data, file = file.path(save_folder, paste0("clean_data_", measurement_date, ".csv")))

#outlier filtering steps ===============
# applying the full outlier filtering makes no sense here, since the
# hartigan and MAD filter need the data from multiple measurement days to function

#Hampel filter function
Hampel <- function(value, mad_mult = 3, bound = "both") {
  
  # bound <- match.arg(bound)
  
  median <- median(value)
  mad <- mad(value, constant = 1)
  
  lower_bound <- median - mad_mult*mad
  upper_bound <- median + mad_mult*mad
  
  if (bound == "upper") {
    which(value > upper_bound)
  } else if (bound == "lower") {
    which(value < lower_bound)
  } else {
    which(value < lower_bound | value > upper_bound)
  }
  
  #might need an update to return a Boolean vector [just remove "which"]
  
}

outlier_plotting <- clean_data |> 
  select(!untreated_od & !relative_od) |> 
  mutate(med = median(od_wo_background), .by = c(strain ,identifier, additive)) |>  
  mutate(norm_od = od_wo_background - med) |> 
  mutate(outlier = case_when(additive == "empty" ~ 2,
                             row_number() %in% Hampel(norm_od, 5) ~ 1, 
                             .default = 0)) 

if (generate_images) {
  plots <- list()
  for (s in strains) {
    plots[[s]] <-  ggplot(data = outlier_plotting |> 
                            filter(strain == s),
                          aes(x = additive, y = norm_od))+
      geom_point(aes(colour = factor(outlier, levels = c(0,1,2), labels = c("Non-outlier", "Hampel", "Control"))),
                 position = position_jitter(width = 0.2)) +
      labs(title = paste("Outliers of ", s),
           x = "Additive",
           y = "OD",
           color = "Outliers") +
      theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
      scale_color_manual(values = c("Non-outlier" = "#00F545", "Hampel" = "#FF6161", "Control" = "#AE6BFF"))
  }
  
  pdf(file = file.path(save_folder, paste0("od_spread_", measurement_date, ".pdf")), width = 7, height = 6)
  for (plot in plots) {
    print(plot)
  }
  dev.off()
}
