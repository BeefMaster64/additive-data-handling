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
#write_csv(x = od_table_wide, file = file.path(save_folder, paste0("plate_", measurement_date, ".csv")))


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
  print(plate_box)
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

#combining the results with all the previous data ============
# importing the previous data title 'merged_data_xxxxxx_yyyyyy'
# where x's are the oldest data date in the set, an y's are the latest data date
message("Please select the merged data of All previous measurements\n
The title should be \'merged_data_xxxxxx-yyyyyy.csv\'\nWhere X and Y are the measurement dates")
previous_data <- file.choose()

#merging the old and new data, solving naming convention issues
merged_data <- read_csv(previous_data) |> 
  bind_rows(clean_data |> 
              mutate(strain = sub(pattern = "_[0-9]+$", "", strain))) |> 
  mutate(identifier = identifier |>
           str_replace("^(id_[0-9]{6}) (.+)_1$", "\\1_A \\2") |>
           str_replace("^(id_[0-9]{6}) (.+)_2$", "\\1_B \\2"), #god the regex is getting complex...
         identifier = sub(pattern = "^([0-9]{6})", "id_\\1", identifier))

#exporting the new merged table to be used for future calculations
#the final name of the newly exported, updated data
export_name <- previous_data |> 
  basename() |> 
  str_remove(pattern = "(?<=-).*$") |> 
  paste0(measurement_date)

write_csv(x = merged_data, 
          file = file.path(output_dir, paste0(export_name, ".csv")))

#combined outlier elimination ---------
#this is pretty dirty code written by me ~a year ago, it works, but needs a revision

#hartigan's dip test of multimodality to remove entire problematic plates
hartigans_test <- merged_data |>
  select(!relative_od & !untreated_od) |>
  filter(additive == "empty") |> 
  summarise(hartigans_dip_p = dip.test(od_wo_background)$p.value, .by = identifier) |> 
  mutate(modality = if_else(hartigans_dip_p <= 0.05, 1, 0)) # 1 = multimodal outlier

# calculating stats for the next filter and joining previous filter results
plate_stats <- merged_data |> 
  filter(additive == "empty") |> 
  summarise(
    med = median(od_wo_background, na.rm = T),
    med_dev = mad(od_wo_background, na.rm = T),
    .by = c(strain, identifier)
  ) |> 
  left_join(hartigans_test, by = "identifier")

# MAD filter for plates that show an excessive spread of measured values
MAD_outliers <- plate_stats |> 
  mutate(madness = if_else(modality == 0,
                           true = if_else(row_number() %in% Hampel(med_dev, 3, "upper"), true = 2, false = 0),
                           false = 1), .by = strain) #2 = MAD filter outlier

#joining the outliers and normalizing the od for the next step
median_frame <- merged_data |> 
  filter(additive == "empty") |> 
  left_join(MAD_outliers, by = c("strain", "identifier")) |> 
  mutate(norm_od = od_wo_background - med) |> 
  select(strain, identifier, additive, od_wo_background, med, med_dev, norm_od, modality, madness)


outliers_combined <- median_frame |> 
  mutate(outlier = if_else(madness == 0,
                           true = if_else(row_number() %in% Hampel(norm_od, 6), 3, 0), # 3 = Hampel Outlier, 0 = non-outlier
                           false =  if_else(madness == 1, true = 1, false = 2)), .by = identifier)

#visualization of the outlier filtering
if (generate_images) {
  strain_list <- unique(merged_data$strain)
  plots <- list()
  
  for (i in strain_list) {
    plots[[i]] <- ggplot(data = outliers_combined |> 
                           filter(strain == i) |> 
                           mutate(identifier = sub(" .*", "", identifier)), 
                         aes(x = identifier, y = od_wo_background)) +
      geom_point(aes(color = factor(outlier, levels = c(0,1,2,3), labels = c("Non-outlier", "Multimodal", "MAD outlier", "Hampel outlier"))),
                 position = position_jitter(width = 0.2, height = 0)) +
      labs(title = paste("Outliers of", i, sep = " "),
           x = "Identifier",
           y = "OD",
           color = "Outliers") +
      theme(axis.text.x = element_text(angle = 40, hjust = 1)) +
      scale_color_manual(values = c("Non-outlier" = "#00F545", "Multimodal" = "#6B72FF", "MAD outlier" = "#AE6BFF", "Hampel outlier" = "#F86BFF"))
  }
  
  #exporting
  pdf(file = file.path(save_folder, paste0("outliers_", measurement_date, ".pdf")), width = 10)
  for (plot in plots) {
    print(plot)
  }
  dev.off()
}

#removing the outlier plates ==================
#such a dumb way to do it....
#but it works..... for now.... old code, needs revision

plates_to_remove <- outliers_combined |> 
  filter(outlier %in% c(1,2)) |> 
  distinct(identifier) |> 
  pull(identifier)

outlier_filtered_frame <- merged_data |> 
  filter(!identifier %in% plates_to_remove)


#filtering out outlier measurements
#applying the Hampel filter to treatments as well
non_outlier_plates <- outliers_combined |> 
  filter(outlier == 0) |> 
  select(!med_dev & !modality & !madness) |> 
  mutate(outlier = 2) |>  #2 = control values in current scope
  relocate(strain, .after = identifier)


additive_means_norms <- outlier_filtered_frame |> 
  select(!untreated_od & !relative_od) |> 
  mutate(med = median(od_wo_background), .by = c(strain ,identifier, additive)) |>  
  mutate(norm_od = od_wo_background - med)

#hampel filter on treated plates
#a bit strange of a solution, but seems to be working
additives_filtered_outliers <- additive_means_norms |> 
  filter(additive != "empty") |>
  mutate(outlier = if_else(row_number() %in% Hampel(norm_od, 5), true = 1, false = 0), .by = c(strain, additive)) |> 
  bind_rows(non_outlier_plates) #these are the controls


#performing the testing =========
#also old code but no way I'm revising this.... unless really needed...

#prparing the data for the test
t_test_data <- additives_filtered_outliers |> 
  group_by(identifier) |> 
  mutate(untreated_od =  #this section filters out the replicates with no control = all controls were outliers
           if (any(additive == "empty")) {
             mean(od_wo_background[additive == "empty"], na.rm = TRUE)
           } else {NA_real_}) |> 
  ungroup() |>  
  select(!med & !norm_od) |> 
  filter(additive != "empty",
         outlier != 1) |> #only here do I remove the Hampel-marked treated values
  mutate(mean_untreated_od = mean(untreated_od),
         mean_treated_od = mean(od_wo_background),
         median_untreated_od = median(untreated_od),
         median_treated_od = median(od_wo_background),
         .by = c("identifier", "additive")) |> 
  distinct(identifier, strain, additive, mean_untreated_od, mean_treated_od, median_untreated_od, median_treated_od)

#prep for the replicate reduction
#extracting the date and also handling it as a date (might be useful later idk?)
dated_data <- t_test_data |> 
  filter(complete.cases(median_untreated_od, median_treated_od)) |> 
  mutate(date_str = str_extract(string = identifier, pattern = "\\d{6}"),
         date = as.Date(date_str, format = '%y%m%d'), .before = identifier)


#creating the working dataframe
date_sorted_data <- dated_data |> 
  group_by(strain, additive) |> 
  
  #ordering newest first (remove 'desc()' to see with oldest-first)
  #identifier is used to handle same-date experiment ties in ordering
  arrange(desc(date), identifier) |>
  
  #counting replicates by group and ordering them based on position in the order
  mutate(replicates = n(),
         repli_order = row_number()) |> 
  
  #selecting only the pairs that have enough replicates to run all tests
  #this is to make sure all tests use the same data
  #replace 3 witha variable to set replicate upper limit, 3 is the standard currently
  filter(replicates >= 3,
         repli_order <= 3) |> 
  ungroup()

#paired t-test, can use both means and medians, but the difference is meager 
#using MEDIANs as of now due to not-insignificant spread in the data
t_test_results <- date_sorted_data |> 
  filter(complete.cases(median_untreated_od, median_treated_od)) |> 
  summarise(paired_p_value =
              if (n() >= 2 && sd(median_untreated_od - median_treated_od) > 0)
              {
                t.test(median_untreated_od, median_treated_od, paired = T)$p.value
              } else
              {
                NA_real_
              }, .by = c("strain", "additive"))

#fdr correction
t_test_results_fdr <- t_test_results |> 
  mutate(paired_FDR_p_value = p.adjust(paired_p_value, method = "fdr"), .by = strain)


#exporting the data  ------------
#preparing the data for exporting
#replicate count
replicate_count <- t_test_data |> 
  summarise(replicate_nr = n(), .by = c("strain", "additive"))

#tidying
final_test_results <- t_test_results_fdr |> 
  filter(additive != "empty") |>
  arrange(strain, additive) |> 
  mutate(p_value_fdr = paired_FDR_p_value) |> 
  select(!paired_FDR_p_value) |> 
  left_join(replicate_count, by = c("strain", "additive"))

#counting significant result 
#only counting strains which have at least 2 replicates
test_result_significance <- final_test_results |> 
  filter(replicate_nr >= 2) |> 
  summarise(total = n(),
            sign = sum(p_value_fdr < 0.1, na.rm = T), #default cutoff 0.05, use 0.1 instead from now
            sign_percent = paste(round(sign/total*100, digits = 3), "%", sep = ''),
            na_count = sum(is.na(p_value_fdr)))

#export
write_csv(x = final_test_results, file = file.path(save_folder, paste0("test_results_", measurement_date, ".csv")))
write_csv(x = test_result_significance, file = file.path(save_folder, paste0("test_significance_", measurement_date, ".csv")))

#exporting with relative OD
OD_values <- outlier_filtered_frame |> 
  select(identifier, strain, additive, relative_od) |> 
  mutate(relative_od = if_else(relative_od < 0, true = 0, false = relative_od)) |> 
  summarise(relative_od_mean  = mean(relative_od), .by = c("strain", "additive"))

results_for_plotting <- t_test_results_fdr |> 
  left_join(OD_values, by = c("strain", "additive"))

write_csv(x = results_for_plotting, file = file.path(save_folder, paste0("realtive_od_results_", measurement_date, ".csv")))

