###

setwd("~/Dropbox/trial_illuminaprocessor/Indexing_input_text/")


library(data.table)

#############################   function     #############################################


illiminaprocessor_prep_function <- function (index_name_value,
                                   location_index_plate_value,
                                               i7_index_value,
                                               i5_index_value,
                                          species_names_value,
                                         path_raw_reads_value,
                                 path_illuminaprocessor_value,
                                        path_clean_data_value,
                                            path_config_value,
                                           number_cores_value = 8,
                                       path_trimmomatic_value) {

# input from user

index_name <- index_name_value
location_index_plate <- location_index_plate_value
i7_index <- i7_index_value
i5_index <- i5_index_value
species_names <- species_names_value
path_raw_reads <- path_raw_reads_value 

path_illuminaprocessor <- path_illuminaprocessor_value
path_clean_data <- path_clean_data_value
path_config <- path_config_value
number_cores <- number_cores_value
path_trimmomatic <- path_trimmomatic_value

master_dir <- getwd()

# read tables

index_name_df <- read.table(index_name, header=T, sep = "\t", stringsAsFactors = F)
location_index_plate_df <- read.table(location_index_plate, header=T, sep = "\t", stringsAsFactors = F)
i7_index_df <- read.table(i7_index, header=T, sep = "\t", stringsAsFactors = F)
i5_index_df <- read.table(i5_index, header=T, sep = "\t", stringsAsFactors = F)
species_names_df <- read.table(species_names, header=T, sep = "\t", stringsAsFactors = F)

# get list of files in pool folder

setwd(path_raw_reads)

pool_input_files <- list.files(path = ".", pattern="*fastq.gz$",full.names=T, ignore.case=T)
name_of_pools_1 <- gsub("./", "",pool_input_files)
name_of_pools_2 <- gsub("_R1_001.fastq.gz", "", name_of_pools_1)
name_of_pools_3 <- unique(gsub("_R2_001.fastq.gz", "", name_of_pools_2))

plate_pool_locations_raw <- gsub("^.*-([A-Z0-9][0-9]+).*$", "\\1", name_of_pools_3)

setwd(master_dir)

# data frame

pool_names_location_DT <-data.table(t(rbind(name_of_pools_3, plate_pool_locations_raw)), key="plate_pool_locations_raw")
pool_names_location_df <- as.data.frame(pool_names_location_DT, stringsAsFactors = F)

# get name of index

i7_i5_list <- list()

for( i in 1:length(plate_pool_locations_raw)){

                                 i7_i5_list[[i]] <- subset(location_index_plate_df, plate_pool_locations == plate_pool_locations_raw[i])

                                             }

i7_i5_selected_df <- do.call(rbind, i7_i5_list)

# reorder and match 

pool_names_location_order_df <- pool_names_location_df[match(i7_i5_selected_df$plate_pool_locations, pool_names_location_df$plate_pool_locations),]
pool_names_location_2_df <- cbind(pool_names_location_order_df, i7_i5_selected_df)
pool_names_location_2_df <- pool_names_location_2_df[,!duplicated(colnames(pool_names_location_2_df))]
pool_names_location_2_DT <- data.table(pool_names_location_2_df)

# get the index sequence

i7_index_in_pool <- as.character(pool_names_location_2_DT$i7)
i5_index_in_pool <- as.character(pool_names_location_2_DT$i5)
pool_location_for_names <- as.character(pool_names_location_2_DT$plate_pool_locations)

## i7

seq_i7_present_in_pool <- list()

for (i in 1:length(i7_index_in_pool)) {

                        seq_i7_present_in_pool[[i]] <-  subset(i7_index_df, i7 == i7_index_in_pool[i], select= c(i7, i7_seq))

                                      }

seq_i7_present_in_pool_df <- do.call(rbind, seq_i7_present_in_pool)
seq_i7_present_in_pool_df <- seq_i7_present_in_pool_df[match(pool_names_location_2_df$i7, seq_i7_present_in_pool_df$i7),]

## i5

seq_i5_present_in_pool <- list()

for (i in 1:length(i5_index_in_pool)) {

                        seq_i5_present_in_pool[[i]] <-  subset(i5_index_df, i5 == i5_index_in_pool[i], select= c(i5, i5_seq))

                                      }

seq_i5_present_in_pool_df <- do.call(rbind, seq_i5_present_in_pool)
seq_i5_present_in_pool_df <- seq_i5_present_in_pool_df[match(pool_names_location_2_df$i5, seq_i5_present_in_pool_df$i5),]


## get species name

species_name_in_pool <- list ()

for (i in 1:length(pool_location_for_names)) {

                        species_name_in_pool[[i]] <-  subset(species_names_df, plate_pool_locations == pool_location_for_names[i], select= c(plate_pool_locations, Genus_and_species_voucher))

                                      }

species_name_in_pool_df <- do.call(rbind, species_name_in_pool)
species_name_in_pool_df <- species_name_in_pool_df[match(pool_names_location_2_df$plate_pool_locations, species_name_in_pool_df$plate_pool_locations),]


## 

pool_names_location_3_DT <- cbind(pool_names_location_2_DT, seq_i7_present_in_pool_df, seq_i5_present_in_pool_df, species_name_in_pool_df)
pool_names_location_3_df <- as.data.frame(pool_names_location_3_DT, stringsAsFactors =F)

## write lines

adapter_line <- "[adapters]"
adapter_i7_line <- "i7:GATCGGAAGAGCACACGTCTGAACTCCAGTCAC*ATCTCGTATGCCGTCTTCTGCTTG"
adapter_i5_line <- "i5:AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT*GTGTAGATCTCGGTGGTCGCCGTATCATT"

tag_line <- "[tag sequences]"

## section for tag sequences

i7_line_1 <- paste("i7-", pool_names_location_3_df$i7, ":", pool_names_location_3_df$i7_seq, sep="")
i5_line_2 <- paste("i5-", pool_names_location_3_df$i5, ":", pool_names_location_3_df$i5_seq, sep="")

##

tag_map_line <- "[tag map]"

## pools line

pools_line <- paste(name_of_pools_3, "_", pool_names_location_3_df$i7_seq, ":i7-", pool_names_location_3_df$i7, ",i5-", pool_names_location_3_df$i5, sep="")

# names line

names_top_line <- "[names]"
names_line <- paste(name_of_pools_3, "_", pool_names_location_3_df$i7_seq, ":", pool_names_location_3_df$Genus_and_species_voucher, sep="")

## create the ext file

setwd(master_dir)

illumina_clean_config_file<-file("config_file_illuminaprocessor.config")
writeLines(          c(adapter_line,
                       adapter_i7_line,
                       adapter_i5_line,
                       "\n",
                       tag_line,
                       i7_line_1,
                       i5_line_2,
                       "\n",
                       tag_map_line,
                       pools_line,
                       "\n",
                       names_top_line,
                       names_line), 
                           illumina_clean_config_file)
close(illumina_clean_config_file)


##########################################################################################
## rename files
## chage the names of the pools by adding i7 sequence

# identify the folders
current.folder <- path_raw_reads
new.folder <- paste(master_dir, "/process_pools", sep = "")

dir.create(new.folder)

# find the files that you want
list.of.files <- list.files(current.folder, "*fastq.gz$",full.names=T)
 
# copy the files to the new folder
file.copy(list.of.files, new.folder)

# new working directory

file_original_names <- list.files(path = new.folder, pattern="*fastq.gz$",full.names=T, ignore.case=T)
file_pool_locations <- gsub("^.*-([A-Z0-9][0-9]+).*$", "\\1", file_original_names)

file_original_names_df <- as.data.frame(cbind(file_original_names, file_pool_locations), stringsAsFactors = F)
file_original_names_df_split <- split(file_original_names_df, file_original_names_df$file_pool_locations)

### for replace names

## R1 and R2

new_names_with_seq_R1 <- list ()
new_names_with_seq_R2 <- list ()

for(i in 1:length(file_original_names_df_split))  {

temp <- file_original_names_df_split[[i]]

name_of_location <- unique(temp$file_pool_locations)
seq_of_location <- as.character(subset(pool_names_location_3_df, plate_pool_locations == name_of_location, select= c(i7_seq)))
plate_code_of_location <- as.character(subset(pool_names_location_3_df, plate_pool_locations == name_of_location, select= c(plate_pool_locations )))
new_name_R1 <- gsub("_R1_001", paste("_", seq_of_location, "_R1_001", sep=""), temp$file_original_names[1])
new_name_R2 <- gsub("_R2_001", paste("_", seq_of_location, "_R2_001", sep=""), temp$file_original_names[2])
new_names_with_seq_R1[[i]] <- rbind(new_name_R1, seq_of_location, plate_code_of_location)
new_names_with_seq_R2[[i]] <- rbind(new_name_R2, seq_of_location, plate_code_of_location)

                                                 }

new_names_with_seq_R1_df <- data.frame(matrix(unlist(new_names_with_seq_R1), nrow=nrow(file_original_names_df)/2, byrow=T),stringsAsFactors=FALSE)
new_names_with_seq_R2_df <- data.frame(matrix(unlist(new_names_with_seq_R2), nrow=nrow(file_original_names_df)/2, byrow=T),stringsAsFactors=FALSE)

names(new_names_with_seq_R1_df) <-  c("new_names", "i7_seq", "plate_pool_locations")
names(new_names_with_seq_R2_df) <-  c("new_names", "i7_seq", "plate_pool_locations")

# split orignal names by pair

file_original_names_R1_df <- subset(file_original_names_df, grepl("R1_001", file_original_names_df$file_original_names))
file_original_names_R2_df <- subset(file_original_names_df, grepl("R2_001", file_original_names_df$file_original_names))

# reorder and match 

file_original_names_R1_df<- file_original_names_R1_df[match(new_names_with_seq_R1_df$plate_pool_locations, file_original_names_R1_df$file_pool_locations),]
file_original_names_R2_df<- file_original_names_R2_df[match(new_names_with_seq_R2_df$plate_pool_locations, file_original_names_R2_df$file_pool_locations),]

# merge together and rename files

file_original_names_R1_merged_df <- cbind(file_original_names_R1_df, new_names_with_seq_R1_df)
file_original_names_R2_merged_df <- cbind(file_original_names_R2_df, new_names_with_seq_R2_df)

file_original_names_final_df <- rbind(file_original_names_R1_merged_df, file_original_names_R2_merged_df)
file.rename(file_original_names_final_df$file_original_names, file_original_names_final_df$new_names)

##########################################################################################
## comands to run 

setwd(master_dir)

first_part_comamd <- "illumiprocessor --input " 

# raw reads path  make a folder per pool 

path_illuminaprocessor_user <- path_illuminaprocessor

out_part_comamd <- " --output "

# clean data path

path_clean_data_user <- path_clean_data

config_part_comamd <- " --config "

# configuration file path a

path_config_user <- path_config

# config file Pool1-B1

name_config_file <- "config_file_illuminaprocessor.config"  

cores_part_comamd <- paste(" --cores ", number_cores, " --trimmomatic ", sep ="")

# path to trimmomatic path to the jar file

path_trimmomatic_user <- path_trimmomatic

last_part_command <- paste(" --r1-pattern ", "&", " --r2-pattern ", " @" , sep = "")


# list of commands

out_comands <- paste(first_part_comamd, path_illuminaprocessor_user, out_part_comamd, path_clean_data_user,
             config_part_comamd, path_config_user, name_config_file, cores_part_comamd,
             path_trimmomatic_user, last_part_command, sep="")


## create the ext file

command_file<-file("illuminaprocessor_commands.txt")
writeLines(out_comands, command_file)
close(command_file)

# write diagnostic tables

write.table(file_original_names_final_df, file="new_and_old_pool_file_names.txt", sep="\t")
write.table(pool_names_location_3_df, file="master_guide_file_names.txt", sep="\t")

}

#############################   END of function     ######################################


illiminaprocessor_prep_function (index_name_value = "~/Dropbox/trial_illuminaprocessor/Indexing_input_text/index_names.txt",
                       location_index_plate_value = "~/Dropbox/trial_illuminaprocessor/Indexing_input_text/location_index_plate.txt",
                                   i7_index_value = "~/Dropbox/trial_illuminaprocessor/Indexing_input_text/i7_index.txt",
                                   i5_index_value = "~/Dropbox/trial_illuminaprocessor/Indexing_input_text/i5_index.txt",
                              species_names_value = "~/Dropbox/trial_illuminaprocessor/Indexing_input_text/species_sequenced_exome_capture.txt",
                             path_raw_reads_value = "~/Dropbox/trial_illuminaprocessor/Indexing_input_text/raw_reads/",
                     path_illuminaprocessor_value = "~/Users/perry.wood/Desktop/trial_illuminaprocessor",
                            path_clean_data_value = "~/Users/perry.wood/Desktop/trial_illuminaprocessor/output",
                                path_config_value = "~/Users/perry.wood/Desktop/trial_illuminaprocessor/",
                               number_cores_value = 8,
                           path_trimmomatic_value = "~/Users/perry.wood/anaconda/bin/trimmomatic-0.35.jar")

# files created:

# config_file_illuminaprocessor.config
# illuminaprocessor_commands.txt
# master_guide_file_names.txt
# new_and_old_pool_file_names.txt


#################  de bug data

index_name_value <- "/Volumes/Genomic2/Exome_capture/input_files/Indexing_input_text/index_names.txt"
location_index_plate_value <- "/Volumes/Genomic2/Exome_capture/input_files/Indexing_input_text/location_index_plate.txt"
i7_index_value <- "/Volumes/Genomic2/Exome_capture/input_files/Indexing_input_text/i7_index.txt"
i5_index_value <- "/Volumes/Genomic2/Exome_capture/input_files/Indexing_input_text/i5_index.txt"
species_names_value <- "/Volumes/Genomic2/Exome_capture/input_files/Indexing_input_text/species_sequenced_exome_capture.txt"
path_raw_reads_value <- "/Volumes/Genomic2/Exome_capture/pools"
path_illuminaprocessor_value <- "/Users/perry.wood/Desktop/trial_illuminaprocessor"
path_clean_data_value <- "/Users/perry.wood/Desktop/trial_illuminaprocessor/output"
path_config_value <- "/Users/perry.wood/Desktop/trial_illuminaprocessor/"
number_cores_value <- 8
path_trimmomatic_value <- "/Users/perry.wood/anaconda/bin/trimmomatic-0.35.jar" 







