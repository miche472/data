
#Group assignment for OXA dose response study (Sept 2026) ####

#we have 8 NZO female mice that have been cannulated in the 3rd ventricle
#Mice will all receive one injection of each dose of OXA and an injection of vehicle
#There is a 72hr washout period between injections
#Mice will be recording in Sable for the duration of the study

#Started: 7-13-26
#Revised:8-27-26 (tailored for second iteration of experiment)

library(tidyverse)

set.seed(123)

#-----------------------------------------------------.#
## 1. Establish variables for IDs and treatments ####
#-----------------------------------------------------.#

# Placeholder IDs and Doses
ID <- c("S1","S2","S3","S4","S5","S6","S7","S8")

DOSE <- c("Dose A",
          "Dose B",
          "Dose C")

# Actual mouse IDs and dose values
actual_ID <- c(3731, 3732, 3733, 3735, 3737, 3738, 3739, 3741)

actual_DOSE <- c(
  "Vehicle",
  "500pmol",
  "2000pmol")

#-----------------------------------------------------.#
## 2. Build standard 3 × 3 Latin square ####
#-----------------------------------------------------.#
latin_square <- tribble(
  ~INJECT_DAY_1, ~INJECT_DAY_2, ~INJECT_DAY_3,
  DOSE[1],  DOSE[2],  DOSE[3],  
  DOSE[2],  DOSE[3],  DOSE[1],
  DOSE[3],  DOSE[1],  DOSE[2])


# Randomize order of the three sequences
latin_square <- latin_square %>%
  slice_sample(n = 3)

#-----------------------------------------------------.#
## 3. Modify latin square for 8 mice ####
#-----------------------------------------------------.#
# Randomly select three sequences to repeat 
    # We have 3 doses so 3 possible sequences, but we have n=8 mice
    #This step is necessary because we don't just have 3x3 design
design <- bind_rows(
  latin_square,
  latin_square[sample(1:3, 5, replace = TRUE), ])

# Shuffle the eight sequences and assign IDs
design <- design %>%
  slice_sample(n = 8) %>%
  mutate(ID = ID, .before = 1)

design

#-----------------------------------------------------.#
## 4. Replace placeholder variable names ####
#-----------------------------------------------------.#
# Replace placeholder ID and DOSE with study specific values
id_key <- setNames(sample(actual_ID), ID)
dose_key <- setNames(actual_DOSE, DOSE)

design_final <- design %>%
  mutate(
    ID = id_key[ID],
    across(starts_with("INJECT"), ~ dose_key[.x]))
design_final

write_csv(x = design_final, "../data/design_final_8-27.csv")
#write_csv(x = design_final, "../data/design_final_8-30.csv")

#-----------------------------------------------------.#
## 5. Change to long format #### 
#-----------------------------------------------------.#
#Allows for later joining with other data frames
# Convert to long format --> One row per mouse × dose combination

design_long <- design_final %>%
  pivot_longer(
    cols = starts_with("INJECT_DAY"),
    names_to = "INJECT_DAY",
    values_to = "DOSE") %>%
  mutate(PERIOD = readr::parse_number(INJECT_DAY)) %>%
  arrange(ID, INJECT_DAY)
design_long

write_csv(x = design_final, "../data/design_long_8-27.csv")
#write_csv(x = design_final, "../data/design_long_8-30.csv")