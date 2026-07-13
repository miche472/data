
#Group assignment for OXA dose response study (July 2026)

#we have 8 NZO female mice that have been cannulated in the 3rd ventricle
#Mice will all receive one injection of each dose of OXA and an injection of vehicle
#There is a 72hr washout period between injections
#Mice will be recording in Sable for the duration of the study

#Started: 7-13-26
#Revised:

library(tidyverse)

set.seed(123)

# Placeholder IDs and Doses
ID <- c("S1","S2","S3","S4","S5","S6","S7","S8")

DOSE <- c("Dose A",
          "Dose B",
          "Dose C",
          "Dose D",
          "Dose E")

# Build standard 5 × 5 Latin square
latin_square <- tribble(
  ~Period1, ~Period2, ~Period3, ~Period4, ~Period5,
  DOSE[1],  DOSE[2],  DOSE[3],  DOSE[4],  DOSE[5],
  DOSE[2],  DOSE[3],  DOSE[4],  DOSE[5],  DOSE[1],
  DOSE[3],  DOSE[4],  DOSE[5],  DOSE[1],  DOSE[2],
  DOSE[4],  DOSE[5],  DOSE[1],  DOSE[2],  DOSE[3],
  DOSE[5],  DOSE[1],  DOSE[2],  DOSE[3],  DOSE[4])


# Randomize order of the five sequences
latin_square <- latin_square %>%
  slice_sample(n = 5)

# Randomly select three sequences to repeat 
    # We have 5 doses so 5 possible sequences, but we have n=8 mice
    #This step is necessary because we don't just have 5x5 design
design <- bind_rows(
  latin_square,
  latin_square[sample(1:5, 3, replace = FALSE), ])

# Shuffle the eight sequences and assign IDs
design <- design %>%
  slice_sample(n = 8) %>%
  mutate(ID = ID, .before = 1)
design

# Actual mouse IDs and dose values
actual_ID <- c(3731, 3732, 3733, 3735,3737, 3739, 3740, 3741)

actual_DOSE <- c(
  "Vehicle",
  "125pmol",
  "250pmol",
  "500pmol",
  "1000pmol")

# Replace placeholder ID and DOSE with study specific values
id_key <- setNames(sample(actual_ID), ID)
dose_key <- setNames(actual_DOSE, DOSE)

design_final <- design %>%
  mutate(
    ID = id_key[ID],
    across(starts_with("Period"), ~ dose_key[.x]))
design_final

# Change format to long to allow for later joining with other data frames
# Convert to long format --> One row per mouse × dose combination

design_long <- design_final %>%
  pivot_longer(
    cols = starts_with("Period"),
    names_to = "PERIOD",
    values_to = "DOSE") %>%
  mutate(PERIOD = readr::parse_number(PERIOD)) %>%
  arrange(ID, PERIOD)
design_long
