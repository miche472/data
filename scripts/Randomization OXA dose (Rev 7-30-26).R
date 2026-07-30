
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
ID <- c("S1","S2","S3","S4","S5","S6","S7","S8", "S9")

DOSE <- c("Dose A",
          "Dose B",
          "Dose C",
          "Dose D",
          "Dose E")

# Build standard 5 × 5 Latin square
latin_square <- tribble(
  ~INJECT_DAY_1, ~INJECT_DAY_2, ~INJECT_DAY_3, ~INJECT_DAY_4, ~INJECT_DAY_5,
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
  latin_square[sample(1:5, 4, replace = FALSE), ])

# Shuffle the eight sequences and assign IDs
design <- design %>%
  slice_sample(n = 9) %>%
  mutate(ID = ID, .before = 1)
design

# Actual mouse IDs and dose values
actual_ID <- c(3731, 3732, 3733, 3735,3737, 3739, 3740, 3741, 3738)

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
    across(starts_with("INJECT"), ~ dose_key[.x]))
design_final

#write_csv(x = design_final, "../data/design_final_7-17.csv")
write_csv(x = design_final, "../data/design_final_7-30.csv")

# Change format to long to allow for later joining with other data frames
# Convert to long format --> One row per mouse × dose combination

design_long <- design_final %>%
  pivot_longer(
    cols = starts_with("INJECT_DAY"),
    names_to = "INJECT_DAY",
    values_to = "DOSE") %>%
  mutate(PERIOD = readr::parse_number(INJECT_DAY)) %>%
  arrange(ID, INJECT_DAY)
design_long

#write_csv(x = design_final, "../data/design_long_7-17.csv")
write_csv(x = design_final, "../data/design_long_7-30.csv")