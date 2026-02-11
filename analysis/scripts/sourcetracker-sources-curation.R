library(tidyverse)
library(janitor)

# As downloaded the 20221111
st_data_raw <- read_tsv("./data/sourcetracker/cdata/kapk-biomes__combined.tsv") |>
  clean_names() |>
  # select(studies, run_accession, biome, fastq_ftp, environment_biome,sample_name,sample_desc, geo_loc_name) |>
  mutate(environment_biome = tolower(environment_biome)) |>
  mutate(environment_biome = case_when(
    environment_biome == "fresh water envo:00002011" ~ "",
    environment_biome == "boreal plains" ~ "boreal_plains",
    environment_biome == "semi-arid plateaux" ~ "semi_arid_plateau",
    environment_biome == "temperate grassland biome" ~ "temperate grassland",
    environment_biome == "woodland" ~ "woodlands",
    environment_biome == "cropland biome" ~ "cropland",
    is.na(environment_biome) ~ "",
    TRUE ~ environment_biome
  ))

# Overview of biomes
st_data_raw |>
  group_by(biome) |>
  count(sort = TRUE) |>
  knitr::kable()

# |biome                                                                  |    n|
# |:----------------------------------------------------------------------|----:|
# |root:Environmental:Aquatic:Marine                                      | 1737|
# |root:Environmental:Aquatic:Marine:Oceanic                              | 1383|
# |root:Host-associated:Animal:Digestive system:Fecal                     |  407|
# |root:Environmental:Aquatic:Freshwater:Lake                             |  375|
# |root:Host-associated:Mammals:Gastrointestinal tract:Intestine:Fecal    |  218|
# |root:Host-associated:Mammals:Digestive system:Fecal                    |  169|
# |root:Environmental:Aquatic:Marine:Intertidal zone:Estuary              |  116|
# |root:Environmental:Terrestrial:Soil                                    |   84|
# |root:Host-associated:Plants                                            |   54|
# |root:Host-associated:Mammals:Digestive system                          |   47|
# |root:Environmental:Aquatic:Marine:Sediment                             |   44|
# |root:Environmental:Aquatic:Marine:Coastal                              |   38|
# |root:Environmental:Terrestrial:Soil:Permafrost                         |   32|
# |root:Host-associated:Mammals:Digestive system:Large intestine:Fecal    |   32|
# |root:Environmental:Aquatic:Freshwater                                  |   30|
# |root:Environmental:Aquatic:Estuary                                     |   24|
# |root:Environmental:Aquatic:Marine:Intertidal zone:Coral reef           |   23|
# |root:Host-associated:Plants:Rhizosphere                                |   23|
# |root:Environmental:Terrestrial:Soil:Forest soil                        |   22|
# |root:Environmental:Aquatic:Marine:Intertidal zone:Salt marsh           |   21|
# |root:Host-associated:Mammals:Digestive system:Stomach:Rumen            |   20|
# |root:Environmental:Terrestrial:Soil:Agricultural                       |   18|
# |root:Environmental:Aquatic:Freshwater:Sediment                         |   16|
# |root:Environmental:Aquatic:Freshwater:Drinking water                   |   15|
# |root:Environmental:Aquatic:Freshwater:Lentic:Sediment                  |   15|
# |root:Environmental:Aquatic:Freshwater:Lotic:Sediment                   |   12|
# |root:Environmental:Terrestrial:Soil:Sand                               |   12|
# |root:Environmental:Aquatic:Marine:Oceanic:Oil-contaminated             |   11|
# |root:Host-associated:Mammals:Digestive system:Foregut:Rumen            |   11|
# |root:Environmental:Aquatic:Freshwater:Drinking water:Delivery networks |   10|
# |root:Environmental:Aquatic:Estuary:Sediment                            |    9|
# |root:Environmental:Terrestrial:Soil:Grasslands                         |    8|
# |root:Environmental:Aquatic:Freshwater:Groundwater                      |    7|
# |root:Environmental:Aquatic:Freshwater:Lotic:Low land river systems     |    7|
# |root:Host-associated:Animal                                            |    7|
# |root:Environmental:Aquatic:Freshwater:Lentic                           |    6|
# |root:Environmental:Aquatic:Marine:Hydrothermal vents:Diffuse flow      |    6|
# |root:Environmental:Aquatic:Marine:Marginal Sea                         |    6|
# |root:Environmental:Aquatic:Marine:Oceanic:Photic zone                  |    6|
# |root:Environmental:Aquatic:Freshwater:Pond:Sediment                    |    5|
# |root:Environmental:Aquatic:Marine:Intertidal zone:Microbialites        |    4|
# |root:Environmental:Terrestrial:Soil:Mine                               |    4|
# |root:Host-associated:Birds:Digestive system:Digestive tube:Cecum       |    4|
# |root:Host-associated:Insecta                                           |    4|
# |root:Host-associated:Plants:Phylloplane                                |    4|
# |root:Environmental:Aquatic:Marine:Pelagic                              |    3|
# |root:Host-associated:Birds:Digestive system                            |    3|
# |root:Host-associated:Plants:Rhizosphere:Soil                           |    3|
# |root:Environmental:Aquatic:Marine:Oceanic:Sediment                     |    2|
# |root:Host-associated:Algae                                             |    2|
# |root:Host-associated:Invertebrates                                     |    2|
# |root:Host-associated:Mammals:Digestive system:Large intestine:Cecum    |    2|
# |root:Environmental:Aquatic:Freshwater:Groundwater:Biofilm              |    1|
# |root:Environmental:Aquatic:Freshwater:Ice:Glacial lake                 |    1|
# |root:Environmental:Aquatic:Marine:Intertidal zone:Oil-contaminated     |    1|
# |root:Environmental:Aquatic:Marine:Oil field:bore hole                  |    1|
# |root:Environmental:Aquatic:Marine:Oil-contaminated sediment            |    1|
# |root:Host-associated:Mammals:Digestive system:Midgut                   |    1|
# |root:Host-associated:Mammals:Respiratory system                        |    1|



# REFINEMENT --------------------------------------------------------------
set.seed(42)
##########################################################
## root:Environmental:Aquatic:Marine
##########################################################

# root:Environmental:Aquatic:Marine  ---------------------
st_data_raw |>
  filter(biome == "root:Environmental:Aquatic:Marine") |>
  group_by(biome, environment_biome, geo_loc_name_1) |>
  select(run_accession, biome, environment_biome, geo_loc_name_1) |>
  add_count() |>
  sample_n(1) |>
  arrange(desc(n)) |>
  filter(n >= 25) |>
  knitr::kable()

# We will pick 50 random samples from each combination
# |run_accession |biome                             |environment_biome            |geo_loc_name_1  |   n|
# |:-------------|:---------------------------------|:----------------------------|:---------------|---:|
# |SRR5788033    |root:Environmental:Aquatic:Marine |ocean_biome                  |Atlantic Ocean  | 286|*
# |SRR12479958   |root:Environmental:Aquatic:Marine |                             |Indian Ocean    | 260|*
# |ERR2752157    |root:Environmental:Aquatic:Marine |marine biome (envo:00000447) |NA              | 213|
# |SRR12479853   |root:Environmental:Aquatic:Marine |                             |Atlantic Ocean  | 202|
# |SRR5788243    |root:Environmental:Aquatic:Marine |ocean_biome                  |Pacific Ocean   | 184|*
# |ERR2762143    |root:Environmental:Aquatic:Marine |polar biome                  |NA              | 121|*
# |SRR12479984   |root:Environmental:Aquatic:Marine |                             |Pacific Ocean   | 106|
# |SRR2135747    |root:Environmental:Aquatic:Marine |                             |not applicable  |  71|
# |ERR2196984    |root:Environmental:Aquatic:Marine |westerlies biome             |NA              |  59|
# |SRR3933188    |root:Environmental:Aquatic:Marine |                             |Canada:Atlantic |  58|*
# |SRR2053279    |root:Environmental:Aquatic:Marine |brackish water               |Baltic Sea      |  25|*


set.seed(42)
sel_aquatic_marine <- st_data_raw |>
  filter(biome == "root:Environmental:Aquatic:Marine") |>
  group_by(biome, environment_biome, geo_loc_name_1) |>
  select(run_accession, biome, environment_biome, geo_loc_name_1) |>
  add_count() |>
  sample_n(1) |>
  arrange(desc(n)) |>
  filter(n >= 25) |>
  ungroup() |>
  mutate(r = row_number()) |>
  filter(r %in% (c(1, 2, 5, 6, 10, 11))) |>
  select(biome, environment_biome, geo_loc_name_1) |>
  inner_join(st_data_raw) |>
  select(run_accession, biome, environment_biome, geo_loc_name_1) |>
  group_by(biome, environment_biome, geo_loc_name_1) |>
  sample_n(25) |>
  distinct() |>
  ungroup() |>
  select(run_accession, biome, environment_biome, geo_loc_name_1) |>
  mutate(biome_class = "root:Environmental:Aquatic:Marine")


# root:Environmental:Aquatic:Marine:Oceanic  ------------------------------
st_data_raw |>
  filter(biome == "root:Environmental:Aquatic:Marine:Oceanic") |>
  group_by(biome, environment_biome, geo_loc_name_1) |>
  select(run_accession, biome, environment_biome, geo_loc_name_1) |>
  add_count() |>
  sample_n(1) |>
  arrange(desc(n)) |>
  filter(n >= 25) |>
  knitr::kable()

# |run_accession |biome                                     |environment_biome            |geo_loc_name_1                                               |    n|
# |:-------------|:-----------------------------------------|:----------------------------|:------------------------------------------------------------|----:|
# |ERR1726828    |root:Environmental:Aquatic:Marine:Oceanic |marine biome (envo:00000447) |NA                                                           | 1127|
# |SRR5720304    |root:Environmental:Aquatic:Marine:Oceanic |ocean biome                  |Pacific Ocean: North Pacific Subtropical Gyre, Station ALOHA |   68|*
# |SRR5720320    |root:Environmental:Aquatic:Marine:Oceanic |ocean biome                  |Atlantic Ocean: Sargasso Sea, BATS                           |   62|*
# |SRR3933388    |root:Environmental:Aquatic:Marine:Oceanic |                             |Canada:Atlantic                                              |   52|
# |ERR2206795    |root:Environmental:Aquatic:Marine:Oceanic |brackish water               |Baltic Sea                                                   |   38|

set.seed(42)
sel_aquatic_marine_oceanic <- st_data_raw |>
  filter(biome == "root:Environmental:Aquatic:Marine:Oceanic") |>
  group_by(biome, environment_biome, geo_loc_name_1) |>
  select(run_accession, biome, environment_biome, geo_loc_name_1) |>
  add_count() |>
  sample_n(1) |>
  arrange(desc(n)) |>
  filter(n >= 25) |>
  ungroup() |>
  mutate(r = row_number()) |>
  filter(r %in% (c(2, 3))) |>
  select(biome, environment_biome, geo_loc_name_1) |>
  inner_join(st_data_raw) |>
  select(run_accession, biome, environment_biome, geo_loc_name_1) |>
  group_by(biome, environment_biome, geo_loc_name_1) |>
  sample_n(25) |>
  distinct() |>
  ungroup() |>
  select(run_accession, biome, environment_biome, geo_loc_name_1) |>
  mutate(biome_class = "root:Environmental:Aquatic:Marine")


# root:Environmental:Aquatic:Marine:Intertidal zone:Estuary ---------------
st_data_raw |>
  filter(biome == "root:Environmental:Aquatic:Marine:Intertidal zone:Estuary") |>
  group_by(biome, environment_biome, geo_loc_name_1) |>
  select(run_accession, biome, environment_biome, geo_loc_name_1) |>
  add_count() |>
  sample_n(1) |>
  arrange(desc(n)) |>
  filter(n >= 25) |>
  knitr::kable()

# |run_accession |biome                                                     |environment_biome |geo_loc_name_1                  |  n|
# |:-------------|:---------------------------------------------------------|:-----------------|:-------------------------------|--:|
# |SRR3719696    |root:Environmental:Aquatic:Marine:Intertidal zone:Estuary |                  |NA                              | 39|
# |SRR5468367    |root:Environmental:Aquatic:Marine:Intertidal zone:Estuary |                  |USA:Columbia River Estuary, USA | 36|


# root:Environmental:Aquatic:Estuary  -------------------------------------

st_data_raw |>
  filter(biome == "root:Environmental:Aquatic:Estuary") |>
  group_by(biome, environment_biome, geo_loc_name_1) |>
  select(run_accession, biome, environment_biome, geo_loc_name_1) |>
  add_count() |>
  sample_n(1) |>
  arrange(desc(n)) |>
  knitr::kable()

# |run_accession |biome                              |environment_biome |geo_loc_name_1                  |  n|
# |:-------------|:----------------------------------|:-----------------|:-------------------------------|--:|
# |SRR5468406    |root:Environmental:Aquatic:Estuary |                  |USA:Columbia River Estuary, USA | 24|

# root:Environmental:Aquatic:Estuary:Sediment -----------------------------
st_data_raw |>
  filter(biome == "root:Environmental:Aquatic:Estuary:Sediment") |>
  group_by(biome, environment_biome, geo_loc_name_1) |>
  select(run_accession, biome, environment_biome, geo_loc_name_1) |>
  add_count() |>
  sample_n(1) |>
  arrange(desc(n)) |>
  knitr::kable()

#   |run_accession |biome                                       |environment_biome |geo_loc_name_1                  |  n|
#   |:-------------|:-------------------------------------------|:-----------------|:-------------------------------|--:|
#   |SRR5808824    |root:Environmental:Aquatic:Estuary:Sediment |marsh             |USA: Chandeleur Islands         |  3|
#   |SRR5808813    |root:Environmental:Aquatic:Estuary:Sediment |subtidal zone     |USA: Chandeleur Islands         |  3|
#   |SRR5468372    |root:Environmental:Aquatic:Estuary:Sediment |                  |USA:Columbia River Estuary, USA |  2|
#   |ERR2730831    |root:Environmental:Aquatic:Estuary:Sediment |estuary           |Canada                          |  1|

set.seed(42)
sel_aquatic_estuary <- st_data_raw |>
  filter(biome == "root:Environmental:Aquatic:Marine:Intertidal zone:Estuary" |
    biome == "root:Environmental:Aquatic:Estuary:Sediment" |
    biome == "root:Environmental:Aquatic:Estuary") |>
  group_by(biome, environment_biome, geo_loc_name_1) |>
  select(run_accession, biome, environment_biome, geo_loc_name_1) |>
  add_count() |>
  sample_n(1) |>
  arrange(desc(n)) |>
  ungroup() |>
  mutate(r = row_number()) |>
  select(biome, environment_biome, geo_loc_name_1) |>
  inner_join(st_data_raw) |>
  select(run_accession, biome, environment_biome, geo_loc_name_1) |>
  group_by(biome, environment_biome, geo_loc_name_1) |>
  sample_n(25, replace = TRUE) |>
  distinct() |>
  ungroup() |>
  select(run_accession, biome, environment_biome, geo_loc_name_1) |>
  mutate(
    biome_class = "root:Environmental:Aquatic:Marine",
    biome = "root:Environmental:Aquatic:Estuary"
  )

# root:Environmental:Aquatic:Marine:Sediment  -----------------------------
st_data_raw |>
  filter(biome == "root:Environmental:Aquatic:Marine:Sediment") |>
  group_by(biome, environment_biome, geo_loc_name_1) |>
  select(run_accession, biome, environment_biome, geo_loc_name_1) |>
  add_count() |>
  sample_n(1) |>
  arrange(desc(n)) |>
  knitr::kable()

# |run_accession |biome                                      |environment_biome |geo_loc_name_1                                            |  n|
# |:-------------|:------------------------------------------|:-----------------|:---------------------------------------------------------|--:|
# |SRR8943094    |root:Environmental:Aquatic:Marine:Sediment |                  |Sweden: Baltic Sea                                        | 12|
# |ERR2564009    |root:Environmental:Aquatic:Marine:Sediment |sediment          |NA                                                        |  9|
# |SRR12059197   |root:Environmental:Aquatic:Marine:Sediment |                  |Denmark                                                   |  7|
# |SRR5469033    |root:Environmental:Aquatic:Marine:Sediment |                  |USA:Hudson Canyon                                         |  5|
# |ERR1333180    |root:Environmental:Aquatic:Marine:Sediment |mine tailing pool |NA                                                        |  5|
# |SRR12059194   |root:Environmental:Aquatic:Marine:Sediment |                  |Sweden                                                    |  3|
# |SRR6193154    |root:Environmental:Aquatic:Marine:Sediment |                  |Puerto Rico:Bioluminescent Bay, La Paraguera, Puerto Rico |  1|
# |SRR6194934    |root:Environmental:Aquatic:Marine:Sediment |                  |USA:Eden Landing Ponds, San Francisco, CA, USA            |  1|
# |SRR5149189    |root:Environmental:Aquatic:Marine:Sediment |marine            |Atlantic Ocean                                            |  1|

set.seed(42)
sel_aquatic_marine_sediment <- st_data_raw |>
  filter(biome == "root:Environmental:Aquatic:Marine:Sediment") |>
  group_by(biome, environment_biome, geo_loc_name_1) |>
  select(run_accession, biome, environment_biome, geo_loc_name_1) |>
  add_count() |>
  sample_n(1) |>
  arrange(desc(n)) |>
  ungroup() |>
  mutate(r = row_number()) |>
  select(biome, environment_biome, geo_loc_name_1) |>
  inner_join(st_data_raw) |>
  select(run_accession, biome, environment_biome, geo_loc_name_1) |>
  group_by(biome, environment_biome, geo_loc_name_1) |>
  sample_n(25, replace = TRUE) |>
  distinct() |>
  ungroup() |>
  select(run_accession, biome, environment_biome, geo_loc_name_1) |>
  mutate(biome_class = "root:Environmental:Aquatic:Marine")


# root:Environmental:Aquatic:Marine:Coastal  ------------------------------
st_data_raw |>
  filter(biome == "root:Environmental:Aquatic:Marine:Coastal") |>
  group_by(biome, environment_biome, geo_loc_name_1) |>
  select(run_accession, biome, environment_biome, geo_loc_name_1) |>
  add_count() |>
  sample_n(1) |>
  arrange(desc(n)) |>
  knitr::kable()

# |run_accession |biome                                     |environment_biome |geo_loc_name_1 |  n|
# |:-------------|:-----------------------------------------|:-----------------|:--------------|--:|
# |ERR2281800    |root:Environmental:Aquatic:Marine:Coastal |aquatic           |Australia      | 12|
# |ERR1698985    |root:Environmental:Aquatic:Marine:Coastal |coastal bay       |NA             | 12|
# |ERR864078     |root:Environmental:Aquatic:Marine:Coastal |marine            |NA             |  6|
# |ERR1698993    |root:Environmental:Aquatic:Marine:Coastal |channel           |NA             |  4|
# |ERR864070     |root:Environmental:Aquatic:Marine:Coastal |estuary           |NA             |  2|
# |ERR864072     |root:Environmental:Aquatic:Marine:Coastal |river             |NA             |  2|

set.seed(42)
sel_aquatic_marine_coastal <- st_data_raw |>
  filter(biome == "root:Environmental:Aquatic:Marine:Coastal") |>
  group_by(biome, environment_biome, geo_loc_name_1) |>
  select(run_accession, biome, environment_biome, geo_loc_name_1) |>
  add_count() |>
  sample_n(1) |>
  arrange(desc(n)) |>
  ungroup() |>
  mutate(r = row_number()) |>
  select(biome, environment_biome, geo_loc_name_1) |>
  inner_join(st_data_raw) |>
  select(run_accession, biome, environment_biome, geo_loc_name_1) |>
  group_by(biome, environment_biome, geo_loc_name_1) |>
  sample_n(25, replace = TRUE) |>
  distinct() |>
  ungroup() |>
  select(run_accession, biome, environment_biome, geo_loc_name_1) |>
  mutate(biome_class = "root:Environmental:Aquatic:Marine")


# root:Environmental:Aquatic:Marine:Intertidal zone:Salt marsh ------------
st_data_raw |>
  filter(biome == "root:Environmental:Aquatic:Marine:Intertidal zone:Salt marsh") |>
  group_by(biome, environment_biome, geo_loc_name_1) |>
  select(run_accession, biome, environment_biome, geo_loc_name_1) |>
  add_count() |>
  sample_n(1) |>
  arrange(desc(n)) |>
  knitr::kable()

# |run_accession |biome                                                        |environment_biome                     |geo_loc_name_1                             |  n|
# |:-------------|:------------------------------------------------------------|:-------------------------------------|:------------------------------------------|--:|
# |ERR257713     |root:Environmental:Aquatic:Marine:Intertidal zone:Salt marsh |marine salt marsh biome envo:01000022 |NA                                         |  7|
# |SRR6201723    |root:Environmental:Aquatic:Marine:Intertidal zone:Salt marsh |                                      |USA:Emeryville, CA                         |  6|
# |ERR2215871    |root:Environmental:Aquatic:Marine:Intertidal zone:Salt marsh |saltmarsh                             |United Kingdom                             |  4|
# |SRR5621770    |root:Environmental:Aquatic:Marine:Intertidal zone:Salt marsh |                                      |USA: Emeryville, California                |  1|
# |SRR6192418    |root:Environmental:Aquatic:Marine:Intertidal zone:Salt marsh |                                      |USA:Alviso Ponds, San Francisco, CA, USA   |  1|
# |ERR257714     |root:Environmental:Aquatic:Marine:Intertidal zone:Salt marsh |sediment                              |United Kingdom                             |  1|
# |SRR948829     |root:Environmental:Aquatic:Marine:Intertidal zone:Salt marsh |temperate shelf and sea biome         |USA: Sippewissett Salt Marsh, Falmouth, MA |  1|


set.seed(42)
sel_aquatic_marine_intertidal_salt_marsh <- st_data_raw |>
  filter(biome == "root:Environmental:Aquatic:Marine:Intertidal zone:Salt marsh" | environment_biome == "marine salt marsh biome") |>
  group_by(biome, environment_biome, geo_loc_name_1) |>
  select(run_accession, biome, environment_biome, geo_loc_name_1) |>
  add_count() |>
  sample_n(1) |>
  arrange(desc(n)) |>
  ungroup() |>
  mutate(r = row_number()) |>
  select(biome, environment_biome, geo_loc_name_1) |>
  inner_join(st_data_raw) |>
  select(run_accession, biome, environment_biome, geo_loc_name_1) |>
  group_by(biome, environment_biome, geo_loc_name_1) |>
  sample_n(25, replace = TRUE) |>
  distinct() |>
  ungroup() |>
  select(run_accession, biome, environment_biome, geo_loc_name_1) |>
  mutate(
    biome_class = "root:Environmental:Aquatic:Marine",
    biome = "root:Environmental:Aquatic:Marine:Intertidal zone:Salt marsh"
  )


##########################################################
## root:Environmental:Aquatic:Freshwater
##########################################################

# root:Environmental:Aquatic:Freshwater:Lake ------------------------------
st_data_raw |>
  filter(biome == "root:Environmental:Aquatic:Freshwater:Lake") |>
  group_by(biome, environment_biome, geo_loc_name_1) |>
  select(run_accession, biome, environment_biome, geo_loc_name_1) |>
  add_count() |>
  sample_n(1) |>
  arrange(desc(n)) |>
  filter(n >= 25) |>
  knitr::kable()

# |:-------------|:------------------------------------------|:------------------|:--------------|--:|
# |ERR2848558    |root:Environmental:Aquatic:Freshwater:Lake |boreal_shield      |Canada         | 48|
# |ERR4869764    |root:Environmental:Aquatic:Freshwater:Lake |boreal_plains      |Canada         | 43|
# |ERR4869818    |root:Environmental:Aquatic:Freshwater:Lake |praires            |Canada         | 40|
# |ERR2848505    |root:Environmental:Aquatic:Freshwater:Lake |atlantic_maritime  |Canada         | 39|
# |ERR2848553    |root:Environmental:Aquatic:Freshwater:Lake |mixedwood_plains   |Canada         | 34|
# |ERR5954490    |root:Environmental:Aquatic:Freshwater:Lake |semi_arid_plateau  |Canada         | 29|
# |ERR5954400    |root:Environmental:Aquatic:Freshwater:Lake |montane cordillera |Canada         | 28|
# |ERR5954445    |root:Environmental:Aquatic:Freshwater:Lake |pacific maritime   |Canada         | 27|
# |ERR2848568    |root:Environmental:Aquatic:Freshwater:Lake |atlantic_highlands |Canada         | 25|


set.seed(42)
sel_aquatic_freshwater_lake <- st_data_raw |>
  filter(biome == "root:Environmental:Aquatic:Freshwater:Lake") |>
  group_by(biome, environment_biome, geo_loc_name_1) |>
  select(run_accession, biome, environment_biome, geo_loc_name_1) |>
  add_count() |>
  sample_n(1) |>
  arrange(desc(n)) |>
  filter(n >= 25) |>
  ungroup() |>
  mutate(r = row_number()) |>
  select(biome, environment_biome, geo_loc_name_1) |>
  inner_join(st_data_raw) |>
  select(run_accession, biome, environment_biome, geo_loc_name_1) |>
  group_by(biome, environment_biome, geo_loc_name_1) |>
  sample_n(25) |>
  distinct() |>
  ungroup() |>
  select(run_accession, biome, environment_biome, geo_loc_name_1) |>
  mutate(biome_class = "root:Environmental:Aquatic:Freshwater")



# root:Environmental:Aquatic:Freshwater -----------------------------------
st_data_raw |>
  filter(biome == "root:Environmental:Aquatic:Freshwater") |>
  group_by(biome, environment_biome, geo_loc_name_1) |>
  select(run_accession, biome, environment_biome, geo_loc_name_1) |>
  add_count() |>
  sample_n(1) |>
  arrange(desc(n)) |>
  knitr::kable()

# |run_accession |biome                                 |environment_biome        |geo_loc_name_1       |  n|
# |:-------------|:-------------------------------------|:------------------------|:--------------------|--:|
# |ERR3440684    |root:Environmental:Aquatic:Freshwater |industrial cooling water |United Arab Emirates | 30| # industrial cooling water

# root:Environmental:Aquatic:Freshwater:Sediment --------------------------
st_data_raw |>
  filter(biome == "root:Environmental:Aquatic:Freshwater:Sediment") |>
  group_by(biome, environment_biome, geo_loc_name_1) |>
  select(run_accession, biome, environment_biome, geo_loc_name_1) |>
  add_count() |>
  sample_n(1) |>
  arrange(desc(n)) |>
  knitr::kable()

# |run_accession |biome                                          |environment_biome |geo_loc_name_1        |  n|
# |:-------------|:----------------------------------------------|:-----------------|:---------------------|--:|
# |ERR3804375    |root:Environmental:Aquatic:Freshwater:Sediment |                  |NA                    | 14|
# |SRR13377137   |root:Environmental:Aquatic:Freshwater:Sediment |                  |Germany: Bremen, pond |  2|

# root:Environmental:Aquatic:Freshwater:Lentic:Sediment  -----------------
st_data_raw |>
  filter(biome == "root:Environmental:Aquatic:Freshwater:Lentic:Sediment") |>
  group_by(biome, environment_biome, geo_loc_name_1) |>
  select(run_accession, biome, environment_biome, geo_loc_name_1) |>
  add_count() |>
  sample_n(1) |>
  arrange(desc(n)) |>
  knitr::kable()

# |run_accession |biome                                                 |environment_biome |geo_loc_name_1 |  n|
# |:-------------|:-----------------------------------------------------|:-----------------|:--------------|--:|
# |ERR4334730    |root:Environmental:Aquatic:Freshwater:Lentic:Sediment |                  |NA             | 12|
# |ERR1725848    |root:Environmental:Aquatic:Freshwater:Lentic:Sediment |boreal            |NA             |  2|
# |SRR7067822    |root:Environmental:Aquatic:Freshwater:Lentic:Sediment |                  |USA: Nevada    |  1|

# root:Environmental:Aquatic:Freshwater:Lotic:Sediment  -------------------
st_data_raw |>
  filter(biome == "root:Environmental:Aquatic:Freshwater:Lotic:Sediment") |>
  group_by(biome, environment_biome, geo_loc_name_1) |>
  select(run_accession, biome, environment_biome, geo_loc_name_1) |>
  add_count() |>
  sample_n(1) |>
  arrange(desc(n)) |>
  knitr::kable()

# |run_accession |biome                                                |environment_biome |geo_loc_name_1                                 |  n|
# |:-------------|:----------------------------------------------------|:-----------------|:----------------------------------------------|--:|
# |SRR5468106    |root:Environmental:Aquatic:Freshwater:Lotic:Sediment |freshwater river  |China: Xiaomei River, Liaocheng                |  9|
# |SRR6284162    |root:Environmental:Aquatic:Freshwater:Lotic:Sediment |                  |Finland: Kymijoki River                        |  1|
# |SRR6284161    |root:Environmental:Aquatic:Freshwater:Lotic:Sediment |                  |USA: Hackensack River                          |  1|
# |SRR6195247    |root:Environmental:Aquatic:Freshwater:Lotic:Sediment |                  |USA:Eden Landing Ponds, San Francisco, CA, USA |  1|


biomes_f_s <- c("root:Environmental:Aquatic:Freshwater:Sediment", "root:Environmental:Aquatic:Freshwater:Lentic:Sediment", "root:Environmental:Aquatic:Freshwater:Lotic:Sediment")
set.seed(42)
sel_aquatic_freshwater_sediment <- st_data_raw |>
  filter(biome %in% biomes_f_s) |>
  group_by(biome, environment_biome, geo_loc_name_1) |>
  select(run_accession, biome, environment_biome, geo_loc_name_1) |>
  add_count() |>
  sample_n(1) |>
  arrange(desc(n)) |>
  ungroup() |>
  mutate(r = row_number()) |>
  select(biome, environment_biome, geo_loc_name_1) |>
  inner_join(st_data_raw) |>
  select(run_accession, biome, environment_biome, geo_loc_name_1) |>
  group_by(biome, environment_biome, geo_loc_name_1) |>
  sample_n(25, replace = TRUE) |>
  distinct() |>
  ungroup() |>
  select(run_accession, biome, environment_biome, geo_loc_name_1) |>
  mutate(biome_class = "root:Environmental:Aquatic:Freshwater")



# root:Environmental:Aquatic:Freshwater:Groundwater -----------------------
st_data_raw |>
  filter(biome == "root:Environmental:Aquatic:Freshwater:Groundwater") |>
  group_by(biome, environment_biome, geo_loc_name_1) |>
  select(run_accession, biome, environment_biome, geo_loc_name_1) |>
  add_count() |>
  sample_n(1) |>
  arrange(desc(n)) |>
  knitr::kable()

#   |run_accession |biome                                             |environment_biome |geo_loc_name_1 |  n|
#   |:-------------|:-------------------------------------------------|:-----------------|:--------------|--:|
#   |SRR5739200    |root:Environmental:Aquatic:Freshwater:Groundwater |                  |Denmark        |  6|
#   |SRR3496377    |root:Environmental:Aquatic:Freshwater:Groundwater |                  |NA             |  1|

# root:Environmental:Aquatic:Freshwater:Lotic:Low land river syste --------
st_data_raw |>
  filter(biome == "root:Environmental:Aquatic:Freshwater:Lotic:Low land river systems") |>
  group_by(biome, environment_biome, geo_loc_name_1) |>
  select(run_accession, biome, environment_biome, geo_loc_name_1) |>
  add_count() |>
  sample_n(1) |>
  arrange(desc(n)) |>
  knitr::kable()

# |run_accession |biome                                                              |environment_biome   |geo_loc_name_1       |  n|
# |:-------------|:------------------------------------------------------------------|:-------------------|:--------------------|--:|
# |SRR1522974    |root:Environmental:Aquatic:Freshwater:Lotic:Low land river systems |freshwater habitats |Brazil: Amazon river |  7|


# root:Environmental:Aquatic:Freshwater:Lentic ----------------------------
st_data_raw |>
  filter(biome == "root:Environmental:Aquatic:Freshwater:Lentic") |>
  group_by(biome, environment_biome, geo_loc_name_1) |>
  select(run_accession, biome, environment_biome, geo_loc_name_1) |>
  add_count() |>
  sample_n(1) |>
  arrange(desc(n)) |>
  knitr::kable()

#   |run_accession |biome                                        |environment_biome      |geo_loc_name_1 |  n|
#   |:-------------|:--------------------------------------------|:----------------------|:--------------|--:|
#   |ERR1474558    |root:Environmental:Aquatic:Freshwater:Lentic |water bacterioplankton |China          |  6|



biomes_f <- c("root:Environmental:Aquatic:Freshwater:Lotic:Low land river systems", "root:Environmental:Aquatic:Freshwater:Lentic")
set.seed(42)
sel_aquatic_freshwater <- st_data_raw |>
  filter(biome %in% biomes_f) |>
  group_by(biome, environment_biome, geo_loc_name_1) |>
  select(run_accession, biome, environment_biome, geo_loc_name_1) |>
  add_count() |>
  sample_n(1) |>
  arrange(desc(n)) |>
  ungroup() |>
  mutate(r = row_number()) |>
  select(biome, environment_biome, geo_loc_name_1) |>
  inner_join(st_data_raw) |>
  select(run_accession, biome, environment_biome, geo_loc_name_1) |>
  group_by(biome, environment_biome, geo_loc_name_1) |>
  sample_n(25, replace = TRUE) |>
  distinct() |>
  ungroup() |>
  select(run_accession, biome, environment_biome, geo_loc_name_1) |>
  mutate(biome_class = "root:Environmental:Aquatic:Freshwater")



##########################################################
## root:Host-associated:Mammals
##########################################################

# root:Host-associated:Animal:Digestive system:Fecal ----------------------
st_data_raw |>
  filter(biome == "root:Host-associated:Animal:Digestive system:Fecal") |>
  group_by(biome, environment_biome, geo_loc_name_1) |>
  select(run_accession, biome, environment_biome, geo_loc_name_1) |>
  add_count() |>
  sample_n(1) |>
  arrange(desc(n)) |>
  filter(n >= 25) |>
  knitr::kable()

# |run_accession |biome                                              |environment_biome |geo_loc_name_1 |  n|
# |:-------------|:--------------------------------------------------|:-----------------|:--------------|--:|
# |ERR2241685    |root:Host-associated:Animal:Digestive system:Fecal |pig               |Germany        | 41|
# |ERR2241945    |root:Host-associated:Animal:Digestive system:Fecal |poultry           |Germany        | 39|
# |ERR2245467    |root:Host-associated:Animal:Digestive system:Fecal |pig               |Belgium        | 33|
# |ERR2241833    |root:Host-associated:Animal:Digestive system:Fecal |pig               |Italy          | 31|
# |ERR2241646    |root:Host-associated:Animal:Digestive system:Fecal |pig               |Netherlands    | 31|
# |ERR2245483    |root:Host-associated:Animal:Digestive system:Fecal |pig               |Poland         | 31|

# root:Host-associated:Mammals:Gastrointestinal tract:Intestine:Fe --------
st_data_raw |>
  filter(biome == "root:Host-associated:Mammals:Gastrointestinal tract:Intestine:Fecal") |>
  group_by(biome, environment_biome, geo_loc_name_1) |>
  select(run_accession, biome, environment_biome, geo_loc_name_1) |>
  add_count() |>
  sample_n(1) |>
  arrange(desc(n)) |>
  filter(n >= 25) |>
  knitr::kable()

# pigs from farms
# |run_accession |biome                                                               |environment_biome |geo_loc_name_1 |   n|
# |:-------------|:-------------------------------------------------------------------|:-----------------|:--------------|---:|
# |ERR2597333    |root:Host-associated:Mammals:Gastrointestinal tract:Intestine:Fecal |                  |Denmark        | 218|

# root:Host-associated:Mammals:Digestive system:Fecal ---------------------
st_data_raw |>
  filter(biome == "root:Host-associated:Mammals:Digestive system:Fecal") |>
  group_by(biome, environment_biome, geo_loc_name_1) |>
  select(run_accession, biome, environment_biome, geo_loc_name_1) |>
  add_count() |>
  sample_n(1) |>
  arrange(desc(n)) |>
  knitr::kable()

# |run_accession |biome                                               |environment_biome |geo_loc_name_1 |  n|
# |:-------------|:---------------------------------------------------|:-----------------|:--------------|--:|
# |ERR1762069    |root:Host-associated:Mammals:Digestive system:Fecal |urban biome       |NA             | 72| # mice
# |ERR2059935    |root:Host-associated:Mammals:Digestive system:Fecal |                  |NA             | 34| # mice
# |ERR318681     |root:Host-associated:Mammals:Digestive system:Fecal |gi                |USA            | 31| # growing kitten
# |ERR1855536    |root:Host-associated:Mammals:Digestive system:Fecal |faeces            |Germany        | 22| # pig
# |ERR2011073    |root:Host-associated:Mammals:Digestive system:Fecal |terrestrial biome |USA            | 10| # mouse


# root:Host-associated:Mammals:Digestive system  --------------------------
st_data_raw |>
  filter(biome == "root:Host-associated:Mammals:Digestive system") |>
  group_by(biome, environment_biome, geo_loc_name_1) |>
  select(run_accession, biome, environment_biome, geo_loc_name_1) |>
  add_count() |>
  sample_n(1) |>
  arrange(desc(n)) |>
  knitr::kable()

# |run_accession |biome                                         |environment_biome           |geo_loc_name_1 |  n|
# |:-------------|:---------------------------------------------|:---------------------------|:--------------|--:|
# |SRR1747062    |root:Host-associated:Mammals:Digestive system |mammalia-associated habitat |Kenya          | 21| # baboon
# |SRR4116665    |root:Host-associated:Mammals:Digestive system |                            |NA             | 14| #mouse
# |SRR11852049   |root:Host-associated:Mammals:Digestive system |                            |Brazil: Tatui  | 12| # capibara

# root:Host-associated:Mammals:Digestive system:Large intestine:Fe --------
st_data_raw |>
  filter(biome == "root:Host-associated:Mammals:Digestive system:Large intestine:Fecal") |>
  group_by(biome, environment_biome, geo_loc_name_1) |>
  select(run_accession, biome, environment_biome, geo_loc_name_1) |>
  add_count() |>
  sample_n(1) |>
  arrange(desc(n)) |>
  knitr::kable()

# |run_accession |biome                                                               |environment_biome |geo_loc_name_1 |  n|
# |:-------------|:-------------------------------------------------------------------|:-----------------|:--------------|--:|
# |ERR1223845    |root:Host-associated:Mammals:Digestive system:Large intestine:Fecal |                  |NA             | 20| # swine
# |ERR1989796    |root:Host-associated:Mammals:Digestive system:Large intestine:Fecal |terrestrial biome |Switzerland    | 12| #mice

# root:Host-associated:Mammals:Digestive system:Stomach:Rumen -------------
st_data_raw |>
  filter(biome == "root:Host-associated:Mammals:Digestive system:Stomach:Rumen") |>
  group_by(biome, environment_biome, geo_loc_name_1) |>
  select(run_accession, biome, environment_biome, geo_loc_name_1) |>
  add_count() |>
  sample_n(1) |>
  arrange(desc(n)) |>
  knitr::kable()

# |run_accession |biome                                                       |environment_biome |geo_loc_name_1 |  n|
# |:-------------|:-----------------------------------------------------------|:-----------------|:--------------|--:|
# |ERR833213     |root:Host-associated:Mammals:Digestive system:Stomach:Rumen |                  |NA             | 20|

# root:Host-associated:Mammals:Digestive system:Foregut:Rumen -------------
st_data_raw |>
  filter(biome == "root:Host-associated:Mammals:Digestive system:Foregut:Rumen") |>
  group_by(biome, environment_biome, geo_loc_name_1) |>
  select(run_accession, biome, environment_biome, geo_loc_name_1) |>
  add_count() |>
  sample_n(1) |>
  arrange(desc(n)) |>
  knitr::kable()

# |run_accession |biome                                                       |environment_biome |geo_loc_name_1 |  n|
# |:-------------|:-----------------------------------------------------------|:-----------------|:--------------|--:|
# |SRR094926     |root:Host-associated:Mammals:Digestive system:Foregut:Rumen |                  |NA             |  8|
# |ERR1278105    |root:Host-associated:Mammals:Digestive system:Foregut:Rumen |moose rumen       |Sweden         |  3|


biomes_h_m_ds <- c(
  "root:Host-associated:Animal:Digestive system:Fecal",
  "root:Host-associated:Mammals:Digestive system:Stomach:Rumen",
  "root:Host-associated:Mammals:Digestive system:Large intestine:Fecal",
  "root:Host-associated:Mammals:Digestive system",
  "root:Host-associated:Mammals:Digestive system:Fecal",
  "root:Host-associated:Mammals:Gastrointestinal tract:Intestine:Fecal"
)

set.seed(42)
sel_host_associated_mammal_digestive <- st_data_raw |>
  filter(biome %in% biomes_h_m_ds) |>
  group_by(biome) |>
  select(run_accession, biome) |>
  add_count() |>
  sample_n(1) |>
  arrange(desc(n)) |>
  ungroup() |>
  mutate(r = row_number()) |>
  select(biome) |>
  inner_join(st_data_raw) |>
  select(run_accession, biome, environment_biome, geo_loc_name_1) |>
  group_by(biome) |>
  sample_n(25, replace = TRUE) |>
  distinct() |>
  ungroup() |>
  select(run_accession, biome, environment_biome, geo_loc_name_1) |>
  mutate(biome_class = "root:Host-associated:Mammals:Digestive")


# root:Host-associated:Plants  --------------------------------------------
st_data_raw |>
  filter(biome == "root:Host-associated:Plants") |>
  group_by(biome, environment_biome, geo_loc_name_1) |>
  select(run_accession, biome, environment_biome, geo_loc_name_1) |>
  add_count() |>
  sample_n(1) |>
  arrange(desc(n)) |>
  knitr::kable()

# |run_accession |biome                       |environment_biome   |geo_loc_name_1            |  n|
# |:-------------|:---------------------------|:-------------------|:-------------------------|--:|
# |ERR3929358    |root:Host-associated:Plants |plant               |Switzerland               | 18|
# |ERR3929370    |root:Host-associated:Plants |temperate grassland |Switzerland               | 18|
# |ERR3929376    |root:Host-associated:Plants |river               |Switzerland               |  8|
# |SRR6981896    |root:Host-associated:Plants |                    |USA: Michigan             |  5|
# |SRR4142418    |root:Host-associated:Plants |                    |Mexico:Guanajuato, Mexico |  2|
# |ERR3929364    |root:Host-associated:Plants |animal manure       |Switzerland               |  2|
# |ERR5526852    |root:Host-associated:Plants |sub-arctic tundra   |Iceland                   |  1|


# root:Host-associated:Plants:Rhizosphere ---------------------------------
st_data_raw |>
  filter(biome == "root:Host-associated:Plants:Rhizosphere") |>
  group_by(biome, environment_biome, geo_loc_name_1) |>
  select(run_accession, biome, environment_biome, geo_loc_name_1) |>
  add_count() |>
  sample_n(1) |>
  arrange(desc(n)) |>
  knitr::kable()

# |run_accession |biome                                   |environment_biome |geo_loc_name_1      |  n|
# |:-------------|:---------------------------------------|:-----------------|:-------------------|--:|
# |SRR6435963    |root:Host-associated:Plants:Rhizosphere |                  |USA: North Carolina | 19|
# |SRR6489884    |root:Host-associated:Plants:Rhizosphere |                  |USA: Nebraska       |  3|
# |SRR444039     |root:Host-associated:Plants:Rhizosphere |                  |NA                  |  1|



# root:Host-associated:Plants:Rhizosphere:Soil  ---------------------------
st_data_raw |>
  filter(biome == "root:Host-associated:Plants:Rhizosphere:Soil") |>
  group_by(biome, environment_biome, geo_loc_name_1) |>
  select(run_accession, biome, environment_biome, geo_loc_name_1) |>
  add_count() |>
  sample_n(1) |>
  arrange(desc(n)) |>
  knitr::kable()

# |run_accession |biome                                        |environment_biome |geo_loc_name_1 |  n|
# |:-------------|:--------------------------------------------|:-----------------|:--------------|--:|
# |SRR13617631   |root:Host-associated:Plants:Rhizosphere:Soil |                  |India:Assam    |  3|

biomes_h_p <- c(
  "root:Host-associated:Plants:Rhizosphere:Soil",
  "root:Host-associated:Plants:Rhizosphere",
  "root:Host-associated:Plants"
)
set.seed(42)
sel_host_associated_plants <- st_data_raw |>
  filter(biome %in% biomes_h_p) |>
  group_by(biome) |>
  select(run_accession, biome) |>
  add_count() |>
  sample_n(1) |>
  arrange(desc(n)) |>
  ungroup() |>
  mutate(r = row_number()) |>
  select(biome) |>
  inner_join(st_data_raw) |>
  select(run_accession, biome, environment_biome, geo_loc_name_1) |>
  distinct() |>
  ungroup() |>
  select(run_accession, biome) |>
  mutate(biome_class = "root:Host-associated:Plants")


##########################################################
## root:Environmental:Terrestrial:Soil
##########################################################

# root:Environmental:Terrestrial:Soil -------------------------------------
st_data_raw |>
  filter(biome == "root:Environmental:Terrestrial:Soil") |>
  group_by(biome, environment_biome, geo_loc_name_1) |>
  select(run_accession, biome, environment_biome, geo_loc_name_1) |>
  add_count() |>
  sample_n(1) |>
  arrange(desc(n)) |>
  knitr::kable()

# |run_accession |biome                               |environment_biome       |geo_loc_name_1            |  n|
# |:-------------|:-----------------------------------|:-----------------------|:-------------------------|--:|
# |ERR687883     |root:Environmental:Terrestrial:Soil |arid                    |NA                        | 24|
# |ERR687896     |root:Environmental:Terrestrial:Soil |temperate               |NA                        | 22|
# |ERR9752757    |root:Environmental:Terrestrial:Soil |farm soil               |USA                       | 12|
# |ERR1811650    |root:Environmental:Terrestrial:Soil |marine salt marsh biome |NA                        |  9|
# |ERR1606220    |root:Environmental:Terrestrial:Soil |field soil              |NA                        |  8|
# |SRR770300     |root:Environmental:Terrestrial:Soil |                        |NA                        |  5|
# |SRR5223441    |root:Environmental:Terrestrial:Soil |                        |Antarctica:Robinson Ridge |  2|
# |ERR5378028    |root:Environmental:Terrestrial:Soil |alpine tundra           |Iceland                   |  1|
# |SRR5223443    |root:Environmental:Terrestrial:Soil |desert                  |Antarctica:Robinson Ridge |  1|


set.seed(42)
sel_terrestrial_soil <- st_data_raw |>
  filter(biome == "root:Environmental:Terrestrial:Soil") |>
  filter(environment_biome != "marine salt marsh biome", environment_biome != "farm soil", environment_biome != "field soil") |>
  group_by(biome, environment_biome, geo_loc_name_1) |>
  select(run_accession, biome, environment_biome, geo_loc_name_1) |>
  add_count() |>
  sample_n(1) |>
  arrange(desc(n)) |>
  ungroup() |>
  mutate(r = row_number()) |>
  filter(r != 7) |>
  select(biome, environment_biome, geo_loc_name_1) |>
  inner_join(st_data_raw) |>
  select(run_accession, biome, environment_biome, geo_loc_name_1) |>
  group_by(biome, environment_biome, geo_loc_name_1) |>
  distinct() |>
  ungroup() |>
  select(run_accession, biome, environment_biome, geo_loc_name_1) |>
  mutate(biome_class = "root:Environmental:Terrestrial:Soil")




# root:Environmental:Terrestrial:Soil:Permafrost  -------------------------
st_data_raw |>
  filter(biome == "root:Environmental:Terrestrial:Soil:Permafrost") |>
  group_by(biome, environment_biome, geo_loc_name_1) |>
  select(run_accession, biome, environment_biome, geo_loc_name_1) |>
  add_count() |>
  sample_n(1) |>
  arrange(desc(n)) |>
  knitr::kable()

# |run_accession |biome                                          |environment_biome |geo_loc_name_1 |  n|
# |:-------------|:----------------------------------------------|:-----------------|:--------------|--:|
# |ERR3890099    |root:Environmental:Terrestrial:Soil:Permafrost |boreal peatlands  |Canada         | 21|
# |ERR1034454    |root:Environmental:Terrestrial:Soil:Permafrost |tundra            |USA            | 10|
# |ERR3890111    |root:Environmental:Terrestrial:Soil:Permafrost |microbial control |United Kingdom |  1|

set.seed(42)
sel_terrestrial_soil_permafrost <- st_data_raw |>
  filter(biome == "root:Environmental:Terrestrial:Soil:Permafrost") |>
  group_by(biome, environment_biome, geo_loc_name_1) |>
  select(run_accession, biome, environment_biome, geo_loc_name_1) |>
  add_count() |>
  sample_n(1) |>
  arrange(desc(n)) |>
  ungroup() |>
  mutate(r = row_number()) |>
  select(biome, environment_biome, geo_loc_name_1) |>
  inner_join(st_data_raw) |>
  select(run_accession, biome, environment_biome, geo_loc_name_1) |>
  group_by(biome, environment_biome, geo_loc_name_1) |>
  distinct() |>
  ungroup() |>
  select(run_accession, biome, environment_biome, geo_loc_name_1) |>
  mutate(biome_class = "root:Environmental:Terrestrial:Soil")



# root:Environmental:Terrestrial:Soil:Forest soil -------------------------
st_data_raw |>
  filter(biome == "root:Environmental:Terrestrial:Soil:Forest soil") |>
  group_by(biome, environment_biome, geo_loc_name_1) |>
  select(run_accession, biome, environment_biome, geo_loc_name_1) |>
  add_count() |>
  sample_n(1) |>
  arrange(desc(n)) |>
  knitr::kable()

# |run_accession |biome                                           |environment_biome |geo_loc_name_1            |  n|
# |:-------------|:-----------------------------------------------|:-----------------|:-------------------------|--:|
# |ERR753912     |root:Environmental:Terrestrial:Soil:Forest soil |forest            |NA                        | 21|
# |SRR11086695   |root:Environmental:Terrestrial:Soil:Forest soil |                  |Germany: Schwaebische Alb |  1|

set.seed(42)
sel_terrestrial_soil_forest <- st_data_raw |>
  filter(biome == "root:Environmental:Terrestrial:Soil:Forest soil") |>
  group_by(biome, environment_biome, geo_loc_name_1) |>
  select(run_accession, biome, environment_biome, geo_loc_name_1) |>
  add_count() |>
  sample_n(1) |>
  arrange(desc(n)) |>
  ungroup() |>
  mutate(r = row_number()) |>
  select(biome, environment_biome, geo_loc_name_1) |>
  inner_join(st_data_raw) |>
  select(run_accession, biome, environment_biome, geo_loc_name_1) |>
  group_by(biome, environment_biome, geo_loc_name_1) |>
  distinct() |>
  ungroup() |>
  select(run_accession, biome, environment_biome, geo_loc_name_1) |>
  mutate(biome_class = "root:Environmental:Terrestrial:Soil")

# root:Environmental:Terrestrial:Soil:Agricultural ------------------------
st_data_raw |>
  filter(biome == "root:Environmental:Terrestrial:Soil:Agricultural") |>
  group_by(biome, environment_biome, geo_loc_name_1) |>
  select(run_accession, biome, environment_biome, geo_loc_name_1) |>
  add_count() |>
  sample_n(1) |>
  arrange(desc(n)) |>
  knitr::kable()

# |run_accession |biome                                            |environment_biome |geo_loc_name_1 |  n|
# |:-------------|:------------------------------------------------|:-----------------|:--------------|--:|
# |ERR2486634    |root:Environmental:Terrestrial:Soil:Agricultural |cropland          |Finland        | 18|

set.seed(42)
sel_terrestrial_soil_agri <- st_data_raw |>
  filter(biome == "root:Environmental:Terrestrial:Soil:Agricultural" | environment_biome == "farm soil" | environment_biome == "field soil") |>
  group_by(biome, environment_biome, geo_loc_name_1) |>
  select(run_accession, biome, environment_biome, geo_loc_name_1) |>
  add_count() |>
  sample_n(1) |>
  arrange(desc(n)) |>
  ungroup() |>
  mutate(r = row_number()) |>
  select(biome, environment_biome, geo_loc_name_1) |>
  inner_join(st_data_raw) |>
  select(run_accession, biome, environment_biome, geo_loc_name_1) |>
  group_by(biome, environment_biome, geo_loc_name_1) |>
  distinct() |>
  ungroup() |>
  select(run_accession, biome, environment_biome, geo_loc_name_1) |>
  mutate(
    biome_class = "root:Environmental:Terrestrial:Soil",
    biome = "root:Environmental:Terrestrial:Soil:Agricultural"
  )

# root:Environmental:Terrestrial:Soil:Grasslands --------------------------
st_data_raw |>
  filter(biome == "root:Environmental:Terrestrial:Soil:Grasslands") |>
  group_by(biome, environment_biome, geo_loc_name_1) |>
  select(run_accession, biome, environment_biome, geo_loc_name_1) |>
  add_count() |>
  sample_n(1) |>
  arrange(desc(n)) |>
  knitr::kable()

# |run_accession |biome                                          |environment_biome   |geo_loc_name_1 |  n|
# |:-------------|:----------------------------------------------|:-------------------|:--------------|--:|
# |ERR1043165    |root:Environmental:Terrestrial:Soil:Grasslands |temperate grassland |USA            |  8|


set.seed(42)
sel_terrestrial_soil_grass <- st_data_raw |>
  filter(biome == "root:Environmental:Terrestrial:Soil:Grasslands") |>
  group_by(biome, environment_biome, geo_loc_name_1) |>
  select(run_accession, biome, environment_biome, geo_loc_name_1) |>
  add_count() |>
  sample_n(1) |>
  arrange(desc(n)) |>
  ungroup() |>
  mutate(r = row_number()) |>
  select(biome, environment_biome, geo_loc_name_1) |>
  inner_join(st_data_raw) |>
  select(run_accession, biome, environment_biome, geo_loc_name_1) |>
  group_by(biome, environment_biome, geo_loc_name_1) |>
  distinct() |>
  ungroup() |>
  select(run_accession, biome, environment_biome, geo_loc_name_1) |>
  mutate(biome_class = "root:Environmental:Terrestrial:Soil")


# root:Environmental:Terrestrial:Soil:Sand --------------------------------
st_data_raw |>
  filter(biome == "root:Environmental:Terrestrial:Soil:Sand") |>
  group_by(biome, environment_biome, geo_loc_name_1) |>
  select(run_accession, biome, environment_biome, geo_loc_name_1) |>
  add_count() |>
  sample_n(1) |>
  arrange(desc(n)) |>
  knitr::kable()


# |run_accession |biome                                    |environment_biome |geo_loc_name_1                     |  n|
# |:-------------|:----------------------------------------|:-----------------|:----------------------------------|--:|
# |SRR1569898    |root:Environmental:Terrestrial:Soil:Sand |                  |USA: Municipal Pensacola Beach, FL | 12|


set.seed(42)
sel_terrestrial_soil_sand <- st_data_raw |>
  filter(biome == "root:Environmental:Terrestrial:Soil:Sand") |>
  group_by(biome, environment_biome, geo_loc_name_1) |>
  select(run_accession, biome, environment_biome, geo_loc_name_1) |>
  add_count() |>
  sample_n(1) |>
  arrange(desc(n)) |>
  ungroup() |>
  mutate(r = row_number()) |>
  select(biome, environment_biome, geo_loc_name_1) |>
  inner_join(st_data_raw) |>
  select(run_accession, biome, environment_biome, geo_loc_name_1) |>
  group_by(biome, environment_biome, geo_loc_name_1) |>
  distinct() |>
  ungroup() |>
  select(run_accession, biome, environment_biome, geo_loc_name_1) |>
  mutate(biome_class = "root:Environmental:Terrestrial:Soil")


###########################################################################
# Combine all selections -------------------------------------------------#
###########################################################################

all_selection <- sel_aquatic_estuary |>
  bind_rows(sel_aquatic_freshwater) |>
  bind_rows(sel_aquatic_freshwater_lake) |>
  bind_rows(sel_aquatic_freshwater_sediment) |>
  bind_rows(sel_aquatic_marine) |>
  bind_rows(sel_aquatic_marine_coastal) |>
  bind_rows(sel_aquatic_marine_intertidal_salt_marsh) |>
  bind_rows(sel_aquatic_marine_oceanic) |>
  bind_rows(sel_aquatic_marine_sediment) |>
  bind_rows(sel_host_associated_mammal_digestive) |>
  bind_rows(sel_host_associated_plants) |>
  bind_rows(sel_terrestrial_soil) |>
  bind_rows(sel_terrestrial_soil_agri) |>
  bind_rows(sel_terrestrial_soil_forest) |>
  bind_rows(sel_terrestrial_soil_grass) |>
  bind_rows(sel_terrestrial_soil_permafrost) |>
  bind_rows(sel_terrestrial_soil_sand)


all_selection |>
  group_by(biome_class, biome, environment_biome, geo_loc_name_1) |>
  count(sort = TRUE) |>
  knitr::kable()

selected_biomes <- all_selection |>
  mutate(biome_subclass = case_when(
    biome == "root:Host-associated:Plants" ~ "root:Host-associated:Plants",
    biome == "root:Environmental:Aquatic:Estuary" ~ "root:Environmental:Aquatic:Estuary",
    biome == "root:Environmental:Aquatic:Freshwater:Lake" & environment_biome == "atlantic_highlands" ~ "root:Environmental:Aquatic:Freshwater:Lake:Atlantic highlands",
    biome == "root:Environmental:Aquatic:Freshwater:Lake" & environment_biome == "atlantic_maritime" ~ "root:Environmental:Aquatic:Freshwater:Lake:Atlantic maritime",
    biome == "root:Environmental:Aquatic:Freshwater:Lake" & environment_biome == "boreal_plains" ~ "root:Environmental:Aquatic:Freshwater:Lake:Boreal plains",
    biome == "root:Environmental:Aquatic:Freshwater:Lake" & environment_biome == "boreal_shield" ~ "root:Environmental:Aquatic:Freshwater:Lake:Boreal shield",
    biome == "root:Environmental:Aquatic:Freshwater:Lake" & environment_biome == "mixedwood_plains" ~ "root:Environmental:Aquatic:Freshwater:Lake:Mixewood plains",
    biome == "root:Environmental:Aquatic:Freshwater:Lake" & environment_biome == "montane cordillera" ~ "root:Environmental:Aquatic:Freshwater:Lake:Montane cordillera",
    biome == "root:Environmental:Aquatic:Freshwater:Lake" & environment_biome == "pacific maritime" ~ "root:Environmental:Aquatic:Freshwater:Lake:Pacific maritime",
    biome == "root:Environmental:Aquatic:Freshwater:Lake" & environment_biome == "praires" ~ "root:Environmental:Aquatic:Freshwater:Lake:Prairies",
    biome == "root:Environmental:Aquatic:Freshwater:Lake" & environment_biome == "semi_arid_plateau" ~ "root:Environmental:Aquatic:Freshwater:Lake:Semi arid plateau",
    biome == "root:Environmental:Aquatic:Marine" & geo_loc_name_1 == "Canada:Atlantic" ~ "root:Environmental:Aquatic:Marine:Atlantic ocean",
    biome == "root:Environmental:Aquatic:Marine" & geo_loc_name_1 == "Indian Ocean" ~ "root:Environmental:Aquatic:Marine:Indian ocean",
    biome == "root:Environmental:Aquatic:Marine" & geo_loc_name_1 == "Baltic Sea" ~ "root:Environmental:Aquatic:Marine:brackis water",
    biome == "root:Environmental:Aquatic:Marine" & geo_loc_name_1 == "Atlantic Ocean" ~ "root:Environmental:Aquatic:Marine:Atlantic ocean",
    biome == "root:Environmental:Aquatic:Marine" & geo_loc_name_1 == "Pacific Ocean" ~ "root:Environmental:Aquatic:Marine:Pacific ocean",
    biome == "root:Environmental:Aquatic:Marine" & environment_biome == "polar biome" ~ "root:Environmental:Aquatic:Marine:Polar",
    biome == "root:Environmental:Aquatic:Marine:Oceanic" & geo_loc_name_1 == "Atlantic Ocean: Sargasso Sea, BATS" ~ "root:Environmental:Aquatic:Marine:Sargasso Sea",
    biome == "root:Environmental:Aquatic:Marine:Oceanic" & geo_loc_name_1 == "Pacific Ocean: North Pacific Subtropical Gyre, Station ALOHA" ~ "root:Environmental:Aquatic:Marine:North Pacific Subtropical Gyre",
    biome == "root:Environmental:Terrestrial:Soil" & environment_biome == "arid" ~ "root:Environmental:Terrestrial:Soil:Arid",
    biome_class == "root:Host-associated:Mammals:Digestive" ~ "root:Host-associated:Mammals:Digestive",
    biome_class == "root:Host-associated:Plants" ~ "root:Host-associated:Plants",
    biome == "root:Environmental:Terrestrial:Soil" & environment_biome == "temperate" ~ "root:Environmental:Terrestrial:Soil:Temperate",
    biome == "root:Environmental:Terrestrial:Soil:Forest soil" ~ "root:Environmental:Terrestrial:Soil:Forest",
    biome == "root:Environmental:Terrestrial:Soil:Permafrost" ~ "root:Environmental:Terrestrial:Soil:Permafrost",
    biome == "root:Environmental:Terrestrial:Soil:Agricultural" ~ "root:Environmental:Terrestrial:Soil:Agricultural",
    biome == "root:Environmental:Aquatic:Marine:Coastal" ~ "root:Environmental:Aquatic:Marine:Coastal",
    biome == "root:Environmental:Aquatic:Marine:Sediment" ~ "root:Environmental:Aquatic:Marine:Sediment",
    biome == "root:Environmental:Terrestrial:Soil:Sand" ~ "root:Environmental:Terrestrial:Soil:Sand",
    biome == "root:Environmental:Aquatic:Freshwater:Lentic:Sediment" ~ "root:Environmental:Aquatic:Freshwater:Sediment",
    biome == "root:Environmental:Aquatic:Freshwater:Sediment" ~ "root:Environmental:Aquatic:Freshwater:Sediment",
    biome == "root:Environmental:Aquatic:Freshwater:Lotic:Sediment" ~ "root:Environmental:Aquatic:Freshwater:Sediment",
    biome == "root:Environmental:Terrestrial:Soil:Grasslands" ~ "root:Environmental:Terrestrial:Soil:Grasslands",
    biome == "root:Environmental:Aquatic:Marine:Intertidal zone:Salt marsh" ~ "root:Environmental:Aquatic:Marine:Salt marsh",
    biome == "root:Environmental:Aquatic:Marine:Salt marsh" ~ "root:Environmental:Aquatic:Marine:Salt marsh",
    biome == "root:Environmental:Aquatic:Freshwater:Lotic:Low land river systems" ~ "root:Environmental:Aquatic:Freshwater",
    biome == "root:Environmental:Aquatic:Freshwater:Lentic" ~ "root:Environmental:Aquatic:Freshwater",
    biome == "root:Environmental:Terrestrial:Soil" & geo_loc_name_1 == "Antarctica:Robinson Ridge" ~ "root:Environmental:Terrestrial:Soil:Antarctica",
    biome == "root:Environmental:Terrestrial:Soil" & environment_biome == "" ~ "root:Environmental:Terrestrial:Soil",
  ))

selected_biomes |>
  group_by(biome_class, biome_subclass) |>
  count(sort = TRUE) |>
  filter(!is.na(biome_subclass)) |>
  knitr::kable()


# |biome_class                            |biome_subclass                                                   |   n|
# |:--------------------------------------|:----------------------------------------------------------------|---:|
# |root:Host-associated:Mammals:Digestive |root:Host-associated:Mammals:Digestive                           | 123|
# |root:Environmental:Aquatic:Marine      |root:Environmental:Aquatic:Estuary                               |  97|
# |root:Host-associated:Plants            |root:Host-associated:Plants                                      |  80|
# |root:Environmental:Aquatic:Marine      |root:Environmental:Aquatic:Marine:Atlantic ocean                 |  50|
# |root:Environmental:Aquatic:Marine      |root:Environmental:Aquatic:Marine:Sediment                       |  44|
# |root:Environmental:Aquatic:Freshwater  |root:Environmental:Aquatic:Freshwater:Sediment                   |  38|
# |root:Environmental:Terrestrial:Soil    |root:Environmental:Terrestrial:Soil:Agricultural                 |  38|
# |root:Environmental:Aquatic:Marine      |root:Environmental:Aquatic:Marine:Coastal                        |  37|
# |root:Environmental:Terrestrial:Soil    |root:Environmental:Terrestrial:Soil:Permafrost                   |  32|
# |root:Environmental:Aquatic:Marine      |root:Environmental:Aquatic:Marine:Salt marsh                     |  29|
# |root:Environmental:Aquatic:Freshwater  |root:Environmental:Aquatic:Freshwater:Lake:Atlantic highlands    |  25|
# |root:Environmental:Aquatic:Freshwater  |root:Environmental:Aquatic:Freshwater:Lake:Atlantic maritime     |  25|
# |root:Environmental:Aquatic:Freshwater  |root:Environmental:Aquatic:Freshwater:Lake:Boreal plains         |  25|
# |root:Environmental:Aquatic:Freshwater  |root:Environmental:Aquatic:Freshwater:Lake:Boreal shield         |  25|
# |root:Environmental:Aquatic:Freshwater  |root:Environmental:Aquatic:Freshwater:Lake:Mixewood plains       |  25|
# |root:Environmental:Aquatic:Freshwater  |root:Environmental:Aquatic:Freshwater:Lake:Montane cordillera    |  25|
# |root:Environmental:Aquatic:Freshwater  |root:Environmental:Aquatic:Freshwater:Lake:Pacific maritime      |  25|
# |root:Environmental:Aquatic:Freshwater  |root:Environmental:Aquatic:Freshwater:Lake:Prairies              |  25|
# |root:Environmental:Aquatic:Freshwater  |root:Environmental:Aquatic:Freshwater:Lake:Semi arid plateau     |  25|
# |root:Environmental:Aquatic:Marine      |root:Environmental:Aquatic:Marine:brackis water                  |  25|
# |root:Environmental:Aquatic:Marine      |root:Environmental:Aquatic:Marine:Indian ocean                   |  25|
# |root:Environmental:Aquatic:Marine      |root:Environmental:Aquatic:Marine:North Pacific Subtropical Gyre |  25|
# |root:Environmental:Aquatic:Marine      |root:Environmental:Aquatic:Marine:Pacific ocean                  |  25|
# |root:Environmental:Aquatic:Marine      |root:Environmental:Aquatic:Marine:Polar                          |  25|
# |root:Environmental:Aquatic:Marine      |root:Environmental:Aquatic:Marine:Sargasso Sea                   |  25|
# |root:Environmental:Terrestrial:Soil    |root:Environmental:Terrestrial:Soil:Arid                         |  24|
# |root:Environmental:Terrestrial:Soil    |root:Environmental:Terrestrial:Soil:Forest                       |  22|
# |root:Environmental:Terrestrial:Soil    |root:Environmental:Terrestrial:Soil:Temperate                    |  22|
# |root:Environmental:Aquatic:Freshwater  |root:Environmental:Aquatic:Freshwater                            |  13|
# |root:Environmental:Terrestrial:Soil    |root:Environmental:Terrestrial:Soil:Sand                         |  12|
# |root:Environmental:Terrestrial:Soil    |root:Environmental:Terrestrial:Soil:Grasslands                   |   8|
# |root:Environmental:Terrestrial:Soil    |root:Environmental:Terrestrial:Soil                              |   5|
# |root:Environmental:Terrestrial:Soil    |root:Environmental:Terrestrial:Soil:Antarctica                   |   3|

selected_biomes |>
  select(run_accession, biome_class, biome_subclass, biome) |>
  inner_join(st_data_raw |> select(-biome)) |>
  write_tsv("./data/sourcetracker/cdata/kapk-biomes-download.txt")
