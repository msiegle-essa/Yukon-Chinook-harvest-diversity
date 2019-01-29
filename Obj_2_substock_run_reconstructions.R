

# ----------------------------------------------------------------- #
# ----------------------------------------------------------------- #
# 2. Sub-stock run reconstructions
#
# code: Matt Siegle, Brendan Connors
# last updated: 25 Jan 2019
# ----------------------------------------------------------------- #
# ----------------------------------------------------------------- #

# create dataframe to get sub-stock specific brood table ####
# set working directory
setwd("~/Dropbox (ESSA Technologies)/EN2438 - Yukon Chinook/Data and analyses/GIT_R.file_organization/input files")

# data input for sub-stock run reconstructions
df <- read.csv(file.choose()) # input_yukon_chin_bt_red.csv
head(df)

require(tidyr)
reg_prop <- gather(df, -year, -sample_size, -run, -spwn, -per_fem,
                   -comm_har, -dom_har, -fsc_har, -rec_har, -yk_rv_har,
                   -por_har, -total_har, -er, -age_4,-age_5, -age_6,
                   -age_7, key = "region", value = "prop") %>%
  arrange(year, region) %>%
  mutate(run_stock = run*(prop/100),
         spwn_stock = spwn*(prop/100),
         tot_har_stock = total_har*(prop/100))


# create data for Objective 3 State Space models ####
age_data <- reg_prop %>%
  select(year, region, age_4, age_5, 
         age_6, age_7) %>%
  mutate(reg_code = case_when(region == "Yukon.Carmacks" ~ "1",
                              region == "Yukon.Lower.Canadian" ~ "2",
                              region == "Yukon.mainstem" ~ "3",
                              region == "Yukon.Pelly" ~ "4",
                              region == "Yukon.Stewart" ~ "5",
                              region == "Yukon.Teslin" ~ "6",
                              region == "Yukon.upper" ~ "7",
                              region == "Yukon.White.Donjek" ~ "8")) %>%
  mutate(gooddata = case_when(year == 1984 ~ 0,
                              year == 1988 ~ 0,
                              year == 1989 ~ 0,
                              year == 1990 ~ 0,
                              year == 1998 ~ 0,
                              year == 2006 ~ 0,
                              year == 2007 ~ 0,
                              year < 1984 ~ 1,
                              year >= 1985 & year <= 1987 ~ 1,
                              year >= 1991 & year <= 1997 ~1,
                              year >= 1999 & year <= 2005 ~1,
                              year >= 2008 ~1))

esc_data <- reg_prop %>%
  select(year, region, spwn_stock) %>%
  mutate(reg_code = case_when(region == "Yukon.Carmacks" ~ "1",
                              region == "Yukon.Lower.Canadian" ~ "2",
                              region == "Yukon.mainstem" ~ "3",
                              region == "Yukon.Pelly" ~ "4",
                              region == "Yukon.Stewart" ~ "5",
                              region == "Yukon.Teslin" ~ "6",
                              region == "Yukon.upper" ~ "7",
                              region == "Yukon.White.Donjek" ~ "8")) %>%
  mutate(truecount = case_when(year == 1984 ~ 0,
                               year == 1988 ~ 0,
                               year == 1989 ~ 0,
                               year == 1990 ~ 0,
                               year == 1998 ~ 0,
                               year == 2006 ~ 0,
                               year == 2007 ~ 0,
                               year < 1984 ~ 1,
                               year >= 1985 & year <= 1987 ~ 1,
                               year >= 1991 & year <= 1997 ~1,
                               year >= 1999 & year <= 2005 ~1,
                               year >= 2008 ~1))

harv_data <- reg_prop %>%
  select(year, region, tot_har_stock) %>%
  mutate(reg_code = case_when(region == "Yukon.Carmacks" ~ "1",
                              region == "Yukon.Lower.Canadian" ~ "2",
                              region == "Yukon.mainstem" ~ "3",
                              region == "Yukon.Pelly" ~ "4",
                              region == "Yukon.Stewart" ~ "5",
                              region == "Yukon.Teslin" ~ "6",
                              region == "Yukon.upper" ~ "7",
                              region == "Yukon.White.Donjek" ~ "8")) %>%
  mutate(truecount = case_when(year == 1984 ~ 0,
                               year == 1988 ~ 0,
                               year == 1989 ~ 0,
                               year == 1990 ~ 0,
                               year == 1998 ~ 0,
                               year == 2006 ~ 0,
                               year == 2007 ~ 0,
                               year < 1984 ~ 1,
                               year >= 1985 & year <= 1987 ~ 1,
                               year >= 1991 & year <= 1997 ~1,
                               year >= 1999 & year <= 2005 ~1,
                               year >= 2008 ~1))


# write csv files for Objective 3 data input
write.table(age_data, file="age_data.csv", sep=" ")
write.table(esc_data, file="esc_data.csv", sep=" ")
write.table(harv_data, file="harv_data.csv", sep=" ")


# export exploitation rates for Objective 3 data input
er <- df %>%
  select(year, er)
write.table(er, file="er.csv", sep=" ")





