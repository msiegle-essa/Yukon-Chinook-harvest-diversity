

# ----------------------------------------------------------------- #
# ----------------------------------------------------------------- #
# 1. Explore Substock Composition
#
# code: Matt Siegle and Brendan Connors
# last updated: 29 Jan 2019
# ----------------------------------------------------------------- #
# ----------------------------------------------------------------- #

# Load necessary libraries and data ####
library(plyr)
library(lme4)
library(ggplot2)


# set working directory ####
setwd("~/Dropbox (ESSA Technologies)/EN2438 - Yukon Chinook/Yukon Chinook harvest-diversity/input files")

#read in border counts (fish wheel and sonar)
border_counts <- read.delim(file="input_bordercounts2.txt",header=TRUE)
head(border_counts)

detach("package:plyr")
require(dplyr)

# sum border counts by year
View(border_counts %>%
       na.omit() %>%
       group_by(year) %>%
       summarise(sum = sum(count)))

# add column of counts per year to new dataframe bc2
bc2 <- border_counts %>% 
  mutate(sum = case_when(year == "1985" ~ 1321,
                         year == "1986" ~ 1998,
                         year == "1987" ~ 938,
                         year == "1988" ~ 976,
                         year == "1989" ~ 1065,
                         year == "1990" ~ 1361,
                         year == "1991" ~ 1726,
                         year == "1992" ~ 1889,
                         year == "1993" ~ 1241,
                         year == "1994" ~ 1290,
                         year == "1995" ~ 2215,
                         year == "1996" ~ 1749,
                         year == "1997" ~ 2221,
                         year == "1998" ~ 1080,
                         year == "1999" ~ 914,
                         year == "2000" ~ 1494,
                         year == "2001" ~ 3986,
                         year == "2002" ~ 1065,
                         year == "2003" ~ 1276,
                         year == "2004" ~ 1361,
                         year == "2005" ~ 81529,
                         year == "2006" ~ 73691,
                         year == "2007" ~ 41697,
                         year == "2008" ~ 38097,
                         year == "2009" ~ 69963,
                         year == "2010" ~ 35074,
                         year == "2011" ~ 51271,
                         year == "2012" ~ 34747,
                         year == "2013" ~ 30725,
                         year == "2014" ~ 63462,
                         year == "2015" ~ 84015,
                         year == "2016" ~ 72329))

# create proportion of total run count per julian day count
bc2$prop <- NA
bc2$prop <- bc2$count/bc2$sum


#read in Genetic Stock ID data from fish collected at border 
# ("probability" is the probability individual fish originated from stock X)
GSI <- read.csv(file="input_Yukon_GSI_82_05_longform2.csv",header=TRUE)
head(GSI)

# remove assignments within each sample with p < 0.5 (effectively generates data frame with single assignment per sample)
detach("package:dplyr")
require(plyr)

GSI.region <- subset(GSI, prob > 0.5)
head(GSI.region.2)

# create data frame of number of samples per day and year 
GSI.counts <- ddply(GSI.region,c("year","julian_date"),function(x){
  count <- dim(x)[1]
  data.frame(count)
})

# add column of gsi counts per year to new dataframe: only keep samples with prob > 0.5
detach("package:plyr")
require(dplyr)

gsi2 <- GSI.counts %>% 
  mutate(year_count = case_when(year == "1982" ~ 124,
                                year == "1983" ~ 141,
                                year == "1985" ~ 149,
                                year == "1986" ~ 149,
                                year == "1987" ~ 148,
                                year == "1991" ~ 147,
                                year == "1992" ~ 149,
                                year == "1993" ~ 149,
                                year == "1994" ~ 149,
                                year == "1995" ~ 149,
                                year == "1996" ~ 144,
                                year == "1997" ~ 150,
                                year == "1998" ~ 0,
                                year == "1999" ~ 147,
                                year == "2000" ~ 148,
                                year == "2001" ~ 149,
                                year == "2002" ~ 150,
                                year == "2003" ~ 148,
                                year == "2004" ~ 131,
                                year == "2005" ~ 142,
                                year == "2006" ~ 150,
                                year == "2007" ~ 149,
                                year == "2008" ~ 452,
                                year == "2009" ~ 646,
                                year == "2010" ~ 467,
                                year == "2011" ~ 497,
                                year == "2012" ~ 344,
                                year == "2013" ~ 290,
                                year == "2014" ~ 708,
                                year == "2015" ~ 1026,
                                year == "2016" ~ 728))

# create proportion of total run count per julian day count
gsi2$prop <- NA
gsi2$prop <- gsi2$count/gsi2$year_count

# create sample size number dataframe for geom_text in plot
dat_text <- data.frame(label = c("124", "141", "149","149", "148", "147",
                                 "149", "149", "149","149", "144", "150", "147", "148",
                                 "149", "150", "148","131", "142","150","149","452","646",
                                 "467","497","344","290","708","1026","728"),
                       year = c(1982,1983,1985,1986,1987,1991,
                                1992,1993,1994,1995,1996,1997, 1999,2000,
                                2001,2002,2003,2004,2005,2006, 2007,2008,2009,
                                2010,2011,2012,2013,2014,2015,2016))

# Run size distribution and GSI size and distribution, memo Fig 2 ####
n <- c(1982:1987,1991:1997,1999:2016)
bc3 <- bc2 %>%
  filter(year %in% n)

gsi3 <- gsi2 %>%
  filter(year %in% n)


# FIGURE: border counts and gsi sample size; memo figure 2 ####
ggplot(bc3, aes(x=julian, y = prop)) +
  geom_bar(stat = "identity") +
  geom_bar(data = gsi3 %>% filter(julian_date < 1000), aes(x=julian_date, y=prop), 
           fill="red", alpha=0.5, stat = "identity") +
  scale_x_continuous(limits=c(150,300), breaks = c(150,175,200,225,250,275,300)) +
  xlab("Julian day") +
  ylab("Run size and GSI sample size") +
  facet_wrap(~year, scales = "free_y") +
  theme_bw() +
  theme(axis.text.x = element_text(size=9, angle = 45, hjust = 1)) +
  theme(strip.text = element_text(size=10)) +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  geom_text(data = dat_text, mapping = aes(x = -Inf, y = -Inf, label = label),
            hjust = -0.3, vjust = -3, size=2.5)


# Figure 2y: 2008/2010/2011/2012 fish wheel test fish run timing compare
bc4 <- bc3 %>%
  filter(year == 2008 | year == 2010 | year == 2011 | year == 2012)
  #filter(!is.na(prop))
  
g <- read.csv(file.choose()) %>%    # in compare FW and test folder, input_wheel_test_jul_comp.csv
  filter(year == 2008 | year == 2011 | year == 2010 | year ==2012) %>%
  filter(prop > 0.5)


detach("package:dplyr")
require(plyr)

gc <- ddply(g,c("year","julian","gear"),function(x){
  count <- dim(x)[1]
  data.frame(count)
})

detach("package:plyr")
require(dplyr)

gc_sum <- gc %>%
  group_by(year,gear) %>%
  summarise(sum = sum(count))  


gjulc <- gc %>%
  group_by(year, gear, julian) %>%
  mutate(sum = case_when(year == 2008 & gear == "Fish Wheel" ~ 139,
                         year == 2008 & gear == "Test Fishery" ~ 324,
                         year == 2011 & gear == "Fish Wheel" ~ 126,
                         year == 2011 & gear == "Test Fishery" ~ 309,
                         year == 2010 & gear == "Fish Wheel" ~ 129,
                         year == 2010 & gear == "Test Fishery" ~ 376,
                         year == 2012 & gear == "Fish Wheel" ~ 140,
                         year == 2012 & gear == "Test Fishery" ~ 273)) %>%
  mutate(freq = count/sum)


# plot figure
ggplot(bc4, aes(x=julian, y = prop*100)) +
  geom_bar(stat = "identity") +
  geom_bar(data = gjulc, aes(x=julian, y=freq*100), 
           fill="red", alpha=0.5, stat = "identity") +
  scale_x_continuous(limits=c(180,240), breaks = c(180,200,220,240)) +
  xlab("Julian day") +
  ylab("Run size and GSI sample size") +
  facet_grid(gear~year, scales = "free_y") +
  theme_bw() +
  theme(axis.ticks.y = element_blank(),
        axis.text.y = element_blank())
ggsave(file = "julian_TFvFW.pdf", width = 7, height = 5)



# Sub-stock proportion test fishery catch by facets sub stock ####
# data input
setwd("~/Dropbox (ESSA Technologies)/EN2438 - Yukon Chinook/Data and analyses/GIT_R.file_organization/input files")
d <- read.csv("input_gsi_error.csv") # input_gsi_error.csv
head(d)
dd <- spread(d, value_level, value)

# FIGURE sub-stock proportion test fishery catch by facets sub stock; memo Figure 3
ggplot(dd, aes(x = year, y = estimate)) +  #%>% filter(region_num =="1")
  geom_bar(stat = "identity", colour = "black", fill="light grey") +
  geom_errorbar(aes(ymin=lower, ymax=upper), width=.2) +
  geom_vline(xintercept = 2008, lty = "dotted") +
  xlab("Year") +
  ylab("Proportion of test fishery catch (%)") +
  scale_x_continuous(breaks = c(1984,1988,1992,1996,2000,2004,2008,2012,2016)) +
  #scale_y_continuous(limits=c(0,40), breaks=c(0,10,20,30,40)) +
  theme_bw() +
  facet_wrap(~Region.Name, ncol=3, scales = "free_y") +
  theme(strip.text = element_text(size=10)) +
  theme(axis.text.x = element_text(angle=45, hjust=1)) +
  theme(axis.text = element_text(size=12)) +
  theme(axis.title = element_text(size=12)) +
  theme(legend.justification = "center")
ggsave(file=" .pdf", width=7,height=5)


# SAMPLING ERROR MEMO ####

# set working directory
setwd("~/Dropbox (ESSA Technologies)/EN2438 - Yukon Chinook/Data and analyses/GIT_R.file_organization/input files")

# sample size effect ####
# input data file and dataframe manipulation
d <- read.csv(file.choose()) # input_2012-2016_fullsampsizetet.csv

dd <- d %>% filter(year >= 2014) # only keep years 2014-2016

x <- as.factor(c(15:350))
df <- dd %>%
  #group_by(year, sample_num, stock, region) %>%
  filter(prob >= 0.51) %>%
  filter(julian_date %in% x) %>%
  filter(julian_date != "unknwn") %>%
  mutate(substock = case_when(region == "30" ~ "Upper Lakes and Mainstem",
                              region == "31" ~ "Teslin River",
                              region == "32" ~ "Carmacks",
                              region == "33" ~ "Middle Mainstem",
                              region == "34" ~ "Pelly",
                              region == "35" ~ "Stewart",
                              region == "36" ~ "Lower Mainstem",
                              region == "38" ~ "White-Donjek"))

# count total samples per year, only keep samples with prob > 0.5
df %>% 
  group_by(year) %>%
  summarise(count = n())

# create sample size bins dataframes

# 150 sample size bin
datalist_150 = list()
for (i in 1:1000) {
  # ... make some data
  dff <- df %>%
    group_by(year) %>%
    sample_n(size=150,replace=FALSE)
  dff$bin <- 150
  dff$i <- i  # maybe you want to keep track of which iteration produced it?
  datalist_150[[i]] <- dff # add it to your list
}
bg_150 <- dplyr::bind_rows(datalist_150)

# 250 sample size bin
datalist_250 = list()
for (i in 1:1000) {
  # ... make some data
  dff <- df %>%
    group_by(year) %>%
    sample_n(size=250,replace=FALSE)
  dff$bin <- 250
  dff$i <- i  # maybe you want to keep track of which iteration produced it?
  datalist_250[[i]] <- dff # add it to your list
}
bg_250 <- dplyr::bind_rows(datalist_250)

# 350 sample size bin
datalist_350 = list()
for (i in 1:1000) {
  # ... make some data
  dff <- df %>%
    group_by(year) %>%
    sample_n(size=350,replace=FALSE)
  dff$bin <- 350
  dff$i <- i  # maybe you want to keep track of which iteration produced it?
  datalist_350[[i]] <- dff # add it to your list
}
bg_350 <- dplyr::bind_rows(datalist_350)


# 450 sample size bin
datalist_450 = list()
for (i in 1:1000) {
  # ... make some data
  dff <- df %>%
    group_by(year) %>%
    sample_n(size=450,replace=FALSE)
  dff$bin <- 450
  dff$i <- i  # maybe you want to keep track of which iteration produced it?
  datalist_450[[i]] <- dff # add it to your list
}
bg_450 <- dplyr::bind_rows(datalist_450)

# join all four data frames together
dd1 <- do.call("rbind", list(bg_150, bg_250, bg_350, bg_450))

# create new dataframe that computes substock proportion for each year, bin, and i
dd2 <- dd1 %>%
  group_by(year, region, bin, i) %>%
  mutate(count = n()) %>%
  mutate(freq = count/bin)

# only keep distinct region proportions per year per sample size bin and per replicate i
dd3 <- dd2 %>%
  distinct(year, substock, bin, i, freq)

dd3$bin <- as.factor(dd3$bin)
dd3$year <- as.factor(dd3$year)

head(dd3)

# plot boxplot with different sample size bins
ggplot(dd3 %>% filter(bin != "50") %>% 
         filter(bin != "450") %>% filter(bin != "350"), 
       aes(x=year, y=freq*100, fill=bin)) +
  geom_boxplot() + # outlier.shape = NA
  xlab("Year") +
  ylab("Substock proportion") +
  facet_wrap(~substock, nrow=3, scales = "free_y") +
  theme_bw()
ggsave(file="samplesizecheck.pdf", height=5, width=7)



# create coef of variation dataframe
ddcv <- dd3 %>%
  distinct(year, substock, bin, i, freq)

ddcv2 <- ddcv %>%
  group_by(year, substock, bin) %>%
  mutate(sd = sd(freq) + 0.00001) %>%
  mutate(mean = mean(freq)) %>%
  mutate(cv = sd/mean) %>%
  arrange(year, substock)

ddcv2$bin <- as.factor(ddcv2$bin)


# FIGURE plot CV for year/substock by sample size bin
names <- factor(ddcv2$substock, levels = c("Lower Mainstem", "White-Donjek", "Upper Lakes and Mainstem", 
                                           "Stewart", "Carmacks", "Middle Mainstem", "Pelly", "Teslin River"))

ggplot(ddcv2, aes(x=names, y=cv, fill=bin)) +
  geom_boxplot() +
  xlab("Sub-stock") +
  ylab("Coefficient of Variation") +
  theme_bw() +
  theme(axis.title = element_text(size=12)) +
  theme(axis.text = element_text(size = 11)) +
  theme(axis.text.x = element_text(angle=30, hjust=1))
ggsave(file="CV_reg.comb.pdf", height=5, width=7)




# CV summary info
ddcv3 <- ddcv2 %>%
  select(-i, -region, -freq, -sd, -mean) %>%
  distinct() %>%
  arrange(year, substock, bin) %>%
  spread(key = bin, value = cv) %>%
  group_by(year, substock) %>%
  mutate(bin_250 = ((`150` - `250`)/`150`)*-100) %>%
  mutate(bin_350 = ((`250` - `350`)/`250`)*-100) %>%
  mutate(bin_450 = ((`350` - `450`)/`350`)*-100) %>%
  ungroup() %>%
  gather(key = bin, value = CV_reduction, bin_250:bin_450) %>%
  gather(key = bin2, value = CV, `150`:`450`)
  

# plot CV reduction
ggplot(ddcv3, aes(x = bin, y = CV_reduction, color = year, group=year)) +
  geom_point(size = 3) +
  geom_line() +
  xlab("Sample size bin") +
  ylab("Percent reduction in CV") +
  scale_y_continuous(limits = c(-100,-20), breaks=c(-20,-40,-60,-80,-100)) +
  facet_wrap(~substock) +
  theme_bw()
ggsave(file="CV_reduction.pdf", width=7, height=5)



MyColour <- c("#ffc431", "#a9fee0", "#ffc7ae", "#288baf", "#f2dc3b", "#dd8071", "#cadeae", "#d2adea")
ggplot(ddcv3, aes(x = bin, y = CV_reduction, color = substock, group=substock)) +
  geom_point(size = 2) +
  geom_line() +
  xlab("Sample size bin") +
  ylab("Percent reduction in Coefficient of Variation") +
  scale_color_manual(values = MyColour) +
  facet_wrap(~year) +
  theme_bw()





# run-timing effect ####
require(tidyverse)

d <- read.csv(file.choose()) # input_2012-2016_fullsampsizetet.csv
x <- as.factor(c(15:350))

df <- d %>%
  filter(prob >= 0.51) %>%
  filter(year >= 2014) %>%
  filter(julian_date %in% x) %>%
  filter(julian_date != "unknwn") %>%
  mutate(substock = case_when(region == "30" ~ "Upper Lakes and Mainstem",
                              region == "31" ~ "Teslin River",
                              region == "32" ~ "Carmacks",
                              region == "33" ~ "Middle Mainstem",
                              region == "34" ~ "Pelly",
                              region == "35" ~ "Stewart",
                              region == "36" ~ "Lower Mainstem",
                              region == "38" ~ "White-Donjek"))

# muck with julian dates to create early and late dataframes
jul <- levels(df$julian_date)

jul <- as.vector(levels(df$julian_date))
jul2 <- jul[-56]

jul.e <- jul2[1:27]
jul.l <- jul2[28:55]

# get values for early julian date substock proportions
df.e <- df %>% filter(julian_date %in% jul.e)
df.e %>% group_by(year) %>% summarise(count = n())
df.e %>% group_by(year, substock) %>% summarise(count = n())

# create early data frame
df.e <- df %>%
  filter(julian_date %in% jul.e) %>%
  group_by(year, substock) %>%
  mutate(count = n()) %>%
  mutate(bin = case_when(year == "2014" ~ 393,
                         year == "2015" ~ 441,
                         year == "2016" ~ 313)) %>%
  mutate(freq = count/bin) %>%
  distinct(year, substock, count, bin, freq) %>%
  mutate(run.timing = "early")


# get values for late julian date substock proportions
df.l <- df %>% filter(julian_date %in% jul.l)
df.l %>% group_by(year) %>% summarise(count = n())
df.l %>% group_by(year, substock) %>% summarise(count = n())

# create early data frame
df.l <- df %>%
  filter(julian_date %in% jul.l) %>%
  group_by(year, substock) %>%
  mutate(count = n()) %>%
  mutate(bin = case_when(year == "2014" ~ 97,
                         year == "2015" ~ 256,
                         year == "2016" ~ 145)) %>%
  mutate(freq = count/bin) %>%
  distinct(year, substock, count, bin, freq) %>%
  mutate(run.timing = "late")

# join early and late dataframes
dff <- full_join(df.e, df.l)


# get values for full run duration julian date substock proportions
df.f <- df %>% filter(julian_date %in% jul2)
df.f %>% group_by(year) %>% summarise(count = n())
df.f %>% group_by(year, substock) %>% summarise(count = n())

# create full run dataframe
df.f <- df %>%
  filter(julian_date %in% jul2) %>%
  group_by(year, substock) %>%
  mutate(count = n()) %>%
  mutate(bin = case_when(year == "2014" ~ 490,
                         year == "2015" ~ 697,
                         year == "2016" ~ 458)) %>%
  mutate(freq = count/bin) %>%
  distinct(year, substock, count, bin, freq) %>%
  mutate(run.timing = "full")

# join dataframes together
df2 <- full_join(dff, df.f)

# re-order run timing factors
df2$run.timing <- factor(df2$run.timing, levels=c("full", "early", "late"))

ddf2 <- df2 %>%
  select(year, substock, freq, run.timing) %>%
  spread(key = run.timing, value = freq) %>%
  mutate(early_diff = (early - full)) %>%
  mutate(late_diff = (late - full)) %>%
  mutate(early_bias = (early_diff/full)*100) %>%
  mutate(late_bias = (late_diff/full)*100) %>%
  gather(key=Run.timing, value = value, early_bias:late_bias)

# FIGURE, run timing percent bias, sampline errors memo
name <- factor(ddf2$substock, c("Lower Mainstem", "White-Donjek", "Pelly", "Stewart",
                                "Carmacks", "Teslin River", "Middle Mainstem", "Upper Lakes and Mainstem"))

ggplot(ddf2, aes(x=name, y=value/100, fill=Run.timing)) +
  geom_boxplot(position = "dodge") +
  xlab("Substock") +
  ylab("Percent bias") +
  scale_y_continuous(labels= percent, breaks = seq(-1,1, by=0.25)) +
  scale_fill_grey(start=0.8, end=0.5) +
  theme_bw() +
  theme(axis.title = element_text(size=12)) +
  theme(axis.text = element_text(size=11)) +
  theme(axis.text.x = element_text(angle=45, hjust=1))
ggsave(file="3_FIG_3.run.time.bias.pdf", width = 7, height = 5)


# FIGURE; combined years substock run timing
dd <- read.csv(file.choose()) # input_GSI_master.csv
ddd <- dd %>%
  filter(prob > 0.5) %>%
  filter(gear == "Fish Wheel") %>% 
  filter(julian != "unknwn") %>%
  mutate(substock = case_when(region == "30" ~ "Upper Lakes and Mainstem",
                              region == "31" ~ "Teslin River",
                              region == "32" ~ "Carmacks",
                              region == "33" ~ "Middle Mainstem",
                              region == "34" ~ "Pelly",
                              region == "35" ~ "Stewart",
                              region == "36" ~ "Lower Mainstem",
                              region == "38" ~ "White-Donjek"))

ggplot(ddd, aes(x=fct_reorder(substock, julian, fun = median, .desc = FALSE), y=julian)) +
  geom_boxplot(size=1.05) +
  xlab("Substock") +
  ylab("Julian Day") +
  theme_bw() +
  theme(axis.title = element_text(size=12)) +
  theme(axis.text = element_text(size=11)) +
  theme(axis.text.x = element_text(angle=45, hjust=1))
ggsave(file="FIG_2_run.timing.all.pdf", width = 7, height = 5)




# gear effect ####
d <- read.csv(file.choose()) # input_master_GSI.csv

master2 <- d %>%
  filter(region != "NA") %>%
  filter(prob > 0.5) %>%
  filter(year == 2008 | year == 2010 | year == 2011 | year == 2012) %>%
  mutate(substock = case_when(region == "30" ~ "Upper Lakes and Mainstem",
                              region == "31" ~ "Teslin River",
                              region == "32" ~ "Carmacks",
                              region == "33" ~ "Middle Mainstem",
                              region == "34" ~ "Pelly",
                              region == "35" ~ "Stewart",
                              region == "36" ~ "Lower Mainstem",
                              region == "38" ~ "White-Donjek"))



# get julian dates for FW samples

jul_dates <- master2 %>%
  filter(gear == "Fish Wheel") %>%
  group_by(year) %>%
  select(julian) %>%
  distinct()

julian_2008 <- jul_dates %>%
  filter(year == 2008) %>%
  pull(julian)

julian_2010 <- jul_dates %>%
  filter(year == 2010) %>%
  pull(julian)

julian_2011 <- jul_dates %>%
  filter(year == 2011) %>%
  pull(julian)

julian_2012 <- jul_dates %>%
  filter(year == 2012) %>%
  pull(julian)


# Test Fishery Full sample substock proportions
TF_full_prop <- master2 %>%
  group_by(year, gear) %>%
  mutate(count = n()) %>%
  group_by(year, gear, substock) %>%
  mutate(count_ind = n()) %>%
  mutate(tf_substock_run_prop = count_ind/count) %>%
  group_by(year, gear) %>% 
  distinct(substock, tf_substock_run_prop) %>%
  #replace(., is.na(.), 0) %>% # change NA to zero
  filter(gear == "Test Fishery") %>%
  arrange(year, substock) %>%
  ungroup() %>%
  select(-gear)



# Fish Wheel Full sample substock proportions
FW_full_prop <- master2 %>%
  group_by(year, gear) %>%
  mutate(count = n()) %>%
  group_by(year, gear, substock) %>%
  mutate(count_ind = n()) %>%
  mutate(fw_substock_run_prop = count_ind/count) %>%
  group_by(year, gear) %>% 
  distinct(substock, fw_substock_run_prop) %>%
  #replace(is.na(.), 0) %>% # change NA to zero
  filter(gear == "Fish Wheel")  %>%
  arrange(year, substock) %>%
  ungroup() %>%
  select(-gear)

# Test Fishery reduced julian day substock proportions TF_red_prop
TF_red_2008 <- master2 %>%
  filter(year == 2008) %>%
  filter(julian %in% julian_2008)
  
TF_red_2010 <- master2 %>%
  filter(year == 2010) %>%
  filter(julian %in% julian_2010)

TF_red_2011 <- master2 %>%
  filter(year == 2011) %>%
  filter(julian %in% julian_2011)

TF_red_2012 <- master2 %>%
  filter(year == 2012) %>%
  filter(julian %in% julian_2012)
  
TF_red_prop1 <- rbind(TF_red_2008, TF_red_2010)
TF_red_prop2 <- rbind(TF_red_prop1, TF_red_2011)
TF_red_prop3 <- rbind(TF_red_prop2, TF_red_2012)


TF_red_prop <- TF_red_prop3 %>%
  group_by(year, gear) %>%
  mutate(count = n()) %>%
  group_by(year, gear, substock) %>%
  mutate(count_ind = n()) %>%
  mutate(tf_red_substock_run_prop = count_ind/count) %>%
  group_by(year, gear) %>% 
  distinct(substock, tf_red_substock_run_prop) %>%
  #replace(., is.na(.), 0) %>% # change NA to zero
  filter(gear == "Test Fishery") %>%
  arrange(year, substock) %>%
  ungroup() %>%
  select(-gear)


compare_df <- full_join(TF_red_prop, FW_full_prop) %>%
  replace(is.na(.), 0) # change NA to zero

# FIGURE: fish wheel and test fishery matching julian dates substock prop
head(compare_df)

# list of colors that matches map colors
# "#ffc431", orange, lower mainstem
# "#ffc7ae", pink, pelly
# "#a9fee0", aquamarine, stewart
# "#dd8071", red, carmacks
# "#f2dc3b", yellow, middle mainstem
# "#cadeae", moss green, upper lakes
# "#288baf", dark blue, white-donjek
# "#d2adea", purple, teslin


# create color scheme
MyColour <- c("#ffc431", "#a9fee0", "#ffc7ae", "#288baf", "#f2dc3b", "#dd8071", "#cadeae", "#d2adea")

ggplot(compare_df, aes(x=tf_red_substock_run_prop, y=fw_substock_run_prop, colour=substock)) +
  geom_point(size=4) +
  geom_abline(slope = 1, intercept = 0, lty = "dotted") +
  scale_color_manual(values = MyColour) +
  xlab("Sub-stock proportion: Test Fishery") +
  ylab("Sub-stock proportion: Fish Wheel") +
  facet_wrap(~year, nrow=2) +
  theme_bw() +
  theme(axis.text = element_text(size=10)) +
  theme(axis.title = element_text(size=12))
ggsave(file="compare_TFred_FW_1_1.pdf", width=7, height=5)



# create bias data frame
bias_df <- full_join(TF_full_prop, FW_full_prop) %>%
  replace(is.na(.), 0) # change NA to zero


bias_df2 <- full_join(bias_df, TF_red_prop) %>%
  replace(is.na(.), 0)  %>% # change NA to zero
  mutate(FW_diff = tf_substock_run_prop - fw_substock_run_prop) %>%
  mutate(TF_red_diff = tf_substock_run_prop - tf_red_substock_run_prop) %>%
  ungroup() %>%
  group_by(year) %>%
  mutate(FW_diff_sum = abs(sum(FW_diff, na.rm = TRUE))) %>%
  mutate(FW_mean_bias = FW_diff_sum/8) %>%
  mutate(TF_red_diff_sum = abs(sum(TF_red_diff, na.rm = TRUE))) %>%
  mutate(TF_red_mean_bias = TF_red_diff_sum/8) %>%
  select(year, substock, FW_diff, TF_red_diff) %>%
  gather(key = label, value = value, FW_diff:TF_red_diff)


# Sole effect of Fish wheel
bias_df3 <- full_join(bias_df, TF_red_prop) %>%
  replace(is.na(.), 0) %>%
  mutate(FW_diff = tf_substock_run_prop - fw_substock_run_prop,
         TF_red_diff = tf_substock_run_prop - tf_red_substock_run_prop,
         FW_TFred_diff = tf_red_substock_run_prop - fw_substock_run_prop,
         FW_bias = (FW_diff/tf_substock_run_prop)*100,
         FW_TFred_bias = (FW_TFred_diff/tf_red_substock_run_prop)*100)


# summary of % bias
bias_summary <- bias_df3 %>%
  group_by(substock) %>%
  mutate(median = median(FW_TFred_diff*100)) %>%
  mutate(min = min(FW_TFred_diff*100)) %>%
  mutate(max = max(FW_TFred_diff*100)) %>%
  mutate(mean_error = mean(abs(FW_TFred_diff*100))) %>%
  select(substock, mean_error, median, min, max) %>%
  distinct(substock, mean_error, median, min, max)



# FIGURE, bias across years by substock sole effect of fish wheel, matching julian dates with FW
ggplot(bias_df3, aes(x=substock, y=FW_TFred_bias/100)) +
  geom_boxplot(size=1.05) +
  xlab("Sub-stock") +
  ylab("Fish Wheel % bias") +
  scale_y_continuous(labels = percent) +
  geom_hline(yintercept = 0, lty= "dotted") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1, size=11)) +
  theme(axis.text.y = element_text(size = 11)) +
  theme(axis.title = element_text(size = 12))
ggsave(file= "FIG_1_samp.memo.FWbias.pdf", width=7, height = 5)



# plot difference (from Test Fishery full sample collection) in substock proportions 
# for Fish Wheel and Test Fishery reduced sample collections
ggplot(bias_df2, aes(x=substock, y=value*100, fill=label)) +
  geom_bar(stat = "identity", position = "dodge") +
  ylab("Diff. in run prop. from full test fishery estimate") +
  xlab("Sub-stock") +
  facet_wrap(~year, nrow=2, scales = "free_y") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust=1))
ggsave(file="FW.julianmatching.TF.sub.diff.pdf", width=7, height=5)


# calculate bias
bias_df2 %>%
  mutate(abs_value = abs(value)) %>%
  group_by(year, label) %>%
  summarise(absolute_error = (sum(abs_value)*100)) %>%
  mutate(mean = absolute_error/8) %>%
  mutate(min = min(absolute_error)) %>%
  mutate(max = max(absolute_error)) %>%
  ungroup() %>%
  select(year, label, min, max) %>%
  distinct(year, label, min, max)


# calculate how well run has been sampled early/late bins ####

# read in border counts input file
border_counts <- read.delim(file="input_bordercounts2.txt",header=TRUE)
head(border_counts)

# get total number of border passage counts for each year, and early and late border passage counts
bord_cnt2 <- border_counts %>%
  replace(is.na(.), 0) %>%
  filter(count > 0) %>%
  group_by(year, count_type) %>%
  mutate(count_tot = sum(count)) %>%
  mutate(count_half = round(count_tot/2, digits = 0)) %>%
  ungroup() %>%
  group_by(year, count_type) %>%
  mutate(cuml_early = cumsum(count)) %>%
  mutate(ce_diff = abs(count_half - cuml_early)) %>%
  mutate(min_diff = min(ce_diff)) %>%
  filter(ce_diff == min_diff) %>%
  mutate(cuml_late = count_tot - cuml_early) %>%
  mutate(early_prop = (cuml_early/count_tot)*100) %>%
  mutate(late_prop = abs(100 - early_prop)) %>%
  select(count_type, year, julian, count_tot, 
         cuml_early, cuml_late, early_prop, late_prop)



# plot border passage counts (as % of total run) per run timing bin
ggplot(bord_cnt2 %>% gather(key = run.timing, value = run.prop, early_prop:late_prop), 
       aes(x=run.timing, y=run.prop, fill=run.timing)) +
  geom_bar(stat="identity", position= "dodge") +
  xlab("Run Timing Bin") +
  ylab("Border Passage Count (%)") +
  facet_wrap(~year) +
  theme_bw()



# read in GSI data input file
data <- read.csv(file.choose()) # input_master_GSI.csv

# filter so only individuals with stock id > 0.5 are included
# additional columns include the count of samples that correspond to
# ‘early’ and ‘late’ run timing bins, the Julian date range for the annual
# ‘early’ and ‘late’ bins were created by taking the julian date from the 
# border_counts dataframe where ~50% of the fish have been counted, 



head(data)

data2 <- data %>%
  filter(prob > 0.5) %>%
  mutate(substock = case_when(region == 30 ~ "Upper Lakes and Mainstem",
                              region == 31 ~ "Teslin River",
                              region == 32 ~ "Carmacks",
                              region == 33 ~ "Middle Mainstem",
                              region == 34 ~ "Pelly",
                              region == 35 ~ "Stewart",
                              region == 36 ~ "Lower Mainstem",
                              region == 38 ~ "White-Donjek")) %>%
  group_by(year, gear) %>%
  mutate(year_samp_count = n()) %>%
  mutate(samp_early_count = case_when(year == 1985 ~ length(julian[julian <= 213]),
                                      year == 1986 ~ length(julian[julian <= 203]), 
                                      year == 1987 ~ length(julian[julian <= 203]), 
                                      year == 1988 ~ length(julian[julian <= 206]), 
                                      year == 1989 ~ length(julian[julian <= 205]), 
                                      year == 1990 ~ length(julian[julian <= 202]), 
                                      year == 1991 ~ length(julian[julian <= 206]), 
                                      year == 1992 ~ length(julian[julian <= 208]), 
                                      year == 1993 ~ length(julian[julian <= 199]), 
                                      year == 1994 ~ length(julian[julian <= 201]), 
                                      year == 1995 ~ length(julian[julian <= 206]), 
                                      year == 1996 ~ length(julian[julian <= 197]), 
                                      year == 1997 ~ length(julian[julian <= 207]), 
                                      year == 1998 ~ length(julian[julian <= 210]), 
                                      year == 1999 ~ length(julian[julian <= 207]),
                                      year == 2000 ~ length(julian[julian <= 206]), 
                                      year == 2001 ~ length(julian[julian <= 207]), 
                                      year == 2002 ~ length(julian[julian <= 208]), 
                                      year == 2003 ~ length(julian[julian <= 196]), 
                                      year == 2004 ~ length(julian[julian <= 200]), 
                                      year == 2005 ~ length(julian[julian <= 203]), 
                                      year == 2006 ~ length(julian[julian <= 207]), 
                                      year == 2007 ~ length(julian[julian <= 203]), 
                                      year == 2008 ~ length(julian[julian <= 209]), 
                                      year == 2009 ~ length(julian[julian <= 202]), 
                                      year == 2010 ~ length(julian[julian <= 208]), 
                                      year == 2011 ~ length(julian[julian <= 205]), 
                                      year == 2012 ~ length(julian[julian <= 210]), 
                                      year == 2013 ~ length(julian[julian <= 205]), 
                                      year == 2014 ~ length(julian[julian <= 199]), 
                                      year == 2015 ~ length(julian[julian <= 203]), 
                                      year == 2016 ~ length(julian[julian <= 199]))) %>%
  mutate(samp_late_count = case_when(year == 1985 ~ length(julian[julian > 213]),
                                     year == 1986 ~ length(julian[julian > 203]), 
                                     year == 1987 ~ length(julian[julian > 203]), 
                                     year == 1988 ~ length(julian[julian > 206]), 
                                     year == 1989 ~ length(julian[julian > 205]), 
                                     year == 1990 ~ length(julian[julian > 202]), 
                                     year == 1991 ~ length(julian[julian > 206]), 
                                     year == 1992 ~ length(julian[julian > 208]), 
                                     year == 1993 ~ length(julian[julian > 199]), 
                                     year == 1994 ~ length(julian[julian > 201]), 
                                     year == 1995 ~ length(julian[julian > 206]), 
                                     year == 1996 ~ length(julian[julian > 197]), 
                                     year == 1997 ~ length(julian[julian > 207]), 
                                     year == 1998 ~ length(julian[julian > 210]), 
                                     year == 1999 ~ length(julian[julian > 207]),
                                     year == 2000 ~ length(julian[julian > 206]), 
                                     year == 2001 ~ length(julian[julian > 207]), 
                                     year == 2002 ~ length(julian[julian > 208]), 
                                     year == 2003 ~ length(julian[julian > 196]), 
                                     year == 2004 ~ length(julian[julian > 200]), 
                                     year == 2005 ~ length(julian[julian > 203]), 
                                     year == 2006 ~ length(julian[julian > 207]), 
                                     year == 2007 ~ length(julian[julian > 203]), 
                                     year == 2008 ~ length(julian[julian > 209]), 
                                     year == 2009 ~ length(julian[julian > 202]), 
                                     year == 2010 ~ length(julian[julian > 208]), 
                                     year == 2011 ~ length(julian[julian > 205]), 
                                     year == 2012 ~ length(julian[julian > 210]), 
                                     year == 2013 ~ length(julian[julian > 205]), 
                                     year == 2014 ~ length(julian[julian > 199]), 
                                     year == 2015 ~ length(julian[julian > 203]), 
                                     year == 2016 ~ length(julian[julian > 199]))) %>%
  mutate(samp_early_prop = (samp_early_count/year_samp_count)*100) %>%
  mutate(samp_late_prop = 100-samp_early_prop) %>%
  mutate(median_julian = case_when(year == 1985 ~  213,
                                   year == 1986 ~  203, 
                                   year == 1987 ~  203, 
                                   year == 1988 ~  206, 
                                   year == 1989 ~  205, 
                                   year == 1990 ~  202, 
                                   year == 1991 ~  206, 
                                   year == 1992 ~  208, 
                                   year == 1993 ~  199, 
                                   year == 1994 ~  201, 
                                   year == 1995 ~  206, 
                                   year == 1996 ~  197, 
                                   year == 1997 ~  207, 
                                   year == 1998 ~  210, 
                                   year == 1999 ~  207,
                                   year == 2000 ~  206, 
                                   year == 2001 ~  207, 
                                   year == 2002 ~  208, 
                                   year == 2003 ~  196, 
                                   year == 2004 ~  200, 
                                   year == 2005 ~  203, 
                                   year == 2006 ~  207, 
                                   year == 2007 ~  203, 
                                   year == 2008 ~  209, 
                                   year == 2009 ~  202, 
                                   year == 2010 ~  208, 
                                   year == 2011 ~  205, 
                                   year == 2012 ~  210, 
                                   year == 2013 ~  205, 
                                   year == 2014 ~  199, 
                                   year == 2015 ~  203, 
                                   year == 2016 ~  199)) %>%
  mutate(deviation_julian = julian - median_julian) %>%
  select(data_label, year, gear, julian, median_julian, deviation_julian, samp_num, 
         region, prob, samp_early_count, samp_late_count, samp_early_prop, samp_late_prop) %>%
  as.data.frame()

head(data2)

# create output file: which is the input for the bias multinomial
write.table(data2, file = "input_bias_multinomial.csv", sep=" ")


# plot sample counts (as % of total run) per run timing bin
data3 <- data2 %>%
  filter(year <= 2005 | year > 2005 & gear == "Test Fishery" ) %>%
  gather(key = run.timing, value = run.prop, samp_early_prop:samp_late_prop)
  
ggplot(data3, aes(x=run.timing, y=run.prop, fill=run.timing)) +
  geom_bar(stat="identity", position= "dodge") +
  xlab("Run Timing Bin") +
  ylab("Sample Passage Count (%)") +
  facet_wrap(~year) +
  theme_bw() +
  theme(axis.text.x = element_text(angle=90, hjust=1))












