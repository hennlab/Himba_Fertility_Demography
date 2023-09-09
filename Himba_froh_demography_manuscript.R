# Code used in Swinford et al: 
# Increased homozygosity due to Endogamy results in fitness consequences in a human population
# NAS

library(dplyr)
library(reshape2)
library(ggplot2)
library(gridExtra)
library(ggplotify)
library(Metrics)
library(tidyverse)
library(brms)
library(tidybayes)
library(modelr)
library(patchwork)

### FROH DISTRIBUTIONS ###
# Read in ROH files
H3Africa_auto_kb500 <- read.table("~/Desktop/H3Africa_n681_NoChr23_thinned_snp50_missing2_het1_kb500.hom.indiv", header=TRUE)
H3Africa_auto_kb1500 <- read.table("~/Desktop/H3Africa_n681_NoChr23_thinned_snp50_missing2_het1_kb1500.hom.indiv", header=TRUE)
H3Africa_auto_kb5000 <- read.table("~/Desktop/H3Africa_n681_NoChr23_thinned_snp50_missing2_het1_kb5000.hom.indiv", header=TRUE)
Mega_auto_kb500 <- read.table("~/Desktop/MegaEx_NoChr23_snp50_missing2_het1_kb500.hom.indiv", header=TRUE)
Mega_auto_kb1500 <- read.table("~/Desktop/MegaEx_NoChr23_snp50_missing2_het1_kb1500.hom.indiv", header=TRUE)
Mega_auto_kb5000 <- read.table("~/Desktop/MegaEx_NoChr23_snp50_missing2_het1_kb5000.hom.indiv", header=TRUE)

# Combine ROH columns
H3Africa_auto_kb500 <- H3Africa_auto_kb500 %>%
  select(IID, NSEG, KB, KBAVG) %>%
  rename(NSEG_500 = NSEG) %>%
  rename(KB_500 = KB) %>%
  rename(KBAVG_500 = KBAVG)
H3Africa_auto_kb1500 <- H3Africa_auto_kb1500 %>%
  select(IID, NSEG, KB, KBAVG) %>%
  rename(IID1 = IID) %>%
  rename(NSEG_1500 = NSEG) %>%
  rename(KB_1500 = KB) %>%
  rename(KBAVG_1500 = KBAVG)
H3Africa_auto_kb5000 <- H3Africa_auto_kb5000 %>%
  select(IID, NSEG, KB, KBAVG) %>%
  rename(IID2 = IID) %>%
  rename(NSEG_5000 = NSEG) %>%
  rename(KB_5000 = KB) %>%
  rename(KBAVG_5000 = KBAVG)
H3Africa_auto_roh <- cbind(H3Africa_auto_kb500, H3Africa_auto_kb1500, H3Africa_auto_kb5000)
head(H3Africa_auto_roh)
H3Africa_auto_roh <- H3Africa_auto_roh %>%
  select(IID, NSEG_500, KB_500, KBAVG_500, NSEG_1500, KB_1500, KBAVG_1500, NSEG_5000, KB_5000, KBAVG_5000)
head(H3Africa_auto_roh)

Mega_auto_kb500 <- Mega_auto_kb500 %>%
  select(IID, NSEG, KB, KBAVG) %>%
  rename(NSEG_500 = NSEG) %>%
  rename(KB_500 = KB) %>%
  rename(KBAVG_500 = KBAVG)
Mega_auto_kb1500 <- Mega_auto_kb1500 %>%
  select(IID, NSEG, KB, KBAVG) %>%
  rename(IID1 = IID) %>%
  rename(NSEG_1500 = NSEG) %>%
  rename(KB_1500 = KB) %>%
  rename(KBAVG_1500 = KBAVG)
Mega_auto_kb5000 <- Mega_auto_kb5000 %>%
  select(IID, NSEG, KB, KBAVG) %>%
  rename(IID2 = IID) %>%
  rename(NSEG_5000 = NSEG) %>%
  rename(KB_5000 = KB) %>%
  rename(KBAVG_5000 = KBAVG)
Mega_auto_roh <- cbind(Mega_auto_kb500, Mega_auto_kb1500, Mega_auto_kb5000)
head(Mega_auto_roh)
Mega_auto_roh <- Mega_auto_roh %>%
  select(IID, NSEG_500, KB_500, KBAVG_500, NSEG_1500, KB_1500, KBAVG_1500, NSEG_5000, KB_5000, KBAVG_5000)
head(Mega_auto_roh)

# Calculate and Add Froh Columns (bim genome length calculated with calcGenomeLength.R)
H3Africa_auto_Froh <- H3Africa_auto_roh %>%
  mutate(Froh_500 = (KB_500*1000)/2791433159) %>%
  mutate(Froh_1500 = (KB_1500*1000)/2791433159) %>%
  mutate(Froh_5000 = (KB_5000*1000)/2791433159)
Mega_auto_Froh <- Mega_auto_roh %>%
  mutate(Froh_500 = (KB_500*1000)/2794139053) %>%
  mutate(Froh_1500 = (KB_1500*1000)/2794139053) %>%
  mutate(Froh_5000 = (KB_5000*1000)/2794139053)

# Combine datasets
Himba_auto_Froh <- rbind(H3Africa_auto_Froh, Mega_auto_Froh)

### PLOT ###
temp <- Himba_auto_Froh %>%
  select(IID, Froh_500, Froh_1500, Froh_5000)
froh.long <- melt(temp, id="IID", variable.name="ROH", value.name="FROH")
levels(froh.long$ROH) <- c("500+ KB", "1500+ KB", "5000+ KB")

cbPalette <- c("#009E73", "#0072B2", "#D55E00") # green, blue, orange-red
# HISTOGRAMS
froh_histograms <- ggplot(data=froh.long, aes(x=FROH, fill=ROH)) + 
  geom_histogram(color="black", alpha=0.5, show.legend=FALSE) + facet_grid(ROH ~ .) + 
  scale_fill_manual(values=cbPalette) 
# SCATTERPLOTS
scatter1 <- Himba_auto_Froh %>% select(IID, NSEG_500, KB_500) %>% 
  mutate(Threshold = "500+ KB") %>% rename(NSEG = NSEG_500, ROH = KB_500)
scatter2 <- Himba_auto_Froh %>% select(IID, NSEG_1500, KB_1500) %>% 
  mutate(Threshold = "1500+ KB") %>% rename(NSEG = NSEG_1500, ROH = KB_1500)
scatter3 <- Himba_auto_Froh %>% select(IID, NSEG_5000, KB_5000) %>% 
  mutate(Threshold = "5000+ KB") %>% rename(NSEG = NSEG_5000, ROH = KB_5000)
nseg_roh_long <- rbind(scatter1, scatter2, scatter3)
nseg_roh_long$Threshold = factor(nseg_roh_long$Threshold, levels=c("500+ KB", "1500+ KB", "5000+ KB"))
nseg_roh_scatterplots <- ggplot(data=nseg_roh_long, aes(x=NSEG, y=ROH, color=Threshold)) + geom_point(alpha=0.5, show.legend=FALSE) + xlab("Number of Segments") +
  facet_grid(Threshold ~ .) + scale_color_manual(values=cbPalette)


#########################################################################################
##### Expectations for offspring of *Himba* First Cousins

Himba_auto_Froh500_belowFirstCuz <- Himba_auto_Froh %>%
  filter(Froh_500 < 0.0625)
Background500 <- mean(Himba_auto_Froh500_belowFirstCuz$Froh_500) 
FirstCuzExp_500 <- Background500 + 0.0625  
Himba_auto_Froh %>%
  filter(Froh_500 >= FirstCuzExp_500) # Only 1 individual meets the expectation (with Froh=0.1082)


#########################################################################################
##### THIN TEST & RMSE #####

setwd("~/Desktop/UC Davis/HENN LAB/Himba Project/F_ROH/DifFroh_Analysis")

thinned_1 <- read.table("1_ThinTest_snp50_missing2_het1_kb1500.hom.indiv", header=TRUE)
thinned_2 <- read.table("2_ThinTest_snp50_missing2_het1_kb1500.hom.indiv", header=TRUE)
thinned_3 <- read.table("3_ThinTest_snp50_missing2_het1_kb1500.hom.indiv", header=TRUE)
thinned_4 <- read.table("4_ThinTest_snp50_missing2_het1_kb1500.hom.indiv", header=TRUE)
thinned_5 <- read.table("5_ThinTest_snp50_missing2_het1_kb1500.hom.indiv", header=TRUE)
thinned_6 <- read.table("6_ThinTest_snp50_missing2_het1_kb1500.hom.indiv", header=TRUE)
thinned_7 <- read.table("7_ThinTest_snp50_missing2_het1_kb1500.hom.indiv", header=TRUE)
thinned_8 <- read.table("8_ThinTest_snp50_missing2_het1_kb1500.hom.indiv", header=TRUE)
thinned_9 <- read.table("9_ThinTest_snp50_missing2_het1_kb1500.hom.indiv", header=TRUE)
thinned_10 <- read.table("10_ThinTest_snp50_missing2_het1_kb1500.hom.indiv", header=TRUE)

# denominators calculated with calcGenomeLength.R
thin1_froh <- thinned_1 %>%
  mutate(Froh1_1500 = (KB*1000)/2787160584)
thin2_froh <- thinned_2 %>%
  mutate(Froh2_1500 = (KB*1000)/2789881009)
thin3_froh <- thinned_3 %>%
  mutate(Froh3_1500 = (KB*1000)/2790651978)
thin4_froh <- thinned_4 %>%
  mutate(Froh4_1500 = (KB*1000)/2790821199)
thin5_froh <- thinned_5 %>%
  mutate(Froh5_1500 = (KB*1000)/2790972299)
thin6_froh <- thinned_6 %>%
  mutate(Froh6_1500 = (KB*1000)/2790291177)
thin7_froh <- thinned_7 %>%
  mutate(Froh7_1500 = (KB*1000)/2787645606)
thin8_froh <- thinned_8 %>%
  mutate(Froh8_1500 = (KB*1000)/2791488892)
thin9_froh <- thinned_9 %>%
  mutate(Froh9_1500 = (KB*1000)/2787721777)
thin10_froh <- thinned_10 %>%
  mutate(Froh10_1500 = (KB*1000)/2790940619)

ThinTest_Froh <- H3Africa_auto_Froh %>%
  select(IID, Froh_1500)
ThinTest_Froh$Froh1_1500 = thin1_froh$Froh1_1500[match(ThinTest_Froh$IID, thin1_froh$IID)]
ThinTest_Froh$Froh2_1500 = thin2_froh$Froh2_1500[match(ThinTest_Froh$IID, thin2_froh$IID)]
ThinTest_Froh$Froh3_1500 = thin3_froh$Froh3_1500[match(ThinTest_Froh$IID, thin3_froh$IID)]
ThinTest_Froh$Froh4_1500 = thin4_froh$Froh4_1500[match(ThinTest_Froh$IID, thin4_froh$IID)]
ThinTest_Froh$Froh5_1500 = thin5_froh$Froh5_1500[match(ThinTest_Froh$IID, thin5_froh$IID)]
ThinTest_Froh$Froh6_1500 = thin6_froh$Froh6_1500[match(ThinTest_Froh$IID, thin6_froh$IID)]
ThinTest_Froh$Froh7_1500 = thin7_froh$Froh7_1500[match(ThinTest_Froh$IID, thin7_froh$IID)]
ThinTest_Froh$Froh8_1500 = thin8_froh$Froh8_1500[match(ThinTest_Froh$IID, thin8_froh$IID)]
ThinTest_Froh$Froh9_1500 = thin9_froh$Froh9_1500[match(ThinTest_Froh$IID, thin9_froh$IID)]
ThinTest_Froh$Froh10_1500 = thin10_froh$Froh10_1500[match(ThinTest_Froh$IID, thin10_froh$IID)]
head(ThinTest_Froh)

# PLOT OF THIN TESTS 
plot(ThinTest_Froh$Froh_1500, ThinTest_Froh$Froh1_1500,
     xlab="FROH (original thinned set)",
     ylab="FROH: Thin Tests (n=10)",
     main="Individual Comparison of Multiple Thin Tests (plink)")
points(ThinTest_Froh$Froh_1500, ThinTest_Froh$Froh2_1500, col="gold")
points(ThinTest_Froh$Froh_1500, ThinTest_Froh$Froh3_1500, col="red")
points(ThinTest_Froh$Froh_1500, ThinTest_Froh$Froh4_1500, col="blue")
points(ThinTest_Froh$Froh_1500, ThinTest_Froh$Froh5_1500, col="purple")
points(ThinTest_Froh$Froh_1500, ThinTest_Froh$Froh6_1500, col="pink")
points(ThinTest_Froh$Froh_1500, ThinTest_Froh$Froh7_1500, col="turquoise")
points(ThinTest_Froh$Froh_1500, ThinTest_Froh$Froh8_1500, col="green")
points(ThinTest_Froh$Froh_1500, ThinTest_Froh$Froh9_1500, col="gray50")
points(ThinTest_Froh$Froh_1500, ThinTest_Froh$Froh10_1500, col="magenta")
legend("topleft", legend=c("test 1", "test 2", "test 3", "test 4", "test 5", "test 6", "test 7", "test 8", "test 9", "test 10"), col=c("black", "gold", "red", "blue", "purple", "pink", "turquoise", "green", "gray50", "magenta"), pch=1, cex=0.8)

# Set up for RMSE Calculation
ThinTest_FrohDifs <- ThinTest_Froh %>%
  mutate(Dif1 = abs(Froh1_1500-Froh_1500)) %>%
  mutate(Dif2 = abs(Froh2_1500-Froh_1500)) %>%
  mutate(Dif3 = abs(Froh3_1500-Froh_1500)) %>%
  mutate(Dif4 = abs(Froh4_1500-Froh_1500)) %>%
  mutate(Dif5 = abs(Froh5_1500-Froh_1500)) %>%
  mutate(Dif6 = abs(Froh6_1500-Froh_1500)) %>%
  mutate(Dif7 = abs(Froh7_1500-Froh_1500)) %>%
  mutate(Dif8 = abs(Froh8_1500-Froh_1500)) %>%
  mutate(Dif9 = abs(Froh9_1500-Froh_1500)) %>%
  mutate(Dif10 = abs(Froh10_1500-Froh_1500)) %>%
  mutate(AvgDif = rowMeans(.[, 13:22]))
ExpectedDif <- rep(0, 504)
ThinTest_FrohDifs <- cbind(ThinTest_FrohDifs, ExpectedDif)
## RMSE Calculation 
rmse(ThinTest_FrohDifs$AvgDif, ThinTest_FrohDifs$ExpectedDif)


##########################################################################################
##### POST-REPRODUCTIVE WOMEN 

PRwomen_autoFroh_allFertility <- read.table("~/Desktop/PRwomen_autoFroh_allFertility.txt", sep=" ", header=TRUE)
# Filter out NA's for model
PRwomen_autoFroh_allFertility_NoNA <- PRwomen_autoFroh_allFertility %>%
  filter(TotalSurvive5_withBirths0 != "NA") %>%
  filter(YOB != "NA") %>%
  filter(NumMarriages != "NA")

stand <- function(x) {(x - mean(x,na.rm=TRUE))/sd(x,na.rm=TRUE)} #standardize function
b_greater<-function(x) {
  out<-(sum(x > 0) / length (x))*100
  round(out,1)}
b_lesser<-function(x) {
  out<-(sum(x < 0) / length (x))*100
  round(out,1)}


cols <- c("#7CAE00","#00B4EF","#E7861B")

#all women
d <- PRwomen_autoFroh_allFertility_NoNA
d$yob_s <- stand(d$YOB)
d$Froh_500_s <- stand(d$Froh_500)
d$Froh_1500_s <- stand(d$Froh_1500)
d$Froh_5000_s <- stand(d$Froh_5000)


#only women with any kids
d1 <- PRwomen_autoFroh_FertileWomen
d1$yob_s <- stand(d1$YOB)
d1$Froh_500_s <- stand(d1$Froh_500)
d1$Froh_1500_s <- stand(d1$Froh_1500)
d1$Froh_5000_s <- stand(d1$Froh_5000)

#### recreating models in BRMS 
# ------------------ Model Parameters --------------------#
iter<-4000
warmup<-2000
chains<-3
cores<-chains
# ------------------ ---------------- --------------------#

#set priors
priors <-  c(set_prior("normal(3,3)",class="Intercept"),
             set_prior("normal(0,1)",class="b"))

#plotting
ggplot(d, aes(x = Froh_500, y = TotalSurvive5_withBirths0)) + geom_jitter()


#### All women models ####

hist(d$TotalSurvive5_withBirths0)

all_b <-brm(TotalSurvive5_withBirths0 ~ yob_s + NumMarriages,
            data = d,
            family = poisson(),
            prior = priors,
            chains = chains, cores = cores, iter = iter, warmup = warmup, refresh = 0,
            control=list(adapt_delta = 0.9))



#___looking at 500 
all_500 <-brm(TotalSurvive5_withBirths0 ~ yob_s + NumMarriages + Froh_500_s,
              data = d,
              family = poisson(),
              prior = priors,
              chains = chains, cores = cores, iter = iter, warmup = warmup, refresh = 0,
              control=list(adapt_delta = 0.9))
#post distribution
post<- posterior_samples(all_500)
b_lesser(post$b_Froh_500_s)


#___looking at 1500 
all_1500 <-brm(TotalSurvive5_withBirths0 ~ yob_s + NumMarriages + Froh_1500_s,
               data = d,
               family = poisson(),
               prior = priors,
               chains = chains, cores = cores, iter = iter, warmup = warmup, refresh = 0,
               control=list(adapt_delta = 0.9))
#post distribution
post<- posterior_samples(all_1500)
b_lesser(post$b_Froh_1500_s)

#___looking at 5000 
all_5000 <-brm(TotalSurvive5_withBirths0 ~ yob_s + NumMarriages + Froh_5000_s,
               data = d,
               family = poisson(),
               prior = priors,
               chains = chains, cores = cores, iter = iter, warmup = warmup, refresh = 0,
               control=list(adapt_delta = 0.9))
#post distribution
post<- posterior_samples(all_5000)
b_lesser(post$b_Froh_5000_s)


#### Only Fertile models ####

f_500 <- update(all_500, newdata = d1)
f_1500 <- update(all_1500, newdata = d1)
f_5000 <- update(all_5000, newdata = d1)

#post distribution
post<- posterior_samples(f_500)
b_lesser(post$b_Froh_500_s)

post<- posterior_samples(f_1500)
b_lesser(post$b_Froh_1500_s)

post<- posterior_samples(f_5000)
b_lesser(post$b_Froh_5000_s)

##### Plot Coefficients  ####
#___all women -----
p<-rbind( #get posteriors
  posterior_samples(all_500) %>% select(-lp__) %>% mutate(type = "500") %>% rename(FROH = b_Froh_500_s),
  posterior_samples(all_1500) %>% select(-lp__) %>% mutate(type = "1500") %>% rename(FROH = b_Froh_1500_s),
  posterior_samples(all_5000) %>% select(-lp__) %>% mutate(type = "5000") %>% rename(FROH = b_Froh_5000_s)) %>%
  rename(Intercept = b_Intercept,
         YOB = b_yob_s,
         "Number of marriages" = b_NumMarriages) %>%
  mutate(type = factor(type, levels = c("500","1500","5000")))

p <- p %>% pivot_longer(cols = Intercept:FROH)

coefp1<-
  p %>% 
  ggplot(aes(x = value, y = name, color = type)) + 
  stat_pointinterval(position=position_dodge(0.3)) +
  scale_y_discrete(limits = c("YOB","Number of marriages","FROH"),
                   labels = c("YOB","Number of\nmarriages","FROH")) +
  scale_x_continuous("Coefficient") +
  geom_vline(xintercept = 0) +
  ggtitle("All Women") +
  labs(color = "ROH") +
  theme(axis.title.y = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(values = cols,
                     aesthetics = c("colour", "fill"))


#___only fertile -----
p2<-rbind( #get posteriors
  posterior_samples(f_500) %>% select(-lp__) %>% mutate(type = "500") %>% rename(FROH = b_Froh_500_s),
  posterior_samples(f_1500) %>% select(-lp__) %>% mutate(type = "1500") %>% rename(FROH = b_Froh_1500_s),
  posterior_samples(f_5000) %>% select(-lp__) %>% mutate(type = "5000") %>% rename(FROH = b_Froh_5000_s)) %>%
  rename(Intercept = b_Intercept,
         YOB = b_yob_s,
         "Number of marriages" = b_NumMarriages) %>%
  mutate(type = factor(type, levels = c("500","1500","5000")))

p2 <- p2 %>% pivot_longer(cols = Intercept:FROH)

coefp2<-
  p2 %>% 
  ggplot(aes(x = value, y = name, color = type)) + 
  stat_pointinterval(position=position_dodge(0.3)) +
  scale_y_discrete(limits = c("YOB","Number of marriages","FROH"),
                   labels = c("YOB","Number of\nmarriages","FROH")) +
  scale_x_continuous("Coefficient") +
  geom_vline(xintercept = 0) +
  ggtitle("Fertile Women") +
  theme(axis.title.y = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(values = cols,
                     aesthetics = c("colour", "fill"))




##### Plot outcomes  ####
#___all women -----
#500
d %>%
  data_grid(yob_s = 0,
            NumMarriages = median(d$NumMarriages),
            Froh_500_s = seq_range(Froh_500_s, n = 100)) %>%
  add_epred_draws(all_500) %>%
  ggplot(aes(x = Froh_500_s, y = TotalSurvive5_withBirths0)) +
  stat_lineribbon(aes(y = .epred), alpha =0.8) +
  geom_point(data = d, alpha = 0.5, size = 2) +
  scale_x_continuous("FROH") + 
  scale_fill_brewer(palette = "Greys") 

#1500
d %>%
  data_grid(yob_s = 0,
            NumMarriages = median(d$NumMarriages),
            Froh_1500_s = seq_range(Froh_1500_s, n = 100)) %>%
  add_epred_draws(all_1500) %>%
  ggplot(aes(x = Froh_1500_s, y = TotalSurvive5_withBirths0)) +
  stat_lineribbon(aes(y = .epred), alpha =0.8) +
  geom_point(data = d, alpha = 0.5, size = 2) +
  scale_fill_brewer(palette = "Greys") 

#5000
d %>%
  data_grid(yob_s = 0,
            NumMarriages = median(d$NumMarriages),
            Froh_5000_s = seq_range(Froh_5000_s, n = 100)) %>%
  add_epred_draws(all_5000) %>%
  ggplot(aes(x = Froh_5000_s, y = TotalSurvive5_withBirths0)) +
  stat_lineribbon(aes(y = .epred), alpha =0.8) +
  geom_point(data = d, alpha = 0.5, size = 2) +
  scale_fill_brewer(palette = "Greys") 

#___fertile women ----- match colors with coef plot
#500
li500<- #labels for plot
  c((0-mean(d1$Froh_500))/sd(d1$Froh_500),
    (0.02-mean(d1$Froh_500))/sd(d1$Froh_500),
    (0.04-mean(d1$Froh_500))/sd(d1$Froh_500),
    (0.06-mean(d1$Froh_500))/sd(d1$Froh_500),
    (0.08-mean(d1$Froh_500))/sd(d1$Froh_500),
    (0.1-mean(d1$Froh_500))/sd(d1$Froh_500))

p500<-
  d1 %>%
  data_grid(yob_s = 0,
            NumMarriages = median(d$NumMarriages),
            Froh_500_s = seq_range(Froh_500_s, n = 100)) %>%
  add_epred_draws(f_500) %>%
  ggplot(aes(x = Froh_500_s, y = TotalSurvive5_withBirths0)) +
  stat_lineribbon(aes(y = .epred), alpha =0.5) +
  geom_point(data = d1, alpha = 0.5, size = 2) +
  scale_y_continuous("Completed Fertility") +
  scale_x_continuous("FROH 500",breaks = li500,
                     labels = c(0,0.02,0.04,0.06,0.08,0.1)) +
  scale_color_manual(values = rep(cols[1], 3),
                     aesthetics = c("colour", "fill")) +
  theme(legend.position = "none")

#1500
li1500<- #labels for plot
  c((0-mean(d1$Froh_1500))/sd(d1$Froh_1500),
    (0.02-mean(d1$Froh_1500))/sd(d1$Froh_1500),
    (0.04-mean(d1$Froh_1500))/sd(d1$Froh_1500),
    (0.06-mean(d1$Froh_1500))/sd(d1$Froh_1500),
    (0.08-mean(d1$Froh_1500))/sd(d1$Froh_1500),
    (0.1-mean(d1$Froh_1500))/sd(d1$Froh_1500))

p1500<-
  d1 %>%
  data_grid(yob_s = 0,
            NumMarriages = median(d$NumMarriages),
            Froh_1500_s = seq_range(Froh_1500_s, n = 100)) %>%
  add_epred_draws(f_1500) %>%
  ggplot(aes(x = Froh_1500_s, y = TotalSurvive5_withBirths0)) +
  stat_lineribbon(aes(y = .epred), alpha =0.5) +
  geom_point(data = d1, alpha = 0.5, size = 2) +
  scale_y_continuous("Completed Fertility") +
  scale_x_continuous("FROH 1500",breaks = li1500,
                     labels = c(0,0.02,0.04,0.06,0.08,0.1)) +
  scale_color_manual(values = rep(cols[2], 3),
                     aesthetics = c("colour", "fill")) +
  theme(legend.position = "none")

#5000
li5000<- #labels for plot
  c((0-mean(d1$Froh_5000))/sd(d1$Froh_5000),
    (0.02-mean(d1$Froh_5000))/sd(d1$Froh_5000),
    (0.04-mean(d1$Froh_5000))/sd(d1$Froh_5000),
    (0.06-mean(d1$Froh_5000))/sd(d1$Froh_5000),
    (0.08-mean(d1$Froh_5000))/sd(d1$Froh_5000),
    (0.1-mean(d1$Froh_5000))/sd(d1$Froh_5000))

p5000<-
  d1 %>%
  data_grid(yob_s = 0,
            NumMarriages = median(d$NumMarriages),
            Froh_5000_s = seq_range(Froh_5000_s, n = 100)) %>%
  add_epred_draws(f_5000) %>%
  ggplot(aes(x = Froh_5000_s, y = TotalSurvive5_withBirths0)) +
  stat_lineribbon(aes(y = .epred), alpha =0.5) +
  geom_point(data = d1, alpha = 0.5, size = 2) +
  scale_y_continuous("Completed Fertility") +
  scale_x_continuous("FROH 5000",breaks = li5000,
                     labels = c(0,0.02,0.04,0.06,0.08,0.1)) +
  scale_color_manual(values = rep(cols[3], 3),
                     aesthetics = c("colour", "fill")) +
  theme(legend.position = "none")


#___plot and save -----
#coefficients
coefp2x<- coefp2 + theme(legend.position = "none", axis.text.y=element_blank())
coef <- coefp1 + coefp2x + plot_layout(guides = "collect")
coef

ggsave(filename = "Coefficient_plot.pdf", plot = last_plot(), dpi = 300, width=6, height = 3, units ="in")

#model predictions
preds <-(p500 | p1500 + ylab(NULL) + theme(axis.title.y = element_blank()) | p5000 + theme(axis.title.y = element_blank()))
preds

ggsave(filename = "model_predictions.pdf", plot = last_plot(), dpi = 300, width = 7, height = 4, units  = "in")


#plot them together - doesn't quite look right, save for later
wrap_plots(coef,preds, ncol = 2, nrow = 1)


###model comparison ####
all_b <- add_criterion(all_b, c("loo","waic"))
all_500 <- add_criterion(all_500, c("loo","waic"))
all_1500 <- add_criterion(all_1500, c("loo","waic"))
all_5000 <- add_criterion(all_5000, c("loo","waic"))

loo_compare(all_b,all_500) %>% print(simplify = F) #adding FROH always increases model fit
loo_compare(all_b,all_1500) %>% print(simplify = F)
loo_compare(all_b,all_5000) %>% print(simplify = F)


bayes_R2(all_500) %>% round(digits = 3)

#looking at R2
rbind(bayes_R2(all_b), 
      bayes_R2(all_500), 
      bayes_R2(all_1500), 
      bayes_R2(all_5000)) %>%
  data.frame() %>%
  mutate(model= c("b","500", "1500", "5000"),
         r_square_posterior_mean = round(Estimate, digits = 3)) %>%
  select(model, r_square_posterior_mean)



####outlier analysis ####
#just going to use all women here
library(loo)
#___500 ----

loo(all_500) #one bad value
loo(all_500) %>% pareto_k_ids(threshold=0.5) #which one?
d %>% slice(loo::loo(all_500) %>% pareto_k_ids(threshold=0.5))  #this one
pareto_k_values(loo(all_500))[loo(all_500) %>% pareto_k_ids(threshold=0.5)] #how bad?

#in the data
d %>%
  data_grid(yob_s = 0,
            NumMarriages = median(d$NumMarriages),
            Froh_500_s = seq_range(Froh_500_s, n = 100)) %>%
  add_epred_draws(all_500) %>%
  ggplot(aes(x = Froh_500_s, y = TotalSurvive5_withBirths0)) +
  stat_lineribbon(aes(y = .epred), alpha =0.8) +
  geom_point(data = d, alpha = 0.5, size = 2) +
  geom_point(data =  d %>% filter(IID == "HMB032"), aes(x= Froh_500_s, y = TotalSurvive5_withBirths0), color = "red") +
  scale_x_continuous("FROH") + 
  scale_fill_brewer(palette = "Greys") 

#plotting pareto_k
tibble(pareto_k = all_500$criteria$loo$diagnostics$pareto_k,
       p_waic   = all_500$criteria$waic$pointwise[, "p_waic"],
       Loc      = pull(d, IID)) %>% 
  ggplot(aes(x = pareto_k, y = p_waic, color = Loc == "IID")) +
  geom_vline(xintercept = .5, linetype = 2, color = "black", alpha = 1/2) +
  geom_point(aes(shape = Loc == "IID")) +
  geom_text(data = . %>% filter(p_waic > 0.5),
            aes(x = pareto_k - 0.03, label = Loc),
            hjust = 1) +
  scale_color_manual(values = cols) +
  labs(subtitle = "FROH 500 Model Diagnostic") +
  theme(legend.position = "none")

#rerun without this point and check if predictor is still meaningful
all_500_2 <- update(all_500, newdata = d %>% filter(IID != "HMB032"))

b_lesser(posterior_samples(all_500)$b_Froh_500_s) #old pr(b)
b_lesser(posterior_samples(all_500_2)$b_Froh_500_s) #new one pr(b) - drops from 98.3 to 95.6, this is fine


#___1500 ----
loo(all_1500) #all good

#___5000 ----
loo(all_5000) #one ok value
loo(all_5000) %>% pareto_k_ids(threshold=0.5) #which one?
d %>% slice(loo::loo(all_5000) %>% pareto_k_ids(threshold=0.5))  #this one
pareto_k_values(loo(all_5000))[loo(all_5000) %>% pareto_k_ids(threshold=0.5)] #how bad?

#rerun without this point and check if predictor is still meaningful
all_5000_2 <- update(all_5000, newdata = d %>% filter(IID != "HMB032"))

b_lesser(posterior_samples(all_5000)$b_Froh_5000_s) #old pr(b)
b_lesser(posterior_samples(all_5000_2)$b_Froh_5000_s) #new one pr(b) - stays the same


###############################################################################
##### African Pops: Compare ROH 
Africans_ROH <- read.table("~/Desktop/Africans_ROH.txt", header=T)
# reorder for plot
Africans_ROH$FID <- factor(Africans_ROH$FID , levels=c("BENCH", "CHABU", "HIMBA", "ZULU", "IGBO", "MANDINKA"))

roh_comp <- ggplot(Africans_ROH, aes(x=FID, y=KB_500, fill=FID)) +
  geom_violin(show.legend=FALSE) + 
  geom_boxplot(width=0.1, show.legend=FALSE) +
  scale_y_continuous(name="Total ROH 500+ KB", limits=c(0,310000)) +
  xlab("Population") + theme(axis.text=element_text(size=7)) +
  labs(tag = "B")
roh_comp

##### PLOT IBDNe #####
himba_unrel <- read.delim("~/Desktop/Himba_H3Africa_unrel4th_GERMLINE2_IBDNe_minibd-4.ne", header=TRUE)
library(dplyr)
himba_unrel_gen5plus <- himba_unrel %>%
  filter(GEN >= 5)

ibdne <- ggplot(himba_unrel_gen5plus, aes(x=GEN, y=NE)) + geom_line(size=1.5) +
  ylab("Effective Population Size") + xlab("Generations Ago") +
  geom_line(aes(x=GEN, y=LWR.95.CI)) +
  geom_line(aes(x=GEN, y=UPR.95.CI)) +
  geom_ribbon(aes(ymin=LWR.95.CI, ymax=UPR.95.CI), fill="blue", alpha=0.5) +
  labs(tag = "A")

##### Plot 2 panel IBDNe & African ROH comparison #####
grid.arrange(ibdne, roh_comp, ncol=2)

################################################################################
##### Plot HapLD
himba_hapld <- read.csv("~/Desktop/himba_hapld_bottleneck.csv", header=TRUE)
zulu_hapld <- read.csv("~/Desktop/zulu_hapld_bottleneck.csv", header=TRUE)

himbaHapLD.plt <- ggplot(himba_hapld, aes(x=TIME, y=Q0.5)) + geom_line(size=1.5) +
  ylab("Effective Population Size") + xlab("Generations Ago") +
  geom_line(aes(x=TIME, y=Q0.975)) +
  geom_line(aes(x=TIME, y=Q0.025)) +
  geom_ribbon(aes(ymin=Q0.025, ymax=Q0.975), fill="blue", alpha=0.5) +
  labs(tag = "C")

temp_hapld_himba <- himba_hapld %>% mutate(Population = "Himba")
temp_hapld_zulu <- zulu_hapld %>% mutate(Population = "Zulu")
hapld_HimbaZulu_long <- rbind(temp_hapld_himba, temp_hapld_zulu)

colVector <- c(rep("#8494FF", 100), rep("#F8766D",100))
zuluHapLD.plt <- ggplot(hapld_HimbaZulu_long, aes(x=TIME, y=Q0.5, color=Population)) + 
  geom_line(size=1.5) + scale_color_manual(values=c("#8494FF", "#F8766D")) +
  ylab("Effective Population Size") + xlab("Generations Ago") +
  geom_line(aes(x=TIME, y=Q0.975)) +
  geom_line(aes(x=TIME, y=Q0.025)) +
  geom_ribbon(aes(ymin=Q0.025, ymax=Q0.975), fill=colVector, alpha=0.5) +
  labs(tag = "D") +
  theme(legend.key.size = unit(0.3, 'cm'))

grid.arrange(ibdne, roh_comp, himbaHapLD.plt, zuluHapLD.plt)

################################################################################
### Plot Simulated Bottleneck IBDNe ###

sim2_gen6_hapibd <- read.delim("~/Desktop/decliningNe_Ne1-4000_Ne2-450_n240_gen6_updated_noDups_hap-ibd_IBDNe_minibd-4_v20.ne", header=TRUE)
sim2_gen60_hapibd <- read.delim("~/Desktop/decliningNe_Ne1-4000_Ne2-450_n240_gen60_updated_noDups_hap-ibd_IBDNe_minibd-4_v20.ne", header=TRUE)

# function for plotting
plot.Ne.traces <- function(Data, CI.col=rgb(140/255, 155/255, 155/255, 0.7), Xlim, Ylim, ...) {
  args <- list(...)
  if(!is.null(args[["xlim"]]) | !is.null(args[["ylim"]])) {
    stop("Use 'Xlim' and 'Ylim' (instead of 'xlim' and 'ylim') to define axes limits")
  }
  
  if (Data[1,1]=="GEN") {
    stop("First line of data is header - re-import data with header = T")
  }
  
  # Initialize plot
  if (missing(Xlim)) {
    Xlim <- c(min(Data$GEN), max(Data$GEN))
  }
  
  if (missing(Ylim)) {
    Ylim <- c(1, max(c(Data$NE, Data$LWR.95.CI, Data$UPR.95.CI)))
  }
  
  plot(1, type="n", xlim=Xlim, ylim=Ylim, ...)
  
  # Draw CI of trace
  polygon(c(Data$GEN,rev(Data$GEN)),
          c(Data$LWR.95.CI,rev(Data$UPR.95.CI)),
          col=CI.col, border=NA)
  
  # Plot trace
  par(new=T)
  plot(Data$GEN, Data$NE, xlim=Xlim, ylim=Ylim, ...)
  par(new=F)
}

Ne <- c()
Gen <- c()
for (t in 0:6) {
  pop <- 4000*(2.71828^(t*-0.364))
  Ne <- c(Ne, pop)
  ga <- 6-t
  Gen <- c(Gen, ga)
}
decline6 <- as.data.frame(cbind(Gen, Ne))

Ne <- c()
Gen <- c()
for (t in 0:60) {
  pop <- 4000*(2.71828^(t*-0.0364))
  Ne <- c(Ne, pop)
  ga <- 60-t
  Gen <- c(Gen, ga)
}
decline60 <- as.data.frame(cbind(Gen, Ne))

par(mfrow=c(1,2))
plot.Ne.traces(sim2_gen6_hapibd,
               CI.col="gray",
               xlab="Generations", 
               ylab=expression('N'[e]),
               Xlim = c(0,100),
               Ylim = c(0,10000),
               pch=20,
               col="gray",
               main="Bottleneck: 6 ga")
lines(sim2_gen6_hapibd$GEN, sim2_gen6_hapibd$NE, lwd=4)
segments(100, 4000, 6, 4000, col="red", lty="dotted", lwd=2)
lines(decline6$Gen, decline6$Ne, col="red", lty="dotted", lwd=2)

plot.Ne.traces(sim2_gen60_hapibd,
               CI.col="gray",
               xlab="Generations", 
               ylab=expression('N'[e]),
               Xlim = c(0,100),
               Ylim = c(0,10000),
               pch=20,
               col="gray",
               main="Bottleneck: 60 ga")
lines(sim2_gen60_hapibd$GEN, sim2_gen60_hapibd$NE, lwd=4)
segments(100, 4000, 60, 4000, col="red", lty="dotted", lwd=2)
lines(decline60$Gen, decline60$Ne, col="red", lty="dotted", lwd=2)

# return to normal single plot
par(mfrow=c(1,1))

################################################################################
### Plot Froh for Simulated and Actual Individuals ###

unrel_1500 <- read.table("~/Desktop/Sim_Froh/Himba_unrel4th_snp50_missing2_het1_kb1500.hom.indiv", header=T)
unrel_Froh1500 <- unrel_1500 %>%
  mutate(Froh_1500 = (KB*1000)/2943912336)

gen60_kb1500 <- read.table("~/Desktop/Sim_Froh/sim2_gen60_snp50_missing2_het1_kb1500.hom.indiv", header=TRUE)
froh1500_gen60_sim2 <- gen60_kb1500 %>%
  mutate(Froh_1500 = (KB*1000)/2789629597)
gen6_kb1500 <- read.table("~/Desktop/Sim_Froh/sim2_gen6_snp50_missing2_het1_kb1500.hom.indiv", header=TRUE)
froh1500_gen6_sim2 <- gen6_kb1500 %>%
  mutate(Froh_1500 = (KB*1000)/2789648850)

FID_6 <- rep("Simulated Gen 6", 120)
gen6_froh500_temp2 <- froh500_gen6_sim2 %>%
  select(IID, Froh_500)
Sim2_Gen6_Froh <- cbind(FID_6, gen6_froh500_temp2, froh1500_gen6_sim2$Froh_1500, froh5000_gen6_sim2$Froh_5000)
head(Sim2_Gen6_Froh)
colnames(Sim2_Gen6_Froh) <- c("FID", "IID", "Froh_500", "Froh_1500", "Froh_5000")

FID_60 <- rep("Simulated Gen 60", 120)
gen60_froh500_temp2 <- froh500_gen60_sim2 %>%
  select(IID, Froh_500)
Sim2_Gen60_Froh <- cbind(FID_60, gen60_froh500_temp2, froh1500_gen60_sim2$Froh_1500, froh5000_gen60_sim2$Froh_5000)
head(Sim2_Gen60_Froh)
colnames(Sim2_Gen60_Froh) <- c("FID", "IID", "Froh_500", "Froh_1500", "Froh_5000")

FID_unrel <- rep("unrelated Himba", 120)
head(unrel_Froh500)
Himba_unrel4th_Froh_temp <- unrel_Froh500 %>%
  select(IID, Froh_500)
Himba_unrel4th_Froh <- cbind(FID_unrel, Himba_unrel4th_Froh_temp, unrel_Froh1500$Froh_1500, unrel_Froh5000$Froh_5000)
colnames(Himba_unrel4th_Froh) <- c("FID", "IID", "Froh_500", "Froh_1500", "Froh_5000")
# Combine all 3
Froh_sim2_comp <- rbind(Sim2_Gen6_Froh, Sim2_Gen60_Froh, Himba_unrel4th_Froh)
# plot
ggplot(Froh_sim2_comp, aes(x=FID, y=Froh_1500, fill=FID)) +
  geom_violin() + 
  geom_boxplot(width=0.025) +
  scale_y_continuous(name="FROH 1500+ KB") +
  xlab("") +
  guides(fill=guide_legend(title="Dataset"))

# sim gen 6
range(froh1500_gen6_sim2$Froh_1500) 
mean(froh1500_gen6_sim2$Froh_1500) 
var(froh1500_gen6_sim2$Froh_1500) 
# sim gen 60
range(froh1500_gen60_sim2$Froh_1500) 
mean(froh1500_gen60_sim2$Froh_1500)
var(froh1500_gen60_sim2$Froh_1500) 
# Unnrel Himba
range(unrel_Froh1500$Froh_1500) 
mean(unrel_Froh1500$Froh_1500) 
var(unrel_Froh1500$Froh_1500)

wilcox.test(unrel_Froh1500$Froh_1500, froh1500_gen60_sim2$Froh_1500)
wilcox.test(unrel_Froh1500$Froh_1500, froh1500_gen6_sim2$Froh_1500) 

### Himba FIS vs FROH ###

library(dplyr)
library(ggplot2)
Himba_auto_Froh <- read.table("~/Desktop/Himba_Froh_files/Himba_auto_froh.txt", header=TRUE)
fis <- read.table("~/Desktop/Himba_Froh_files/FIS_megaex_h3africa.het", header=TRUE)

froh_fis <- Himba_auto_Froh
froh_fis$FIS = fis$F[match(froh_fis$IID, fis$IID)]
head(froh_fis)
froh_fis$FIS <- as.numeric(froh_fis$FIS)

meanpt <- data.frame(avg.Froh_1500 = mean(froh_fis$Froh_1500), avg.FIS = mean(froh_fis$FIS))
ggplot(froh_fis, aes(x=Froh_1500, y=FIS)) + 
  geom_point(fill="#0072B2", shape = 21, alpha = 0.5) + xlab("FROH 1500") +
  geom_abline(slope=1, intercept = 0, linetype="dashed") +
  geom_point(data=meanpt, mapping=aes(x = avg.Froh_1500, y = avg.FIS), col="red")

cor(froh_fis$Froh_1500, froh_fis$FIS) # 0.9188784 ~ 0.92
# The base R function cor() takes a matrix or data.frame and computes the correlation between all the column pairs.

froh_fis %>% filter(FIS < 0) %>% nrow() # n=457 --> 67% of people
t.test(froh_fis$FIS) # not sig dif from 0, p=1.343e-10


######################################################################################
##### IBD-sharing between partners
omoka_wgermIBD <- read.table("~/Desktop/UC Davis/HENN LAB/Himba Project/F_ROH/Omoka_Froh/omoka_Froh1500_wgermIBD.txt", header=TRUE)
non_omoka_relatedParents_IBD <- read.table("~/Desktop/UC Davis/HENN LAB/Himba Project/F_ROH/Omoka_Froh/non_omoka_relatedParents_IBD.txt", header=TRUE)
non_omoka_unrelParents_IBD <- read.table("~/Desktop/UC Davis/HENN LAB/Himba Project/F_ROH/Omoka_Froh/non_omoka_unrelParents_IBD.txt", header=TRUE)

wilcox.test(non_omoka_unrelParents_IBD$IBD_SF_Mom, non_omoka_relatedParents_IBD$IBD_SF_Mom) 
wilcox.test(omoka_wgermIBD$IBD_BioF_Mom, non_omoka_unrelParents_IBD$IBD_SF_Mom)
wilcox.test(omoka_wgermIBD$IBD_BioF_Mom, non_omoka_relatedParents_IBD$IBD_SF_Mom) 




