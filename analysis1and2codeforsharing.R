##Does social relationship quality predict coalition formation among adult female chimpanzees? 
#Fox et al. 2022

# Load Tools & Libraries ####

library(tidyverse)
library(lme4)
library(DHARMa)
library(effects)
library(ggeffects)
library(MuMIn)
library(gridExtra)
library(MuMIn)
library(grid)
library(cowplot)
library(emmeans)
library(progress)

ymd <- lubridate::ymd
mdy <- lubridate::mdy
month <- lubridate::month
year <- lubridate::year
ymd_hms <- lubridate::ymd_hms
hms <- lubridate::hms

select <- dplyr::select
filter <- dplyr::filter

contr <- glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun=2e5))

## vif.mer function - this is for testing vif within glmer models, is adaptation from rms::vif
## available here: https://osf.io/p6ahy/?pid=ezkpa 
vif.mer <- function (fit) {
  ## adapted from rms::vif
  
  v <- vcov(fit)
  nam <- names(fixef(fit))
  
  ## exclude intercepts
  ns <- sum(1 * (nam == "Intercept" | nam == "(Intercept)"))
  if (ns > 0) {
    v <- v[-(1:ns), -(1:ns), drop = FALSE]
    nam <- nam[-(1:ns)]
  }
  
  d <- diag(v)^0.5
  v <- diag(solve(v/(d %o% d)))
  names(v) <- nam
  v
}




# ANALYSIS 1 ####
# Load Data ####
analysis1_datashare <- read.csv("analysis1_datashare.csv", header = TRUE)

# Variable key ####

# ID1 & ID2: females 
# period: A-E
# partycount_d: number of scans dyad was in the same party
# coal_count_bi: number of coalitions dyad participated in during the two year period
# SRI_bi_z: party association index z-score
# index_5m_bi_z: five-meter index z-score
# gmdur_index_z: grooming duration index z-score
# gm_yn_bi: N = did not groom during period, Y = groomed during period
# ranksumbi_z: summed dyad dominance rank during period (sum of each individual's average during the period)
# immdyad: dyad type regarding female residency, 0 = resident-resident, 1 = immigrant-resident, 2 = immigrant-immigrant
# kinship: 0 = unrelated, 1 = kin

# Data Subsets ####
Cdata <- analysis1_datashare %>% #full data, 811 rows
  filter(partycount_d > 200) %>% #must have > 200 party scans together, 713 rows
  filter(total_AB_party > 0) %>% #must be in party together during a focal, 703 rows
  filter(!is.na(ranksumbi_z)) #must have dominance data for both IDs, rows 689
#to run models on data subset without kin dyads, use this: 
Cdata.kin <- Cdata %>% filter(kinship != "1")
#to run models on data subset with only resident females, use this:
Cdata.imm <- Cdata %>% filter(immdyad == "0")


# Example null model ####

C.mod0 <-
  glmer(coal_count_bi ~ 1 
        + ranksumbi_z
        + kinship
        + offset(log(partycount_d))
        + (1|ID1)
        + (1|ID2)
        + (1|period),
        nAGQ = 0,
        family = poisson,
        control = contr,
        data = Cdata)
summary(C.mod0)
testDispersion(C.mod0)
plotResiduals(simulateResiduals(C.mod0))
testResiduals(C.mod0)

# Example full model ####

C.mod15 <- 
  glmer(coal_count_bi ~ 1 
        + ranksumbi_z
        + gmdur_index_z
        + gm_yn_bi
        + index_5m_bi_z
        + SRI_bi_z
        + kinship
        + offset(log(partycount_d))
        + (1|ID1)
        + (1|ID2)
        + (1|period),
        nAGQ = 0,
        family = poisson,
        control=contr,
        data = Cdata)
summary(C.mod15)
plotResiduals(simulationOutput= (simulateResiduals(C.mod15)))
testResiduals(C.mod15)
vif.mer(C.mod15)

# Example model selection ####
#we used dredge to run the model comparison but we also ran each model individually to test model assumptions using code similar to above for each model  
options(na.action = "na.fail")
#if you want to run this in the dataset without kin, remove kinship here
C15dredge <- dredge(C.mod15, fixed = c("ranksumbi_z", "offset(log(partycount_d))", "kinship"))
C15.95 <- subset(C15dredge, cumsum(weight) <= .95, recalc.weights = FALSE)
C15modavg <- model.avg(subset(C15dredge, cumsum(weight) <= .95, recalc.weights = FALSE))
summary(C15modavg)
get.models(C15.95, subset = TRUE)

#r.squaredGLMM is from MuMin package, see Nakagawa et al. 2017 http://doi.org/10.1098/rsif.2017.0213 
r.squaredGLMM(C.mod15)

# Example graphs ####
#PARTY ASSOC GRAPH
set.seed(123)
predict1 <- ggemmeans(C.mod15, 
                      terms = c("SRI_bi_z"), 
                      #type = "zero_inflated", 
                      back.transform = TRUE)
p1 <- ggplot(predict1, aes(x, predicted)) +
  theme_bw(base_size = 18) + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        legend.position = "none") +
  scale_y_continuous(breaks = c(0, 2, 4, 6, 8, 10), trans = "log1p") +
  geom_jitter(data = Cdata,
              aes(x = SRI_bi_z, y = coal_count_bi), 
              height = 0.05, color = "grey64", alpha = 0.5, size = 2) +
  xlab("Party association index z-score") + ylab("Number of Coalitions") +
  geom_line(colour = "purple") + 
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, alpha = 0.1), fill = "purple")
print(p1)

#FIVE M GRAPH
set.seed(123)
predict2 <- ggemmeans(C.mod15, 
                      terms = c("index_5m_bi_z"), 
                      #type = "random", 
                      back.transform = TRUE)
p2 <- ggplot(predict2, aes(x, predicted)) +
  theme_bw(base_size = 18) + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        legend.position = "none") +
  scale_y_continuous(breaks = c(0, 2, 4, 6, 8, 10), trans = "log1p") + 
  geom_jitter(data = Cdata,
              aes(x = index_5m_bi_z, y = coal_count_bi), 
              height = 0.05, color = "grey64", alpha = 0.5, size = 2) + 
  xlab("Five meter association index z-score") + ylab("") +
  geom_line(colour = "purple") + 
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, alpha = 0.1), fill = "purple")
print(p2)

#GROOMING GRAPH *need to make sure it is adding actual data not predicted data?
predict3 <- ggemmeans(C.mod15, 
                      terms = c("gmdur_index_z[all]"), 
                      back.transform = TRUE)  

p3 <- ggplot(predict3, aes(x, predicted)) +
  theme_bw(base_size = 18) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none") +
  scale_y_continuous(breaks = c(0, 2, 4, 6, 8, 10), trans = "log1p") +
  scale_x_continuous(limits = c(-0.5, 16)) +
  geom_jitter(data = Cdata, 
              aes(x = gmdur_index_z, y = coal_count_bi),
              height = 0.0, width = 0, color = "grey64", alpha = 0.5, size = 2) +
  xlab("Grooming duration index z-score") + 
  ylab("Number of Coalitions") +
  geom_line(colour = "purple") +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, alpha = 0.1), fill = "purple")
print(p3)

#make p3a to make an insert
p3a <- ggplot(predict3, aes(x, predicted)) +
  theme_bw(base_size = 18) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none") +
  scale_y_continuous(breaks = c(0, 2, 10), trans = "log1p") +
  scale_x_continuous(limits = c(-.5, 1.5), breaks = c(0, 1)) +
  geom_jitter(data = Cdata, 
              aes(x = gmdur_index_z, y = coal_count_bi),
              height = 0.0, width = 0, color = "grey64", alpha = 0.5, size = 2) +
  xlab(" ") + ylab(" ") +
  geom_line(colour = "purple") +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, alpha = 0.1), fill = "purple")
print(p3a)

p3_inset <- p3 + annotation_custom(ggplotGrob(p3a), xmin = 5, xmax = 16, ymin = 1.08)
print(p3_inset)

predict4 <- ggemmeans(C.mod15,
                      terms = "gm_yn_bi",
                      back.transform = FALSE)
pred4 <- predict4 %>% data.frame()
plot(predict4)
#Plot the real data with mean and CI manually added
p4 <- ggplot() +
  geom_boxplot(data = Cdata,
               mapping = aes(gm_yn_bi, y = coal_count_bi),
               outlier.shape = NA,
               colour = "grey") +
  geom_jitter(data = Cdata,
              mapping = aes(gm_yn_bi, y = coal_count_bi),
              color = "gray64", height = 0.1, width = 0.3, alpha = 0.5, size = 2) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        plot.title = element_blank(),
        text = element_text(size = 18)) +
  xlab("Groomed") + ylab("") +
  scale_y_continuous(breaks= c(0, 2, 4, 6, 8, 10), trans = "log1p", limits = c(-0.1, NA)) +
  geom_pointrange(data = pred4, mapping = aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high),
                  colour = "purple", size = 1, shape = "diamond")

print(p4)
plot_grid(p1, p2, p3_inset, p4, nrow = 2, labels = c("A", "B", "C", "D"))


# Example permutation ####
# Main data
# Isolate IDs
# Sample IDs and turn back into df
# Rejoin them to data
# 
# Create Storage for each parameter (empty vector, df, or list)
# Run model on data
# Extract coefficient, statistic, and p value and save to storage
# Do this 1000 times
# 
# Summarise output as averages of each column and proportion of pvals that are signif
# Note that results may shift slightly because IDs were shuffled once for starting dataset provided here (ie. 1 permutation different from paper)
# Example for Model 0 - ie the NULL Model
set.seed(5)
list.coef <- vector("list", length = 1000)
list.AIC <- vector("list", length = 1000)
full_data_noID <- Cdata %>% select(-ID1, -ID2)

pb <- progress_bar$new(format = " running [:bar] :percent eta: :eta",
                       total = 1000, clear = FALSE, width= 60)
t <- Sys.time()
for(i in seq(1000)){
  pb$tick()
  Sys.sleep(1/1000)
  
  IDs <- Cdata %>% 
    select(ID1, ID2) %>%
    apply(., 1, sample) %>% t() %>% data.frame() %>%
    rename(ID1 = X1, ID2 = X2)
  
  full_data_rID <- full_data_noID %>%
    cbind(., IDs) 
  
  C.mod0 <-
    glmer(coal_count_bi ~ 1 
          + ranksumbi_z
          + kinship
          + offset(log(partycount_d))
          + (1|ID1)
          + (1|ID2)
          + (1|period),
          nAGQ = 0,
          family = "poisson",
          control = contr,
          data = full_data_rID)
  
  list.coef[[i]] <- coef(summary(C.mod0)) %>% 
    data.frame() %>% 
    rownames_to_column()
  
  list.AIC[[i]] <- llikAIC(C.mod0) %>% 
    data.frame() %>% 
    rownames_to_column()
}
Sys.time() - t
r.coefs.0 <- do.call("rbind", list.coef)
r.AIC.0 <- do.call("rbind", list.AIC) %>% filter(rowname == "AIC") %>% summarise(n = n(), AIC_mean = mean(AICtab))
#make these into a table 
coef.avgd.0 <- r.coefs.0 %>%
  group_by(rowname) %>%
  summarise(n = n(), 
            Est_mean = mean(Estimate), 
            SE_mean = mean(Std..Error), 
            z_mean = mean(z.value),
            pval_prop = length(which(Pr...z.. < 0.05))/n*100) %>% 
  ungroup() %>%
  mutate(model_num = "model_0")
d1 <- coef.avgd.0 %>% filter(rowname == "(Intercept)") %>% 
  rename_with(~str_c("Int_", .), .cols = c(Est_mean, SE_mean, z_mean, pval_prop)) %>%
  select(-rowname)
d2 <- coef.avgd.0 %>% filter(rowname == "ranksumbi_z")%>% 
  rename_with(~str_c("Rank_", .), .cols = c(Est_mean, SE_mean, z_mean, pval_prop)) %>%
  select(-rowname)
d3 <- coef.avgd.0 %>% filter(rowname == "kinship1")%>% 
  rename_with(~str_c("Kin_", .), .cols = c(Est_mean, SE_mean, z_mean, pval_prop)) %>%
  select(-rowname)
mod0_results <- d1 %>%
  left_join(d2) %>%
  left_join(d3) %>%
  left_join(r.AIC.0)
#View(mod0_results)

# Permutation of model averaging ####

# Run dredge 1000 times and take model averaged results
# Each iteration, we remove the IDs, randomly shuffle positions of ID1 and ID2, then reattach  
# This took my computer about 40 minutes

list.coef <- vector("list", length = 1000)
list.msTable <- vector("list", length = 1000)
list.dredge <- vector("list", length = 1000)
options(na.action = "na.fail")
full_data_noID <- Cdata %>% select(-ID1, -ID2)
set.seed(5)
t <- Sys.time()
for(i in seq(1000)){
  Sys.sleep(1/1000)
  
  IDs <- Cdata %>% 
    select(ID1, ID2) %>%
    apply(., 1, sample) %>% t() %>% data.frame() %>%
    rename(ID1 = X1, ID2 = X2)
  
  full_data_rID <- full_data_noID %>%
    cbind(., IDs) 
  
  C.mod15 <- 
    glmer(coal_count_bi ~ 1 
          + ranksumbi_z
          + gmdur_index_z
          + gm_yn_bi
          + index_5m_bi_z
          + SRI_bi_z
          + kinship
          + offset(log(partycount_d))
          + (1|ID1)
          + (1|ID2)
          + (1|period),
          nAGQ = 0,
          control=contr,
          family = "poisson",
          data = full_data_rID)
  
  C15dredge <- dredge(C.mod15, fixed = c("ranksumbi_z", "offset(log(partycount_d))", "kinship"))
  C15.95 <- subset(C15dredge, cumsum(weight) <= .95, recalc.weights = FALSE)
  #in run 442 only one model is in the 95% subset, so it can't be averaged and causes an error
  #so we just extract the information on run 442 (don't add it to the list)
  # below we extract the information from it outside the loop and add it back in
  if(i == 442) {
    data442 <- full_data_rID
    dredge442 <- C15dredge
    C15.95442 <- C15.95
  } 
  else{
    C15modavg <- model.avg(subset(C15dredge, cumsum(weight) <= .95, recalc.weights = FALSE))
    
    list.coef[[i]] <- summary(C15modavg)$coefmat.subset %>% 
      data.frame() %>%
      mutate(run = i) %>%
      rownames_to_column()
    
    #don't use these but saving just in case  
    list.msTable[[i]] <- C15modavg$msTable %>% 
      data.frame() %>%
      mutate(run = i) %>%
      rownames_to_column(var = "modelcontents")
    
    list.dredge[[i]] <- C15dredge %>%
      data.frame() %>%
      mutate(run = i) %>%
      rownames_to_column(var = "model_num")
  }
}
Sys.time() - t
coefs.modavg <- do.call("rbind", list.coef)
#add back in the results from model 441 
msTable.modavg <- do.call("rbind", list.msTable)
dredge.modavg <- do.call("rbind", list.dredge)


# save(data442,
#      dredge442,
#      C15.95442,
#      coefs.modavg,
#      msTable.modavg,
#      dredge.modavg,
#      list.coef,
#      list.msTable,
#      list.dredge,
#      file = "avg of the avg.Rdata")


#extract what you want from run 442
dredge442a <- dredge442 %>%
  rownames_to_column(var = "model_num") %>%
  mutate(run = 442) %>%
  rename(X.Intercept. = "(Intercept)", 
         offset.log.partycount_d.. = "offset(log(partycount_d))")

# View(C15.95442)
#the top model of run 442 contained gmyn, 5m, SRI 
#We need ot re run that model with the 442 data to get the ceof we want to add to the coef table
C.mod15.442 <- 
  glmer(coal_count_bi ~ 1 
        + ranksumbi_z
        + gm_yn_bi
        + index_5m_bi_z
        + SRI_bi_z
        + kinship
        + offset(log(partycount_d))
        + (1|ID1)
        + (1|ID2)
        + (1|period),
        nAGQ = 0,
        control=contr,
        family = "poisson",
        data = data442)
coefs.442 <- summary(C.mod15.442)$coef %>%
  data.frame() %>%
  mutate(run = 442) %>%
  rownames_to_column() 


#summarise these results: coefs from top models and AIC/weights from all models

#we want to average these coef to get the average B/se/pval *from the models in the 95% subset*
coefs.modavg.summary <- coefs.modavg %>%
  select(-Adjusted.SE) %>%
  #add on the data from 442
  rbind(coefs.442) %>%
  group_by(rowname) %>%
  summarise(n = n(), #just want to check n to make sure we are getting all the data
            Est_mean = mean(Estimate), 
            SE_mean = mean(Std..Error), 
            z_mean = mean(z.value),
            pval_prop = length(which(Pr...z.. < 0.05))/n*100) %>% 
  ungroup() %>%
  column_to_rownames(., var = "rowname") %>%
  select(-z_mean, -n) %>%
  round(., digits = 3) %>%
  rotate_df()

write.csv(coefs.modavg.summary, file = "coefsmodavgsummary.csv")

#we want to average the AIC and weights from each run of each model
dredge.modavg.summary <- dredge.modavg %>%
  #add in model 442 info
  rbind(dredge442a) %>%
  group_by(model_num) %>%
  summarise(AICc_mean = mean(AICc), 
            df = mean(df), 
            delta = mean(delta), 
            weight = mean(weight)) %>%
  column_to_rownames(var = "model_num") %>%
  round(., digits = 3) %>%
  arrange(AICc_mean) %>%
  rotate_df()

write.csv(dredge.modavg.summary, file = "dredgemodavgsummary.csv")


# ANALYSIS 2 ####

# Load Data ####
analysis2_datashare <- read.csv("analysis2_datashare.csv", header = TRUE)

# Variable key ####

# ID1: female who was part of coalitions
# party_member: potential coalition partner (present in party at time of coalition)
# coal_id: unique coalition event identified
# event_factor: 1 = partner was chosen for coalition, 0 = partner was not chosen 
# SRI_bi_z: party association index z-score
# index_5m_bi_z: five-meter index z-score
# gmdur_index_z: grooming duration index z-score
# gm_yn_bi: N = did not groom during period, Y = groomed during period
# ranksumbi_z: summed dyad dominance rank during period (sum of each individual's average during the period)
# immdyad: dyad type regarding female residency, 0 = resident-resident, 1 = immigrant-resident, 2 = immigrant-immigrant
# kinship: 0 = unrelated, 1 = kin
# period: A-E
# party_size_fem_z: z-score number of females in the party (in addition to ID1, so really n+1 = # females in the party)

# Data subsets ####
#to run without kin dyads use this
log15kin <- analysis2_datashare %>% filter(kinship != "1")
#to run with only resident females use this
log15imm <- analysis2_datashare %>% filter(immdyad == "0")

# Example null model ####
log.opp0 <- glmer(event_factor ~ 1 
                  + party_size_fem_z
                  + ranksumbi_z
                  + kinship
                  + (1|coal_id)
                  + (1|party_member)
                  + (1|period), 
                  family = binomial(logit), 
                  data = analysis2_datashare,
                  control = contr)
summary(log.opp0)
testDispersion(log.opp0)
plot(simulateResiduals(log.opp0))
vif.mer(log.opp0)

# Example full model ####
log.opp15 <- glmer(event_factor ~ 1 
                   + index_5m_bi_z
                   + SRI_bi_z
                   + gmdur_index_z
                   + gm_yn_bi
                   + kinship
                   + ranksumbi_z
                   + party_size_fem_z 
                   + (1|coal_id)
                   + (1|party_member)
                   + (1|period), 
                   family = binomial(logit), 
                   data = analysis2_datashare,
                   control = contr)
summary(log.opp15)
testDispersion(log.opp15)
plot(simulateResiduals(log.opp15))
vif.mer(log.opp15)

# Example model selection ####
#we used dredge to run the model comparison but we also ran each model individually to test model assumptions using code similar to above for each model  
options(na.action = "na.fail")
#if you want to run this in the dataset without kin, remove kinship here
log15dredge <- dredge(log.opp15, fixed = c("party_size_fem_z", "ranksumbi_z", "kinship"))
log15.95 <- subset(log15dredge, cumsum(weight) <= .95, recalc.weights = FALSE)
log15modavg <- model.avg(subset(log15dredge, cumsum(weight) <= .95, recalc.weights = FALSE))
summary(log15modavg)
get.models(log15.95, subset = TRUE)
r.squaredGLMM(log.opp15)

# Example Graphs #### 

#PARTY ASSOC GRAPH
predict5 <- ggeffect(log.opp15, terms = "SRI_bi_z[all]", back.transform = TRUE)
p5 <- ggplot(predict5, aes(x, predicted)) +
  geom_line() + 
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, alpha = 0.1), fill ="purple") +
  theme_bw(base_size = 18) + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        legend.position = "none") +
  geom_point(data = analysis2_datashare, aes(x = SRI_bi_z, y = event_factor), alpha =0.2, size = 2, colour = "grey64") +
  xlab("Party association index z-score") + ylab("Probability of \n partner selection") 
print(p5)

#FIVE M GRAPH
predict6 <- ggeffect(log.opp15, terms = "index_5m_bi_z[all]", back.transform = TRUE)
p6 <- ggplot(predict6, aes(x, predicted)) +
  geom_line() + 
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, alpha = 0.1, size = 2), fill = "purple") +
  theme_bw(base_size = 18) + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        legend.position = "none") +
  geom_point(data = analysis2_datashare, aes(x = index_5m_bi_z, y = event_factor), alpha = 0.2, colour = "grey64") +
  xlab("Five meter association index z-score") + ylab("") 
print(p6)

#GROOMING GRAPHS
predict7 <- ggeffect(log.opp15, terms = "gmdur_index_z[all]", back.transform = TRUE)
p7 <- ggplot(predict7, aes(x, predicted)) +
  geom_line() + 
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, alpha = 0.1), fill = "purple") +
  theme_bw(base_size = 18) + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        legend.position = "none") +
  geom_point(data = analysis2_datashare, aes(x = gmdur_index_z, y = event_factor), alpha = 0.2, colour = "grey64") +
  xlab("Grooming duration index z-score") +
  ylab("Probability of \n partner selection") +
  scale_x_continuous(limits = c(-1, 12), breaks = c(-1, 0, 2, 4, 6, 8, 10, 12))
print(p7)

predict8 <- ggeffect(log.opp15, terms = "gm_yn_bi", back.transform = TRUE)
pred8 <- predict8 %>% data.frame()
p8 <- ggplot() +
  geom_pointrange(data = pred8, mapping = aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high), 
                  colour = "purple", size = 1, shape = "diamond") +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        legend.position = "none", 
        plot.title = element_blank(), 
        text = element_text(size = 18)) +
  xlab("Groomed") + ylab("")
print(p8)

plot_grid(p5, p6, p7, p8, nrow = 2, labels = c("A", "B", "C", "D"))
