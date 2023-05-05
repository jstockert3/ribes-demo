## Notes from Natali:
## Need confidence intervals
## elasticity analysis very important for conservation
## simulate some environmental variable like impact of pollination limitation
## sensitivity and elasticity both in popbio
## elasticity easier to interpret because they sum to 1
## life table response experiments
## need more years for stochastic simulations
## Look for Conradina glabra demography Natural Areas Journal


library(tidyverse)
library(popbio)
setwd("~/Documents/ribes")
source('./scripts/load_data.R')
sdct <- read_excel('./data/seed_count.xlsx')
seeds <- read_excel('./data/seedtrials-all.xlsx')

#### Some helper functions for data cleaning ####
##add binomial survival variable based stem count in 2023
add_survival_repr_factor <- function(x) {
  df <- x
  df <- df %>% mutate(surv = case_when(
    Sum_stems.23 > 0 ~ 1,
    Sum_stems.23 == 0 ~ 0
  ), repr = case_when(
    N_flowers.23 > 0 | N_flowers.22 > 0 ~ 1,
    N_flowers.23 == 0 & N_flowers.22 ==0 ~ 0,
  ), site = case_when(
    Plot == 'hp1a' | Plot == 'hp1b' ~ 'hp1',
    Plot == 'hp2a' | Plot == 'hp2b' ~ 'hp2',
    Plot == 'nf7a' | Plot == 'nf7b' ~ 'nf7',
    Plot == 'nf11a' | Plot == 'nf11b' ~ 'nf11',
    Plot == 'nna' | Plot == 'nnb' ~ 'nn',
    Plot == 'nsa' | Plot == 'nsb' ~ 'ns',
    Plot == 'mpna' | Plot == 'mpnb' ~ 'mpn',
    Plot == 'mpsa' | Plot == 'mpsb' ~ 'mps',
  ), repr1 = case_when(
    N_flowers.22 > 0 ~ 1,
    N_flowers.22 == 0 ~ 0
  ), repr2 = case_when(
    N_flowers.23 > 0 ~ 1,
    N_flowers.23 == 0 ~ 0
  ))
  df$Surv <- as.factor(df$surv)
  df$Population <- as.factor(df$Population)
  df$Treatment <- as.factor(df$Treatment)
  df$Plot <- as.factor(df$Plot)
  return(df)
}
#remove new individuals and missing tags DO NOT USE for fertility
rm_missing <- function(x){
  df <- x %>% filter(!is.na(surv))
  return(df)
}
rm_new <- function(x){
  df <- x %>% filter(!is.na(len))
}

#### dataset cleanup ####
lookup <- c(len = 'Len_LS_cm.22', len_next = 'Len_LS_cm.23', 
            flw = 'N_flowers.22', flw_next = 'N_flowers.23',
            woody = 'N_Wstem.22', woody_next = 'N_Wstem.23',
            green = 'N_Gstem.22', green_next = 'N_Gstem.23',
            stems = 'Sum_stems.22', stems_next = 'Sum_stems.23',
            h = 'Height_cm.22', h_next = 'Height_cm.23',
            w = 'Width_cm.22', w_next = 'Width_cm.23',
            vol = 'VOLUME.22', vol_next='VOLUME.23', Site='site')
ribes22_23 <- add_survival_repr_factor(ribes22_23) %>% 
  rename(all_of(lookup))
ribes22_23 <- ribes22_23 %>%
  select(Population, Treatment, Plot, TagNum, Site,
         flw, len, flw_next, len_next,
         woody, green, stems, h, w, vol,
         woody_next, green_next, stems_next, h_next, w_next, vol_next, surv,
         repr, repr1, repr2)
ribes22_23 <- rm_missing(ribes22_23)
ribes_noNew <- rm_new(ribes22_23)
ribes_noDead <- ribes_noNew %>% filter(surv != 0)
ribes_rec <- ribes22_23 %>% filter(is.na(len))
seeds <- seeds %>% mutate(notGerm=5-Ngerm)
ribes<-ribes22_23
ribes_sc <- ribes %>% filter(Population=='sc')
ribes_fl <- ribes %>% filter(Population=='fl')



#### Define Size Classes ####

# 1. Select state variable
cor.test(ribes$h, ribes$surv) # r = 0.21
cor.test(ribes$len, ribes$surv) # r = 0.18
cor.test(ribes$woody, ribes$surv) # r = 0.16
cor.test(ribes$green, ribes$surv) # r = 0.05
cor.test(ribes$w, ribes$surv) # r = 0.16
cor.test(ribes$vol, ribes$surv) # r = 0.09
cor.test(ribes$flw, ribes$surv) # NS, p = 0.14


# 2. Explore data for size cutoffs
ggplot(ribes, aes(x=h)) + 
  geom_histogram(fill='lightblue', color='black') +
  geom_histogram(data=ribes %>% filter(surv==1), aes(x=h),
                 fill='red', color='black') +
  geom_histogram(data=ribes %>% filter(repr1==1), aes(x=h),
                 fill='orange', color='black') +
  scale_x_continuous(breaks=seq(0,150,5))

# 3. Separate into size classes
## pre-reproductive: height from 0 cm to 22 cm
## small: from 22 cm to 45 cm
## medium: from 45 cm to 78 cm
## large: greater than 78 cm

ribes_manual_classes <- ribes %>%
  mutate(class=case_when(h < 22 ~ "pre",
                         h >= 22 & h <= 45 ~ "sm",
                         h > 45 & h <= 78 ~ 'md',
                         h > 78 ~ 'lg'))

ggplot(ribes_manual_classes, aes(x=h, fill=class)) + 
  geom_histogram(color='black') +
  scale_fill_brewer(palette = 'Dark2')

ggplot(ribes_manual_classes, aes(x=h, fill=class)) + 
  geom_histogram(alpha=0.5) +
  scale_fill_brewer(palette = 'Dark2') +
  geom_histogram(data=ribes_manual_classes %>% filter(surv==1),
                 aes(x=h, fill=class), alpha=0.5) +
  geom_histogram(data=ribes_manual_classes %>% filter(repr1==1),
                 aes(x=h, fill=class))

ribes_manual_classes <- ribes_manual_classes %>%
  mutate(class_next = case_when(h_next > 0 & h_next < 22 ~ "pre",
                              h_next >= 22 & h_next <= 45 ~ "sm",
                              h_next > 45 & h_next <= 78 ~ 'md',
                              h_next > 78 ~ 'lg',
                              h_next == 0 ~ 'dead'))


#### Create Transition Matrix ####
### South Carolina Population
totals_sc <- ribes_manual_classes %>% filter(!is.na(class)) %>%
  filter(Population=='sc') %>%
  summarise(pre=length(class[class=='pre']), 
            sm=length(class[class=='sm']),
            md=length(class[class=='md']), 
            lg=length(class[class=='lg'])) 

View(totals_sc)

transition_sc <- ribes_manual_classes %>%
  filter(!is.na(class)) %>%
  filter(Population=='sc') %>%
  group_by(class_next) %>%
  summarise(pre=length(class[class=='pre']) / totals_sc$pre, 
            sm=length(class[class=='sm']) / totals_sc$sm,
            md=length(class[class=='md']) / totals_sc$md, 
            lg=length(class[class=='lg']) / totals_sc$lg ) %>%
  arrange(factor(class_next, c('pre', 'sm', 'md', 'lg')))

View(transition_sc)
sum(transition_sc$pre)


### Florida Population
totals_fl <- ribes_manual_classes %>% filter(!is.na(class)) %>%
  filter(Population=='fl') %>%
  summarise(pre=length(class[class=='pre']), 
            sm=length(class[class=='sm']),
            md=length(class[class=='md']), 
            lg=length(class[class=='lg']))

transition_fl <- ribes_manual_classes %>%
  filter(!is.na(class)) %>%
  filter(Population=='fl') %>%
  group_by(class_next) %>%
  summarise(pre=length(class[class=='pre']) / totals_fl$pre, 
            sm=length(class[class=='sm']) / totals_fl$sm,
            md=length(class[class=='md']) / totals_fl$md, 
            lg=length(class[class=='lg']) / totals_fl$lg ) %>%
  arrange(factor(class_next, c('pre', 'sm', 'md', 'lg')))

View(transition_fl)
View(totals_fl)

#### Adjust size classes ####
ribes_fl <- ribes %>% filter(Population=='fl')

ggplot(ribes_fl, aes(x=h)) + 
  geom_histogram(fill='lightblue', color='black') +
  geom_histogram(data=ribes_fl %>% filter(surv==1), aes(x=h),
                 fill='red', color='black') +
  geom_histogram(data=ribes_fl %>% filter(repr1==1), aes(x=h),
                 fill='orange', color='black') +
  scale_x_continuous(breaks=seq(0,150,5))

ribes_manual_classes_fl <- ribes_fl %>%
  mutate(class=case_when(h > 0 & h < 23 ~ "pre",
                         h >= 23 & h <= 42 ~ "sm",
                         h > 42 ~ 'lg'))

ggplot(ribes_manual_classes_fl, aes(x=h, fill=class)) + 
  geom_histogram(color='black') +
  scale_fill_brewer(palette = 'Dark2')

ggplot(ribes_manual_classes_fl, aes(x=h, fill=class)) + 
  geom_histogram(alpha=0.5) +
  scale_fill_brewer(palette = 'Dark2') +
  geom_histogram(data=ribes_manual_classes_fl %>% filter(surv==1),
                 aes(x=h, fill=class), alpha=0.5) +
  geom_histogram(data=ribes_manual_classes_fl %>% filter(repr1==1),
                 aes(x=h, fill=class))

ribes_manual_classes_fl <- ribes_manual_classes_fl %>%
  mutate(class_next = case_when(h_next > 0 & h_next < 23 ~ "pre",
                                h_next >= 23 & h_next <= 42 ~ "sm",
                                h_next > 42 ~ 'lg',
                                h_next == 0 ~ 'dead'))

##Re-calculate transition matrix
totals_fl <- ribes_manual_classes_fl %>% 
  filter(!is.na(class)) %>%
  summarise(pre=length(class[class=='pre']), 
            sm=length(class[class=='sm']), 
            lg=length(class[class=='lg']))

transition_fl <- ribes_manual_classes_fl %>%
  filter(!is.na(class)) %>%
  group_by(class_next) %>%
  summarise(pre=length(class[class=='pre']) / totals_fl$pre, 
            sm=length(class[class=='sm']) / totals_fl$sm, 
            lg=length(class[class=='lg']) / totals_fl$lg ) %>%
  arrange(factor(class_next, c('pre', 'sm', 'lg', 'dead')))

View(transition_fl)
View(totals_fl)


#### Create fecundity matrix ####

####South Carolina Population####

#Pull out recruits from dataset by filtering for indiv with size 0 in year 1
#Get number of recruits in each size class in year 2, then divide
#by total number of recruits to get proportions
#Do some formatting to coerce proportions to a matrix-like structure
recruits_sc <- ribes_manual_classes %>%
  filter(is.na(class)) %>% 
  filter(Population=='sc') %>% 
  summarise(pre=length(class_next[class_next=='pre'])/length(class_next), 
            sm=length(class_next[class_next=='sm'])/length(class_next),
            md=length(class_next[class_next=='md'])/length(class_next), 
            lg=length(class_next[class_next=='lg'])/length(class_next)) %>%
  mutate(var='p_size') %>% 
  pivot_longer(cols=!var, names_to = 'class') %>%
  select(!var) %>% 
  mutate(pre=value, sm=value, md=value, lg=value) %>%
  select(!value) %>% 
  column_to_rownames(var='class')
View(recruits_sc)

#Get the total number of flowers produced by individuals in each size class
# in year 1, as well as number of flowering individuals
flowers_sc_1 <- ribes_manual_classes %>%
  filter(!is.na(class)) %>%
  filter(Population=='sc') %>%
  group_by(class) %>% 
  summarise(flowers = sum(flw),
            n=sum(repr1))

#Same for year 2
flowers_sc_2 <- ribes_manual_classes %>%
  filter(!is.na(class_next), class_next!='dead') %>%
  filter(Population=='sc') %>%
  group_by(class_next) %>% 
  summarise(flowers = sum(flw_next),
            n=sum(repr2)) %>%
  rename(class=class_next)

#Add flower sums together, and then take an average
flowers_sc <- rbind(flowers_sc_2, flowers_sc_1) %>% 
  group_by(class) %>% 
  summarise(flowers=sum(flowers), n=sum(n),
            avg = flowers/n) %>%
  select(class,avg) %>%
  arrange(factor(class, c('pre', 'sm', 'md', 'lg')))

#Turn NA values into zeros
flowers_sc$avg[is.nan(flowers_sc$avg)]<-0

#Pivot flower data
flowers_sc <- flowers_sc %>% 
  pivot_wider(names_from = class,
              values_from = avg)

#Get average average germination rate per seed
germ_sc <- seeds %>% #data from in situ germination trials
  filter(Population=='sc') %>% 
  summarise(germ=sum(Ngerm),
            not=sum(notGerm),
            mean=sum(Ngerm)/(sum(Ngerm)+sum(notGerm)))
germ_sc = germ_sc$mean

#Get average number of seeds per fruit
avseeds_sc <- sdct %>% 
  filter(Population=='sc') %>%
  summarise(n=length(NumSeeds),
            sum=sum(NumSeeds),
            avg=sum(NumSeeds)/length(NumSeeds))
avseeds_sc = avseeds_sc$avg

#Get total number of flowers
flower_total_sc <- ribes_manual_classes %>%
  filter(!is.na(class)) %>% 
  filter(Population=='sc') %>%
  summarise(flowers=sum(flw))
flower_total_sc = flower_total_sc$flowers

#Get total number of recruits
recruit_n_sc <- ribes_manual_classes %>% 
  filter(is.na(class)) %>% 
  filter(Population=='sc') %>%
  summarise(n=sum(surv))
recruit_n_sc = recruit_n_sc$n

#Calculate overall recruitment rate = 
# number of recruits / number of flowers * seeds * germination rate
recruit_rate_sc = recruit_n_sc / (flower_total_sc * avseeds_sc * germ_sc)

#Multiply average number of flowers produced in each size class by
#average seeds, germination rate, and recruitment rate
#to get the number of recruits produced in each size class
flowers_sc <- flowers_sc * avseeds_sc * germ_sc * recruit_rate_sc

#Multiply average number of recruits produced by size proportions
#to get a matrix of contributions to each size class in year 2 by
#flowering individuals in year 1
fecundity_sc <- data.frame(pre=recruits_sc$pre*flowers_sc$pre,
                           sm=recruits_sc$sm*flowers_sc$sm,
                           md=recruits_sc$md*flowers_sc$md,
                           lg=recruits_sc$lg*flowers_sc$lg)

fecundity_sc <- fecundity_sc %>%
  mutate(class=c('pre', 'sm', 'md', 'lg')) %>% 
  relocate(class)

#### Do the same for Florida Population ####

recruits_fl <- ribes_manual_classes_fl %>%
  filter(is.na(class)) %>%
  summarise(pre=length(class_next[class_next=='pre'])/length(class_next), 
            sm=length(class_next[class_next=='sm'])/length(class_next), 
            lg=length(class_next[class_next=='lg'])/length(class_next)) %>%
  mutate(var='p_size') %>% 
  pivot_longer(cols=!var, names_to = 'class') %>%
  select(!var) %>% 
  mutate(pre=value, sm=value, lg=value) %>%
  select(!value) %>% 
  column_to_rownames(var='class')
View(recruits_fl)

#Get the total number of flowers produced by individuals in each size class
# in year 1, as well as number of flowering individuals
flowers_fl_1 <- ribes_manual_classes_fl %>%
  filter(!is.na(class)) %>%
  group_by(class) %>% 
  summarise(flowers = sum(flw),
            n=sum(repr1))

#Same for year 2
flowers_fl_2 <- ribes_manual_classes_fl %>%
  filter(!is.na(class_next), class_next!='dead') %>%
  group_by(class_next) %>% 
  summarise(flowers = sum(flw_next),
            n=sum(repr2)) %>%
  rename(class=class_next)

#Add flower sums together, and then take an average
flowers_fl <- rbind(flowers_fl_2, flowers_fl_1) %>% 
  group_by(class) %>% 
  summarise(flowers=sum(flowers), n=sum(n),
            avg = flowers/n) %>%
  select(class,avg) %>%
  arrange(factor(class, c('pre', 'sm', 'lg')))

#Turn NA values into zeros
flowers_fl$avg[is.nan(flowers_fl$avg)]<-0

#Pivot flower data
flowers_fl <- flowers_fl %>% 
  pivot_wider(names_from = class,
              values_from = avg)

#Get average average germination rate per seed
germ_fl <- seeds %>% #data from in situ germination trials
  filter(Population=='fl') %>% 
  summarise(germ=sum(Ngerm),
            not=sum(notGerm),
            mean=sum(Ngerm)/(sum(Ngerm)+sum(notGerm)))
germ_fl = germ_fl$mean

#Get average number of seeds per fruit
avseeds_fl <- sdct %>% 
  filter(Population=='fl') %>%
  summarise(n=length(NumSeeds),
            sum=sum(NumSeeds),
            avg=sum(NumSeeds)/length(NumSeeds))
avseeds_fl = avseeds_fl$avg

#Get total number of flowers
flower_total_fl <- ribes_manual_classes_fl %>%
  filter(!is.na(class)) %>% 
  filter(Population=='fl') %>%
  summarise(flowers=sum(flw))
flower_total_fl = flower_total_fl$flowers

#Get total number of recruits
recruit_n_fl <- ribes_manual_classes_fl %>% 
  filter(is.na(class)) %>% 
  filter(Population=='fl') %>%
  summarise(n=sum(surv))
recruit_n_fl = recruit_n_fl$n

#Calculate overall recruitment rate = 
# number of recruits / number of flowers * seeds * germination rate
recruit_rate_fl = recruit_n_fl / (flower_total_fl * avseeds_fl * germ_fl)

#Multiply average number of flowers produced in each size class by
#average seeds, germination rate, and recruitment rate
#to get the number of recruits produced by each size class
flowers_fl <- flowers_fl * avseeds_fl * germ_fl * recruit_rate_fl

#Multiply average number of recruits produced by size proportions
#to get a matrix of contributions to each size class in year 2 by
#flowering individuals in year 1
fecundity_fl <- data.frame(pre=recruits_fl$pre*flowers_fl$pre,
                           sm=recruits_fl$sm*flowers_fl$sm,
                           lg=recruits_fl$lg*flowers_fl$lg)

fecundity_fl <- fecundity_fl %>%
  mutate(class=c('pre', 'sm', 'lg')) %>% 
  relocate(class)

#### Create Projection Matrix ####

###South Carolina Population

#Drop the row for dead individuals
transition_sc <- transition_sc %>% 
  filter(class_next!='dead') %>% 
  rename(class=class_next)

#Add together values from transition matrix and fecundity matrix
projection_sc <- data.frame(pre=transition_sc$pre+fecundity_sc$pre,
                            sm=transition_sc$sm+fecundity_sc$sm,
                            md=transition_sc$md+fecundity_sc$md,
                            lg=transition_sc$lg+fecundity_sc$lg)
rownames(projection_sc) <- c('pre', 'sm', 'md', 'lg')
View(projection_sc)

### Florida Population
transition_fl <- transition_fl %>% 
  filter(class_next!='dead') %>% 
  rename(class=class_next)
projection_fl <- data.frame(pre=transition_fl$pre+fecundity_fl$pre,
                            sm=transition_fl$sm+fecundity_fl$sm,
                            lg=transition_fl$lg+fecundity_fl$lg)
rownames(projection_fl) <- c('pre', 'sm', 'lg')
View(projection_fl)

#### Population Projection ####

##Create starting vectors by taking the total number of individuals in
##each size class. We did this earlier but need to pivot the dataframes
totals_sc <- totals_sc %>% 
  mutate(var='class') %>% 
  pivot_longer(cols=!var, names_to='class') %>% 
  select(class, value) %>% 
  rename(n=value)

totals_fl <- totals_fl %>% 
  mutate(var='class') %>% 
  pivot_longer(cols=!var, names_to='class') %>% 
  select(class, value) %>% 
  rename(n=value)

project_result_sc <- pop.projection(data.matrix(projection_sc), totals_sc$n)
eigen.analysis(data.matrix(projection_sc))

project_result_fl <- pop.projection(data.matrix(projection_fl), totals_fl$n)
lambda(data.matrix(projection_fl))

#### More popbio ####

## Create stage-fate dataframe
sc_trans <- ribes_manual_classes %>%
  filter(Population=='sc') %>%
  filter(!is.na(class)) %>% 
  mutate(pre=flw*avseeds_sc*germ_sc*recruit_rate_sc*recruits_sc$pre[1],
         sm=flw*avseeds_sc*germ_sc*recruit_rate_sc*recruits_sc$pre[2],
         md=flw*avseeds_sc*germ_sc*recruit_rate_sc*recruits_sc$pre[3],
         lg=flw*avseeds_sc*germ_sc*recruit_rate_sc*recruits_sc$pre[4])

stages<-c('pre', 'sm', 'md', 'lg')
proj_mat_sc <- projection.matrix(data.frame(sc_trans), 
                             stage=class, 
                             fate=class_next, 
                             sort=stages)
proj_sc <- pop.projection(proj_mat_sc, n=totals_sc$n)
boot_sc <- boot.transitions(data.frame(sc_trans), 
                 iterations=200, 
                 stage=class, 
                 fate=class_next,
                 sort=stages)

##Find 95% confidence intervals around the mean bootstrap estimate of lambda
mean_k_sc = mean(boot_sc$lambda)
sd_k_sc = sd(boot_sc$lambda)
se_k_sc = sd_k_sc / sqrt(100)
alpha = 0.05
df_k_sc = 99
t_k_sc = qt(p=alpha/2, df=df_k_sc,lower.tail = F)
error_margin_k_sc = t_k_sc * se_k_sc
lower_k_sc = mean_k_sc - error_margin_k_sc
upper_k_sc = mean_k_sc + error_margin_k_sc

print(c(lower_k_sc, mean_k_sc, upper_k_sc))
print(proj_sc$lambda)
quantile(boot_sc$lambda, c(0.025,0.975))

              

####End####

#### Vandermeer-Moloney ####
# 3. Implement Vandermeer-Moloney algorithm
#    ***IMPORTANT*** Give Moloney (1986) a read!!! I am not a math whiz!!!
#                            THIS COULD ALL BE VERY WRONG :O
ribes_size <- ribes %>%
  select(Population, TagNum, h, h_next) %>%
  filter(h>0)

## Aim: estimate probability P that individual starting in size class at 
##      first census will remain in same size class next census

### 3.1: define change in size
ribes_size <- ribes_size %>% mutate(dh = h_next - h)

### 3.2: dummy variable 't' records whether individual started in size class
ribes_size <- ribes_size %>% mutate(t = 0)

### 3.3: dummy variable 'r' records whether individual remained
ribes_size <- ribes_size %>% mutate(r=0)

### 3.4: dummy variable r* or 'r_s' is used to calculate distribution error
###      which is a a measure of the change in P that occurs as the size 
###      of an individual moves away from the midpoint size within a size class
ribes_size <- ribes_size %>% mutate(r_s=0)

### 3.5: dummy variable r** or 'r_ss' is used to calculate sample error
###      which is a measure of the inaccuracy in estimates of P that increases
###      as the width of a size interval approaches zero
ribes_size <- ribes_size %>% mutate(r_ss=0)

### 3.6: calculate values for t, r, and r* for size class with
###      user defined upper and lower bounds min and max

calc_dummies <- function(df, min, max) {
  #m star gives the size class midpoint
  m_s = (min+max)/2
  #define values for dummy variables
  df <- df %>% mutate(t = case_when(min <= h & max >= h ~ 1,
                                    .default = 0),
                      r = case_when(min <= h & max >= h &
                                      min <= h_next & max >= h_next ~ 1,
                                    .default = 0),
                      r_s = case_when(min <= h & max >= h &
                                        min <= m_s & max >= m_s ~ 1,
                                      .default = 0))
  return(df)
}

### 3.7 calculate values for r** for a number of subpopulations
calc_r_ss <- function(df, min, max, ns) {
  num_resample <- round(length(df$TagNum)/ns)
  df <- calc_dummies(df, min, max) %>% filter(t == 1)
  if(length(df$t) == 0){
    return(data.frame(r_ss=c(0), t=c(1), r=c(0), r_s=(0)))
  }
  df <- sample_n(df, size=num_resample, replace=T)
  df <- df %>% mutate(r_ss = case_when(min <= h & max >= h &
                                         min <= h_next & max >= h_next ~ 1,
                                       .default = 0))
  return(df)
}

### 3.8 calculate probabilities P, Q, and P*
calc_pqs <- function(df, min, max){
  df <- calc_dummies(df, min, max)
  pqs <- df %>% summarise(p = sum(r)/sum(t),
                          q = 1-(sum(r)/sum(t)),
                          ps = sum(r_s)/sum(t))
  return(pqs)
}

### 3.9 calculate distribution error
calc_de <- function(df, min, max, nj, nk){
  pqs <- calc_pqs(df, min, max)
  c <- (pqs$q^2+pqs$p^2)/((2*(nj-1)*nk)*(pqs$p*pqs$q))
  de <- c * (pqs$ps - pqs$p)^2
  return(de)
}

### 3.10 calculate probability P**
calc_pss <- function(df, min, max, ns){
  pss <- data.frame(pss=c(0))
  for(x in 1:ns){
    df <- calc_r_ss(df, min, max, ns)
    rbind(pss, df %>% summarise(pss = sum(r_ss)/sum(t)))
  }
  return(pss)
}

### 3.11 calculate sample error
calc_se <- function(df, min, max, nj, nk){
  pss <- calc_pss(df, min, max, nk)
  pqs <- calc_pqs(df, min, max)
  css <- (pqs$q^2+pqs$p^2)/((2*(nj-1)*nk*nk)*(pqs$p*pqs$q))
  se <- css * (sum(pss$pss)-pqs$p)^2
  return(se)
}

### 3.12 calculate distribution and probability error
###      for a range of max size class bounds
test_bounds <- function(df, min, max, nj, nk) {
  error <- data.frame(de=c(0), se=c(0), max=0)
  for(x in min:max){
    e <- data.frame(de=calc_de(df, min, min+x, nj, nk),
                    se=calc_se(df,min, min+x,nj,nk),
                    max=x)
    error <- rbind(error, e)
  }
  return(error)
}

### 3.13 Iteratively check for the max bound where DE and SE are minimized.
###      Probably best accomplished graphically.
###      Start with the smallest size class, where the minimum class bound is 
###      absoulte minimum size for the species)
a<-test_bounds(ribes_size, 0, 100, 2, 2)
a <- a %>% na.omit() %>% filter(max!=0) 
ggplot(a, aes(y=se, x=max)) + geom_line() + 
  geom_line(data=a, aes(y=de, x=max))

### 3.14 Set new mini