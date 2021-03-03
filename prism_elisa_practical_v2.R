#############################
## LOAD DATA AND FUNCTIONS ##
#############################

# Set this path to the folder you have saved the files for this practical #
data_path <- '/Users/lsh289938/Documents/Teaching and Supervising/UCL/Session 10 Serology Practical/'

# This command loads the background functions #
source(paste0(data_path,'glbh0034_functions_v2.R'))

# This command loads the PRISM dataset #
data <- read.csv(paste0(data_path,"prism_elisa.csv"))

# This command allows you to view the data in a separate tab #
View(data)

############################
## DESCRIPTIVE STATISTICS ##
############################

# Age distribution #
hist(data$age,main="Population Age Distribution",xlab="Age, years",
     col=alpha("dodgerblue3",0.8),border="white",breaks=35)

# Malaria indicators and epidemiological variables #
table_sex <- addmargins( with(data, table(district,sex)), margin=1 )
print(table_sex)
round( prop.table(table_sex, margin=1)*100, digits=2)

table_age <- addmargins( with(data, table(district,age_group)), margin=1 )
print(table_age)
round( prop.table(table_age, margin=1)*100, digits=2)

table_bednet <- addmargins( with(data, table(district,bednet)), margin=1 )
print(table_bednet)
round( prop.table(table_bednet, margin=1)*100, digits=2)

table_microscopy <- addmargins( with(data, table(district,microscopy)), margin=1 )
print(table_microscopy)
round( prop.table(table_microscopy, margin=1)*100, digits=2)

table_rdt <- addmargins( with(data, table(district,rdt)), margin=1 )
print(table_rdt)
round( prop.table(table_rdt, margin=1)*100, digits=2)

table_hb <- addmargins( with(data, table(district,hb_category)), margin=1 )
print(table_hb)
round( prop.table(table_hb, margin=1)*100, digits=2)

#########################################
## ODDS OF SERO-POSITIVITY BY VARIABLE ##
#########################################

sero_sex <- addmargins( with(data, table(seropositive,sex)), margin=1 )
print(sero_sex)
sero_sex_prop <- round( prop.table(sero_sex, margin=1)*100, digits=2)
print(sero_sex_prop)

sero_age <- addmargins( with(data, table(seropositive,age_group)), margin=1 )
print(sero_age)
sero_age_prop <- round( prop.table(sero_age, margin=1)*100, digits=2)
print(sero_age_prop)

sero_bednet <- addmargins( with(data, table(seropositive,bednet)), margin=1 )
print(sero_bednet)
sero_bednet_prop <- round( prop.table(sero_bednet, margin=1)*100, digits=2)
print(sero_bednet_prop)

sero_microscopy <- addmargins( with(data, table(seropositive,microscopy)), margin=1 )
print(sero_microscopy)
sero_mic_prop <- round( prop.table(sero_microscopy, margin=1)*100, digits=2)
print(sero_mic_prop)

sero_rdt <- addmargins( with(data, table(seropositive,rdt)), margin=1 )
print(sero_rdt)
sero_rdt_prop <- round( prop.table(sero_rdt, margin=1)*100, digits=2)
print(sero_rdt_prop)

sero_hb <- addmargins( with(data, table(seropositive,hb_category)), margin=1 )
print(sero_hb)
sero_hb_prop <- round( prop.table(sero_hb, margin=1)*100, digits=2)
print(sero_hb_prop)

# Using logistic regression, we can estimate the odds of an individual being sero-positive #
# dependent on the epidemiological variables described above #

fit_all <- glm(seropositive ~ sex + age_group + hemoglobin + bednet + district, data=data,family="binomial")
summary(fit_all)

fit_jinja <- glm(seropositive ~ sex + age_group + hemoglobin + bednet, data=subset(data,district=="1: Jinja"),family="binomial")
summary(fit_jinja)

fit_kanungu <- glm(seropositive ~ sex + age_group + hemoglobin + bednet, data=subset(data,district=="2: Kanungu"),family="binomial")
summary(fit_kanungu)

fit_tororo <- glm(seropositive ~ sex + age_group + hemoglobin + bednet, data=subset(data,district=="3: Tororo"),family="binomial")
summary(fit_tororo)

##############################
## DEFINING SERO-POSITIVITY ##
##############################

hist(log(data$optical_density),xlab="optical density (OD)", main="Distribution Ab Response\nAntigen: PfAMA1",
     breaks=35,col=alpha("mediumpurple4",0.8),border="white",xlim=c(-8,2),freq=F)

## The normalmixEM function runs the finite mixture model for your data ##
## This estimates a bi-model distribution (i.e. assumes two normal distributions in your data) ##
fmm_fit <- normalmixEM(na.omit(log(data$optical_density)))
fmm_fit_test <- normalmixEM2comp(na.omit(log(data$optical_density)))

# mu1, sig1 - mean and standard deviation (SD) of normal distribution 1 #
# mu2, sig2 - mean and standard deviation (SD) of normal distribution 2 #
mu1 <- fmm_fit$mu[1]
mu2 <- fmm_fit$mu[2]
sig1 <- fmm_fit$sigma[1]
sig2 <- fmm_fit$sigma[2]

# Mean and SD of the lower distribution is your sero-negative population #
# used to defined your sero-positivity threshold #
min_comp1 <- which(fmm_fit$mu == min(fmm_fit$mu))

# Individuals with OD values greater than cufoff are sero-positive #
cutoff <- fmm_fit$mu[min_comp1] + sqrt(fmm_fit$sigma[min_comp1])

# Plot the distribution of antibodies responses and sero-positivity threshold #
plot_fmm(data)

#####################
## SERO-PREVALENCE ##
#####################

## JINJA ##

jinja.data <- subset(data,district=="1: Jinja")
jinja.data.sero <- create.data.object(round(jinja.data$age),jinja.data$seropositive)
jinja.age.profile <- create.age.profile(jinja.data.sero,lag=0.15,analysis='overall')
jinja.age.profile$age.profiles

plot_seroprev(plot.age.profile=jinja.age.profile, district="Jinja",colour="#0144D2")

## KANUNGU ##

kanungu.data <- subset(data,district=="2: Kanungu")
kanungu.data.sero <- create.data.object(round(kanungu.data$age),kanungu.data$seropositive)
kanungu.age.profile <- create.age.profile(kanungu.data.sero,lag=0.15,analysis='overall')
kanungu.age.profile$age.profiles

plot_seroprev(plot.age.profile=kanungu.age.profile, district="Kanungu",colour="#AB1779")

## TORORO ##

tororo.data <- subset(data,district=="3: Tororo")
tororo.data.sero <- create.data.object(round(tororo.data$age),tororo.data$seropositive)
tororo.age.profile <- create.age.profile(tororo.data.sero,lag=0.15,analysis='overall')
tororo.age.profile$age.profiles

plot_seroprev(plot.age.profile=tororo.age.profile, district="Tororo",colour="#E4A804")

###############################################################
## ESTIMATING FORCE OF INFECTION USING SERO-CONVERSION RATES ##
###############################################################

# Individuals ages less than  1 year excluded due to maternal antibodies #
sero.data <- subset(data,age>=1)

# Estimate sero-conversion rate (force of infection) by study site #
# lambda - sero-conversion rate with 95% confidence intervals #
# rho - sero-reversion rate with 95% confidence intervals #
data2 <- create.data.object(round(sero.data$age),sero.data$seropositive,sero.data$district)
scr_fit <- simple.rcm.analysis(data2,analysis='split-unshared-rho',int.rho=c(0.001,0.250))
print(scr_fit$estimates)

# PLOT SERO-CONVERSION RATE BY STUDY SITE #
par(mfrow=c(1,3))
plot_scr(plot.age.profile=jinja.age.profile,district="Jinja",colour="#0144D2")
plot_scr(plot.age.profile=kanungu.age.profile,district="Kanungu",colour="#AB1779")
plot_scr(plot.age.profile=tororo.age.profile,district="Tororo",colour="#E4A804")

# Fit sero-conversion rate with change in SCR #
 
rho1 <- scr_fit$estimates[1,6] # lower rho estimate for Jinja
rho2 <- scr_fit$estimates[1,7] # upper rho estimate for Jinja

# set time.common to the age where you can see a potential change point in your data
# in this example, I put age 9 years
# lambda1 represents the historical transmission intensity (estimated for individuals born before time of change)
# lambda 2 represents the current transmission intensity (estimated for individuals born since time of change)
scr_fit2 <- two.rcm.analysis(jinja.data.sero, analysis="split-shared-rho",int.rho=c(0.001,0.05),time.common=9)
scr_fit2$estimates
plot_scr2(plot.age.profile=jinja.age.profile,district="Jinja",colour="#0144D2")
