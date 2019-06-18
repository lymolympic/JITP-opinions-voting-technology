# Dependencies: recdata.Rdata

setwd("C:/Users/ines/Dropbox/Public Opinion US/VotingTechExp/Analysis/JITP Replication Package/Replication Package")
#setwd("C:/Dropbox/VotingTechExp/Analysis/JITP Replication Package/Replication Package")

library(survey)
library(weights)
library(boot)
library(MatchIt)
library(cem)
library(broom)
library(xtable)

#*******************************************************************
# Load Recoded Data
#*******************************************************************

load("recdata.Rdata")

#*******************************************************************
# Create Survey Design
#*******************************************************************

recdata.svydesign <- svydesign(ids = ~1, data = recdata, weights = recdata$weight)
recdata.svydesign.trim = trimWeights(recdata.svydesign, lower = 0.3, upper = 3)
recdata$weight.trim <- weights(recdata.svydesign.trim)
recdata.svydesign.trim$variables$weight.trim <- recdata$weight.trim

#*******************************************************************
# Summary Stats
#*******************************************************************

# Make an appendix table with summary statistics for the entire sample (N=982)

TableA1.1 <- svymean(~educcat+homeown+female+nonwhite+disabled+region+VBM+NoExcAbs+PermAbs+EarlyVoting+novote12+obama12+obama12Cty, recdata.svydesign.trim) * 100
TableA1.2 <- svymean(~age+lengthres+pid7, recdata.svydesign.trim)
#! Warning message is normal
TableA1.3 <- sqrt(svyvar(~age+lengthres+pid7, recdata.svydesign.trim))

TableA1 <- as.table(c(TableA1.1,TableA1.2))

rownames(TableA1) <- c("HS or less", "Some college", "College or more", "Homeowner", "Female", "Non-White", "Disabled", "Northeast", "Midwest", "South", "West", "VBM", "No-excuse absentee", "Permanent absentee", "Early voting", "Didn't vote 2012", "Voted Obama 2012", "Obama 2012 county", "Age", "Length of residency", "Party ID")

# Table A1
TableA1.sorted <- TableA1[c("Age","HS or less", "Some college", "College or more", "Homeowner", "Length of residency", "Female", "Non-White","Disabled", "Party ID", "Didn't vote 2012", "Voted Obama 2012", "Obama 2012 county","Northeast","Midwest","South","West","VBM","No-excuse absentee","Permanent absentee","Early voting")]

xtable(TableA1.sorted, digits = 1)

#*******************************************************************
# Balance fraud/convenience samples, including inexperienced
#*******************************************************************

# Make an appendix table with summary statistics for fraud/convenience samples
# THIS VERSION OF THE TABLE INCLUDES INEXPERIENCED RESPONDENTS

TableA4.1 <- t(svyby(~educcat+homeown+female+nonwhite+disabled+region+VBM+NoExcAbs+PermAbs+EarlyVoting+novote12+obama12+obama12Cty, ~FraudCue, recdata.svydesign.trim, svymean, keep.var=FALSE)) * 100
TableA4.2 <- t(svyby(~age+lengthres+pid7, ~FraudCue, recdata.svydesign.trim, svymean, keep.var=FALSE))

weighted_t_test <- function(x){
  test_output <- wtd.t.test(x[recdata$FraudCue == 1], x[recdata$FraudCue == 0], weight = recdata$weight.trim[recdata$FraudCue == 1], weighty = recdata$weight.trim[recdata$FraudCue == 0], alternative="two.tailed", samedata = FALSE)
  test_output$coefficients
}

TableA4.3 <- t(apply(recdata[,c("hs","somecol","college", "homeown", "female","nonwhite","disabled", "northeast","south","midwest","west", "VBM","NoExcAbs","PermAbs", "EarlyVoting","novote12","obama12","obama12Cty", "age", "lengthres", "pid7")], 2, weighted_t_test))

TableA4 <- cbind(rbind(TableA4.1[-1,],TableA4.2[-1,]), TableA4.3[,-2])

rownames(TableA4) <- c("HS or less", "Some college", "College or more", "Homeowner", "Female", "Non-White", "Disabled", "Northeast", "Midwest", "South", "West", "VBM", "No-excuse absentee", "Permanent absentee", "Early voting", "Didn't vote 2012", "Voted Obama 2012", "Obama 2012 county", "Age", "Length of residency", "Party ID")

TableA4.sorted <- TableA4[c("Age","HS or less", "Some college", "College or more", "Homeowner", "Length of residency", "Female", "Non-White","Disabled", "Party ID", "Didn't vote 2012", "Voted Obama 2012", "Obama 2012 county","Northeast","Midwest","South","West","VBM","No-excuse absentee","Permanent absentee","Early voting"),]

xtable(TableA4.sorted, digits = 2)

#**************************************************************************
# Effects of fraud/convenience cues, raw sample, including inexperienced
#**************************************************************************

# This is our old analysis of the effect of fraud/convenience cues
# It includes respondents without previous experiences
# We now report this in a footnote and in the appendix

Effects0.cues <- wtd.t.test(recdata$TVMChoice[recdata$FraudCue == 1], recdata$TVMChoice[recdata$FraudCue == 0], weight = recdata$weight.trim[recdata$FraudCue == 1], weighty = recdata$weight.trim[recdata$FraudCue == 0], alternative="two.tailed", samedata = FALSE)$additional

set.seed(1000)

func.cues.boot <- function(data, indices) {
  d <- data[indices,] # allows boot to select sample
  prop.treat <- mean(d$TVMChoice[d$FraudCue == 1])
  prop.control <- mean(d$TVMChoice[d$FraudCue == 0])
  diff <- prop.treat - prop.control
  return(c(prop.treat, prop.control, diff))
}

props0.boot <- boot(data = recdata, statistic= func.cues.boot, R=1000, parallel = "multicore", ncpus  = 6, sim = "ordinary", simple = TRUE, weights = recdata$weight.trim)

sd0.treats <- apply(props0.boot$t, 2, sd)

# Old Figure 2
ggplot(data=as.data.frame(props0.boot$t)) +
  geom_density(mapping=aes(x=V1), linetype="solid") +
  geom_density(mapping=aes(x=V2), linetype="dashed") +
  geom_hline(yintercept=0, colour="grey", size=1) +
  geom_vline(xintercept = Effects0.cues[2], linetype="solid") +
  geom_vline(xintercept = Effects0.cues[3], linetype="dashed") +
  annotate("text", x = 0.385, y = 18, label = paste("Fraud frame, mean = ",round(Effects0.cues[2]*100, 2),"%", sep = ""), size=2) +
  annotate("text", x = 0.505, y = 17, label = paste("Convenience frame, mean = ",round(Effects0.cues[3]*100, 2),"%", sep = ""), size=2) +
  annotate("text", x = 0.525, y = 19, label = paste("Difference: ", round(Effects0.cues[1]*100, 1), sep = ""), size=2) +
  annotate("text", x = 0.525, y = 18.5, label = paste("(standard error = ", round(sd0.treats[3]*100, 1), ")", sep = ""), size=2) +
  coord_cartesian(ylim = c(0, 20)) +
  labs(x = "Proportion of respondents who prefer TVMs", y = "Density") +
  theme_bw() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme(text = element_text(size=8))

ggsave("Fig2old.png", height=4, width=6, units='in', dpi=600)

#*******************************************************************
# Balance previous voting experiences
#*******************************************************************

# Make an appendix table with summary statistics for previous experiences

TableA2.1 <- svyby(~educcat+homeown+female+nonwhite+disabled+region+VBM+NoExcAbs+PermAbs+EarlyVoting+novote12+obama12+obama12Cty, ~experience, recdata.svydesign.trim, svymean, keep.var=FALSE)
TableA2.1 <- t(TableA2.1[,-1]) * 100

TableA2.2 <- svyby(~age+lengthres+pid7, ~experience, recdata.svydesign.trim, svymean, keep.var=FALSE)
TableA2.2 <- t(TableA2.2[,-1])

TableA2 <- rbind(TableA2.1,TableA2.2)

rownames(TableA2) <- c("HS or less", "Some college", "College or more", "Homeowner", "Female", "Non-White", "Disabled", "Northeast", "Midwest", "South", "West", "VBM", "No-excuse absentee", "Permanent absentee", "Early voting", "Didn't vote 2012", "Voted Obama 2012", "Obama 2012 county", "Age", "Length of residency", "Party ID")

# Table A2
TableA2.sorted <- TableA2[c("Age","HS or less", "Some college", "College or more", "Homeowner", "Length of residency", "Female", "Non-White","Disabled", "Party ID", "Didn't vote 2012", "Voted Obama 2012", "Obama 2012 county","Northeast","Midwest","South","West","VBM","No-excuse absentee","Permanent absentee","Early voting"),]

xtable(TableA2.sorted, digits = 2)

#*******************************************************************
# Relationship between voting experiences and tech preferences
#*******************************************************************

# Table 1
Table1 <- prop.table(svytable(~experience+f.Choice, design = recdata.svydesign.trim), 1)[c(3, 4, 1, 2),c(2, 3, 1)] * 100

xtable(Table1, digits = 2)

Table1N <- table(recdata$experience) # Number of respondents for each row
Table1N

#*******************************************************************
# Drop people who report no experiences with PB or TVM technologies
#*******************************************************************

datavoters <- subset(recdata, experience != "No Experience", select = c("weight", "region","northeast","south","midwest","west", "age", "educcat","hs","somecol","college", "female", "nonwhite", "disabled", "lengthres", "homeown", "novote12","obama12", "obama12Cty", "pid7", "VBM", "NoExcAbs", "PermAbs", "EarlyVoting","TVMexperience", "FraudCue", "TVMChoice", "PBChoice"))
datavoters.svydesign <- svydesign(ids = ~1, data = datavoters, weights = datavoters$weight)
datavoters.svydesign.trim <- trimWeights(datavoters.svydesign, lower = 0.3, upper = 3)
datavoters$weight.trim <- weights(datavoters.svydesign.trim)
datavoters.svydesign.trim$variables$weight.trim <- datavoters$weight.trim

#*******************************************************************
# Balance fraud/convenience samples, excluding inexperienced
#*******************************************************************

# Make an appendix table with summary statistics for fraud/convenience samples
# THIS VERSION OF THE TABLE EXCLUDES INEXPERIENCED RESPONDENTS

TableA3.1 <- t(svyby(~educcat+homeown+female+nonwhite+disabled+region+VBM+NoExcAbs+PermAbs+EarlyVoting+novote12+obama12+obama12Cty, ~FraudCue, datavoters.svydesign.trim, svymean, keep.var=FALSE)) * 100
TableA3.2 <- t(svyby(~age+lengthres+pid7, ~FraudCue, datavoters.svydesign.trim, svymean, keep.var=FALSE))

weighted_t_test2 <- function(x){
  test_output <- wtd.t.test(x[datavoters$FraudCue == 1], x[datavoters$FraudCue == 0], weight = datavoters$weight.trim[datavoters$FraudCue == 1], weighty = datavoters$weight.trim[datavoters$FraudCue == 0], alternative="two.tailed", samedata = FALSE)
  test_output$coefficients
}

TableA3.3 <- t(apply(datavoters[,c("hs","somecol","college", "homeown", "female","nonwhite","disabled", "northeast","south","midwest","west", "VBM","NoExcAbs","PermAbs", "EarlyVoting","novote12","obama12","obama12Cty", "age", "lengthres", "pid7")], 2, weighted_t_test2))

TableA3 <- cbind(rbind(TableA3.1[-1,],TableA3.2[-1,]), TableA3.3[,-2])

rownames(TableA3) <- c("HS or less", "Some college", "College or more", "Homeowner", "Female", "Non-White", "Disabled", "Northeast", "Midwest", "South", "West", "VBM", "No-excuse absentee", "Permanent absentee", "Early voting", "Didn't vote 2012", "Voted Obama 2012", "Obama 2012 county", "Age", "Length of residency", "Party ID")

TableA3.sorted <- TableA3[c("Age","HS or less", "Some college", "College or more", "Homeowner", "Length of residency", "Female", "Non-White","Disabled", "Party ID", "Didn't vote 2012", "Voted Obama 2012", "Obama 2012 county","Northeast","Midwest","South","West","VBM","No-excuse absentee","Permanent absentee","Early voting"),]

xtable(TableA3.sorted, digits = 2)

#**************************************************************************
# Effects of fraud/convenience cues, raw sample, excluding inexperienced
#**************************************************************************

# This is our new analysis of the effect of fraud/convenience cues

Effects.cues <- wtd.t.test(datavoters$TVMChoice[datavoters$FraudCue == 1], datavoters$TVMChoice[datavoters$FraudCue == 0], weight = datavoters$weight.trim[datavoters$FraudCue == 1], weighty = datavoters$weight.trim[datavoters$FraudCue == 0], alternative="two.tailed", samedata = FALSE)$additional

set.seed(1000)

props.boot <- boot(data = datavoters, statistic= func.cues.boot, R=1000, parallel = "multicore", ncpus  = 6, sim = "ordinary", simple = TRUE, weights = datavoters$weight.trim)

sd.treats <- apply(props.boot$t, 2, sd)

# New Figure 2
ggplot(data=as.data.frame(props.boot$t)) +
  geom_density(mapping=aes(x=V1), linetype="solid") +
  geom_density(mapping=aes(x=V2), linetype="dashed") +
  geom_hline(yintercept=0, colour="grey", size=1) +
  geom_vline(xintercept = Effects.cues[2], linetype="solid") +
  geom_vline(xintercept = Effects.cues[3], linetype="dashed") +
  annotate("text", x = 0.39, y = 17, label = paste("Fraud frame, mean = ",round(Effects.cues[2]*100, 2),"%", sep = ""), size=2) +
  annotate("text", x = 0.545, y = 16, label = paste("Convenience frame, mean = ",round(Effects.cues[3]*100, 2),"%", sep = ""), size=2) +
  annotate("text", x = 0.525, y = 19, label = paste("Difference: ", round(Effects.cues[1]*100, 1), sep = ""), size=2) +
  annotate("text", x = 0.525, y = 18.5, label = paste("(standard error = ", round(sd.treats[3]*100, 1), ")", sep = ""), size=2) +
  coord_cartesian(ylim = c(0, 20)) +
  labs(x = "Proportion of respondents who prefer TVMs", y = "Density") +
  theme_bw() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme(text = element_text(size=8))

ggsave("Fig2.png", height=4, width=6, units='in', dpi=600)

#**************************************************************************
# Effects of previous experiences, raw sample
#**************************************************************************

Effects.exps.raw <- wtd.t.test(datavoters$TVMChoice[datavoters$TVMexperience == 1], datavoters$TVMChoice[datavoters$TVMexperience == 0], weight = datavoters$weight.trim[datavoters$TVMexperience == 1], weighty = datavoters$weight.trim[datavoters$TVMexperience == 0], alternative="two.tailed", samedata = FALSE)$additional

set.seed(1000)

# Bootstrapping

func.props.boot <- function(data, indices) {
  d <- data[indices,] # allows boot to select sample
  prop.treat <- mean(d$TVMChoice[d$TVMexperience == 1])
  prop.control <- mean(d$TVMChoice[d$TVMexperience == 0])
  diff <- prop.treat - prop.control
  return(c(prop.treat, prop.control, diff))
}

props2.boot <- boot(data = datavoters, statistic = func.props.boot, R = 1000, parallel = "multicore", ncpus  = 6, sim = "ordinary", simple = TRUE, weights = datavoters$weight.trim)

sd.exps.raw <- apply(props2.boot$t, 2, sd)

Table.exp.effects.raw <- c(Effects.exps.raw[1] * 100, sd.exps.raw[3] * 100, dim(datavoters)[1])

# Figure 1 (Effect of previous voting experiences)
ggplot(data=as.data.frame(props2.boot$t)) +
  geom_density(mapping=aes(x=V1), linetype="solid") +
  geom_density(mapping=aes(x=V2), linetype="dashed") +
  geom_hline(yintercept=0, colour="grey", size=1) +
  geom_vline(xintercept = Effects.exps.raw[2], linetype="solid") +
  geom_vline(xintercept = Effects.exps.raw[3], linetype="dashed") +
  annotate("text", x = 0.365, y = 16, label = paste("PB experience, mean = ",round(Effects.exps.raw[3]*100, 1), "%", sep = ""), size=2) +
  annotate("text", x = 0.555, y = 15, label = paste("TVM experience, mean = ",round(Effects.exps.raw[2]*100, 1), "%", sep = ""), size=2) +
  annotate("text", x = 0.575, y = 17, label = paste("Difference: ", round(Effects.exps.raw[1]*100, 1), sep = ""), size=2) +
  annotate("text", x = 0.575, y = 16.5, label = paste("(standard error = ", round(sd.exps.raw[3]*100, 1), ")", sep = ""), size=2) +
  coord_cartesian(ylim = c(0, 18)) +
  labs(x = "Proportion of respondents who prefer TVMs", y = "Density") +
  theme_bw() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme(text = element_text(size=8))

ggsave("Fig1.png", height=4, width=6, units='in', dpi=600)

#************************************************************************************
# Effects of previous experiences, raw sample, regression estimate
#***********************************************************************************

# OLS

OLS.raw <- svyglm(TVMChoice ~ 1 + TVMexperience + age + female + nonwhite + disabled + lengthres + homeown + obama12 + obama12Cty + pid7 + educcat + region + VBM + NoExcAbs + PermAbs + EarlyVoting + FraudCue, design = datavoters.svydesign.trim)
summary(OLS.raw)

Effects.OLS.raw <- summary(OLS.raw)$coefficients["TVMexperience", 1:2]

Table.exp.effects.reg <- c(Effects.OLS.raw[1] * 100, Effects.OLS.raw[2] * 100, dim(datavoters)[1])

TableA5.OLScoefs <- summary(OLS.raw)$coefficients[,1:2]

# Logit

Logit.raw <- svyglm(TVMChoice ~ 1 + TVMexperience + age + female + nonwhite + disabled + lengthres + homeown + obama12 + obama12Cty + pid7 + educcat + region + VBM + NoExcAbs + PermAbs + EarlyVoting + FraudCue, design = datavoters.svydesign.trim, family = binomial(link = "logit"))
summary(Logit.raw)

TableA5.Logitcoefs <- summary(Logit.raw)$coefficients[,1:2]

TableA5 <- cbind(TableA5.OLScoefs, TableA5.Logitcoefs)

rownames(TableA5) <- c("Intercept", "TVM Experience", "Age", "Female", "Non-White", "Disabled", "Length of residency (years)", "Homeowner", "Voted for Obama", "Obama won in county",  "Party identification (1-7 scale)", "Education: Some college", "Education: College or more", "Region: Midwest", "Region: South", "Region: West", "All vote by mail", "No-excuse absentee voting", "Permanent absentee voting", "Early voting", "Fraud cue")

# Table A5
TableA5.sorted <- TableA5[c("Intercept", "TVM Experience", "Fraud cue", "Age", "Education: Some college", "Education: College or more", "Homeowner", "Length of residency (years)", "Female", "Non-White", "Disabled", "Voted for Obama", "Obama won in county",  "Party identification (1-7 scale)", "Region: Midwest", "Region: South", "Region: West", "All vote by mail", "No-excuse absentee voting", "Permanent absentee voting", "Early voting"),]
xtable(TableA5.sorted, digits = 3)

#************************************************************************************
# Table 2b: Difference in TVM Choice by Experiences, matched sample, bootstrapping
#************************************************************************************

#-------------- Subclassification Matching ------------#

set.seed(1000)

n.sub <- 9

m.out<-matchit(TVMexperience ~ age + female + nonwhite + disabled + lengthres + homeown + obama12 + obama12Cty + pid7 + educcat + region + VBM + NoExcAbs + PermAbs + EarlyVoting + FraudCue + weight.trim, distance = "logit", data = datavoters, method="subclass", subclass = n.sub, sub.by = "all")

matched.data <- match.data(m.out)

# Adjusted weights

subclass.weights <- aggregate(matched.data$weight[matched.data$TVMexperience == 1], by =list(matched.data$subclass[matched.data$TVMexperience == 1]), sum)$x

matched.data$adjweight <- matched.data$weight * matched.data$weights

# Treatment effects

Effects.subclass <- wtd.t.test(matched.data$TVMChoice[matched.data$TVMexperience == 1], matched.data$TVMChoice[matched.data$TVMexperience == 0], weight = matched.data$adjweight[matched.data$TVMexperience == 1], weighty = matched.data$adjweight[matched.data$TVMexperience == 0], alternative="two.tailed", samedata = FALSE)$additional

# Bootstrapped standard errors

func.submatch.boot <- function(data, indices) {
  n.sub <- 9
  d <- data[indices,] # allows boot to select sample
  m.out<-matchit(TVMexperience ~ age + female + nonwhite + disabled + lengthres + homeown + obama12 + obama12Cty + pid7 + educcat + region + VBM + NoExcAbs + PermAbs + EarlyVoting + FraudCue + weight.trim, distance = "logit", data = d, method="subclass", subclass = n.sub, sub.by = "all")
  matched.data <- match.data(m.out)
  prop.treat <- weighted.mean(matched.data$TVMChoice[matched.data$TVMexperience == 1], matched.data$weights[matched.data$TVMexperience == 1])
  prop.control <- weighted.mean(matched.data$TVMChoice[matched.data$TVMexperience == 0], matched.data$weights[matched.data$TVMexperience == 0])
  diff <- prop.treat - prop.control
  return(diff)
}

subclass.boot <- boot(data = datavoters, statistic= func.submatch.boot, R=1000, parallel = "multicore", ncpus  = 6, sim = "ordinary", simple = TRUE, weights = datavoters$weight.trim)

sd.exps.sub <- apply(subclass.boot$t, 2, sd)

Table.exp.effects.sub <- c(Effects.subclass[1] * 100, sd.exps.sub * 100, dim(matched.data)[1])

#-------------- Nearest-Neighbor Matching ------------#

set.seed(1000)

m.out2 <- matchit(TVMexperience ~ age + female + nonwhite + disabled + lengthres + homeown + obama12 + obama12Cty + pid7 + educcat + region + VBM + NoExcAbs + PermAbs + EarlyVoting + FraudCue + weight.trim, distance = "logit", data = datavoters,method="nearest", m.order="largest", caliper = 0.05, discard="control", ratio = 1, replace = FALSE)

matched.data2 <- rbind(datavoters[rownames(m.out2$match.matrix),],
                       datavoters[m.out2$match.matrix[,1],])

matched.data2$corr.rownames <- c(rownames(m.out2$match.matrix), m.out2$match.matrix[,1])

matched.data2$adjweight2 <- rep(datavoters[rownames(m.out2$match.matrix),]$weight,2) * rep(m.out2$weights[m.out2$treat == 1],2)

matched.data2 <- subset(matched.data2, adjweight2 > 0)

Effects.nn <- wtd.t.test(matched.data2$TVMChoice[matched.data2$TVMexperience == 1], matched.data2$TVMChoice[matched.data2$TVMexperience == 0], weight = matched.data2$adjweight2[matched.data2$TVMexperience == 1], weighty = matched.data2$adjweight2[matched.data2$TVMexperience == 0], alternative="two.tailed", samedata = FALSE)$additional

# Bootstrapped standard errors

func.nnmatch.boot <- function(data, indices) {
  d <- data[indices,] # allows boot to select sample
  m.out2 <- matchit(TVMexperience ~ age + female + nonwhite + disabled + lengthres + homeown + obama12 + obama12Cty + pid7 + educcat + region + VBM + NoExcAbs + PermAbs + EarlyVoting + FraudCue + weight.trim, distance = "logit", data = d,method="nearest", m.order="largest", caliper = 0.05, discard="control", ratio = 1, replace = FALSE)
  matched.data2 <- match.data(m.out2)
  prop.treat <- weighted.mean(matched.data2$TVMChoice[matched.data2$TVMexperience == 1], matched.data2$weights[matched.data2$TVMexperience == 1])
  prop.control <- weighted.mean(matched.data2$TVMChoice[matched.data2$TVMexperience == 0], matched.data2$weights[matched.data2$TVMexperience == 0])
  diff <- prop.treat - prop.control
  return(diff)
}

nnmatch.boot <- boot(data = datavoters, statistic= func.nnmatch.boot, R=1000, parallel = "multicore", ncpus  = 6, sim = "ordinary", simple = TRUE, weights = datavoters$weight.trim)

sd.exps.nn <- apply(nnmatch.boot$t, 2, sd)

Table.exp.effects.nn <- c(Effects.nn[1] * 100, sd.exps.nn * 100, dim(matched.data2)[1])

#-------------- Coarsened Exact Matching ------------#

set.seed(1000)

age.cut <- svyquantile(~age, datavoters.svydesign.trim, c(0,0.5,1))
lr.cut <- svyquantile(~lengthres, datavoters.svydesign.trim, c(0,0.5,1))
educ.group <- list(c("HS or less", "Some college"), c("College or more"))
pid7.group <- list(c(1, 2, 3), 4, c(5, 6,7))

datavoters$permissive.abs <- ifelse(datavoters$VBM == 1| datavoters$NoExcAbs == 1|datavoters$PermAbs == 1, 1, 0)

m.out3 <- cem(treatment="TVMexperience", data = datavoters, drop=c("TVMChoice","PBChoice","weight","weight.trim","VBM","PermAbs","NoExcAbs","region","northeast","south","midwest","west","hs","somecol","college","novote12"), cutpoints = list(age=age.cut,lengthres=lr.cut), grouping = list(educcat=educ.group,pid7=pid7.group)); m.out3

weights.by.strata <- aggregate(datavoters$weight[datavoters$TVMexperience == 1], by = list(m.out3$mstrata[datavoters$TVMexperience == 1]), mean)

TVMChoice.by.strata <- aggregate(datavoters$TVMChoice, by = list(m.out3$mstrata, datavoters$TVMexperience), mean)

prop.control.cem <- weighted.mean(TVMChoice.by.strata$x[TVMChoice.by.strata$Group.2 == 0], weights = weights.by.strata$x)
prop.treat.cem <- weighted.mean(TVMChoice.by.strata$x[TVMChoice.by.strata$Group.2 == 1], weights = weights.by.strata$x)

Effects.cem <- prop.treat.cem - prop.control.cem

# Bootstrapped standard errors

func.cem.boot <- function(data, indices) {
  d <- data[indices,] # allows boot to select sample
  
  age.cut <- quantile(d$age, c(0,0.5,1))
  lr.cut <- quantile(d$lengthres, c(0,0.5,1))
  educ.group <- list(c("HS or less", "Some college"), c("College or more"))
  pid7.group <- list(c(1, 2, 3), 4, c(5, 6,7))
  
  m.out3 <- cem(treatment="TVMexperience", data = d, drop=c("TVMChoice","PBChoice","weight","weight.trim","VBM","PermAbs","NoExcAbs","region","northeast","south","midwest","west","hs","somecol","college","novote12"), cutpoints = list(age=age.cut,lengthres=lr.cut), grouping = list(educcat=educ.group,pid7=pid7.group))
  matched.data3 <- subset(d, m.out3$mstrata != "NA")
  prop.treat <- mean(matched.data3$TVMChoice[matched.data3$TVMexperience == 1])
  prop.control <- mean(matched.data3$TVMChoice[matched.data3$TVMexperience == 0])
  diff <- prop.treat - prop.control
  return(diff)
}

cem.boot <- boot(data = datavoters, statistic= func.cem.boot, R=1000, parallel = "multicore", ncpus  = 6, sim = "ordinary", simple = TRUE, weights = datavoters$weight.trim)

sd.exps.cem <- apply(cem.boot$t, 2, sd)

Table.exp.effects.cem <- c(Effects.cem * 100, sd.exps.cem * 100, sum(m.out3$mstrata != "NA", na.rm = T))

#******************************************************************************************
# Effects of previous experiences
#******************************************************************************************

Table2 <- rbind(Table.exp.effects.raw, Table.exp.effects.reg, Table.exp.effects.sub, Table.exp.effects.nn, Table.exp.effects.cem)
rownames(Table2) <- c("Naive", "Regression", "Matching (subclass)", "Matching (nearest)", "Matching (CEM)")
colnames(Table2) <- c("Effect", "Std. Dev.", "N")

xtable(Table2, digits = 1)

#******************************************************************************************
# Create srvey design w/ adjusted weights from subclassification matching
#******************************************************************************************

matchdata.svydesign <- svydesign(ids = ~1, data = matched.data, weights = matched.data$adjweight)
matchdata.svydesign.trim = trimWeights(matchdata.svydesign, lower = 0.3, upper = 3)
matched.data$adjweight.trim <- weights(matchdata.svydesign.trim)
matchdata.svydesign.trim$variables$adjweight.trim <- matched.data$adjweight.trim

#******************************************************************************************
# Heterogeneous effects of previous experiencs, by experimental cues
#******************************************************************************************

# OLS

OLS.matched <- svyglm(TVMChoice ~ 1 + TVMexperience + age + female + nonwhite + disabled + lengthres + homeown + obama12 + obama12Cty + pid7 + educcat + region + VBM + NoExcAbs + PermAbs + EarlyVoting + FraudCue + FraudCue * TVMexperience, design = matchdata.svydesign.trim)
summary(OLS.matched)

TableA6.OLScoefs <- summary(OLS.matched)$coefficients[,1:2]

# Logit

Logit.matched <- svyglm(TVMChoice ~ 1 + TVMexperience + age + female + nonwhite + disabled + lengthres + homeown + obama12 + obama12Cty + pid7 + educcat + region + VBM + NoExcAbs + PermAbs + EarlyVoting + FraudCue + FraudCue * TVMexperience, design = matchdata.svydesign.trim, family = binomial(link = "logit"))
summary(Logit.matched)

TableA6.Logitcoefs <- summary(Logit.matched)$coefficients[,1:2]

TableA6 <- cbind(TableA6.OLScoefs, TableA6.Logitcoefs)

rownames(TableA6) <- c("Intercept", "TVM Experience", "Age", "Female", "Non-White", "Disabled", "Length of residency (years)", "Homeowner", "Voted for Obama", "Obama won in county",  "Party identification (1-7 scale)", "Education: Some college", "Education: College or more", "Region: Midwest", "Region: South", "Region: West", "All vote by mail", "No-excuse absentee voting", "Permanent absentee voting", "Early voting", "Fraud cue", "TVM experience x Fraud cue")

# Table A6
TableA6.sorted <- TableA6[c("Intercept", "TVM Experience", "Fraud cue", "TVM experience x Fraud cue", "Age", "Education: Some college", "Education: College or more", "Homeowner", "Length of residency (years)", "Female", "Non-White", "Disabled", "Voted for Obama", "Obama won in county",  "Party identification (1-7 scale)", "Region: Midwest", "Region: South", "Region: West", "All vote by mail", "No-excuse absentee voting", "Permanent absentee voting", "Early voting"),]
xtable(TableA6.sorted, digits = 3)

#******************************************************************************************
# Simulation of effects of previous experiences, by experimental cues
#******************************************************************************************

# Define a typical voter to be one with weighted modal value for categorical charactersitics
# and weighted mean value for other characteristics
# Simulate probabilities of selecting TVMs as technology for future elections

age.mean <- weighted.mean(matched.data$age, matched.data$adjweight.trim)
lengthres.mean <- weighted.mean(matched.data$lengthres, matched.data$adjweight.trim)

newdat <- data.frame(TVMexperience = c(1, 1, 0, 0), age = rep(age.mean, 4), female = rep(1, 4), nonwhite = rep(0, 4), disabled = rep(0, 4), lengthres = rep(lengthres.mean, 4), homeown = rep(1, 4), obama12 = rep(0, 4), obama12Cty = rep(1, 4), pid7 = rep(3, 4), educcat = factor(rep("HS or less", 4), levels(datavoters$educcat)), region = factor(rep("South", 4), levels(datavoters$region)), VBM = rep(0, 4), NoExcAbs = rep(1, 4), PermAbs = rep(0, 4), EarlyVoting = rep(1, 4), FraudCue = c(0, 1, 0, 1))

res.predicted <- predict(Logit.matched, newdata = newdat, type = "response")

T3 <- cbind(newdat$TVMexperience, newdat$FraudCue, res.predicted)
colnames(T3) <- c("TVM experience", "Fraud Cue", "P(TVM Choice")

Table3 <- matrix(c(T3[2,3], T3[4,3], T3[1,3], T3[3,3]), ncol = 2, byrow = T)
Table3 <- cbind(Table3, Table3[,1] - Table3[,2]) 
Table3 <- rbind(Table3, Table3[1,] - Table3[2,]) * 100

colnames(Table3) <- c("TVM exp.", "PB exp. only", "Difference")
rownames(Table3) <- c("Fraud cue", "Convenience cue", "Difference")

# Bootstrapped standard errors

set.seed(1000)

func.logit.boot <- function(data, indices) {
  d <- data[indices,] # allows boot to select sample 
  n.sub <- 9
  m.out<-matchit(TVMexperience ~ age + female + nonwhite + disabled + lengthres + homeown + obama12 + obama12Cty + pid7 + educcat + region + VBM + NoExcAbs + PermAbs + EarlyVoting + FraudCue + weight.trim, distance = "logit", data = d, method="subclass", subclass = n.sub, sub.by = "all")
  matched.data <- match.data(m.out)
  Logit.matched <- glm(TVMChoice ~ 1 + TVMexperience + age + female + nonwhite + disabled + lengthres + homeown + obama12 + as.numeric(obama12Cty) + pid7 + educcat + region + VBM + NoExcAbs + PermAbs + EarlyVoting + FraudCue + FraudCue * TVMexperience, data = matched.data, family = binomial(link = "logit"), weights = weights)
  res.predicted <- predict(Logit.matched, newdata = newdat, type = "response")
  prob.predicted <- as.vector(res.predicted)
  TVMPBdiff_under_fraud <- prob.predicted[2] - prob.predicted[4]
  TVMPBdiff_under_conv <- prob.predicted[1] - prob.predicted[3]
  FCdiff_for_TVMvoters <- prob.predicted[2] - prob.predicted[1]
  FCdiff_for_PBvoters <- prob.predicted[4] - prob.predicted[3]
  Diff_in_Diff <- TVMPBdiff_under_fraud - TVMPBdiff_under_conv
  return(c(prob.predicted, TVMPBdiff_under_fraud, TVMPBdiff_under_conv, FCdiff_for_TVMvoters, FCdiff_for_PBvoters, Diff_in_Diff))
}

# bootstrapping with 1000 replications 
logit.boot <- boot(data = datavoters, statistic= func.logit.boot, R=1000, parallel = "multicore", ncpus  = 6, sim = "ordinary", simple = TRUE, weights = datavoters$weight.trim)

sd.logit <- apply(logit.boot$t, 2, sd)

Table3.boot.se <- matrix(c(sd.logit[2], sd.logit[4], sd.logit[5], sd.logit[1], sd.logit[3], sd.logit[6], sd.logit[7], sd.logit[8], sd.logit[9]), ncol = 3, byrow = T) * 100

colnames(Table3.boot.se) <- c("TVM exp.", "PB exp. only", "Difference")
rownames(Table3.boot.se) <- c("Fraud cue s.e.", "Convenience cue s.e.", "Difference s.e.")

Table3.with.se <- rbind(Table3, Table3.boot.se)[c("Fraud cue", "Fraud cue s.e.", "Convenience cue", "Convenience cue s.e.", "Difference", "Difference s.e."),]

xtable(Table3.with.se, digits = 1)