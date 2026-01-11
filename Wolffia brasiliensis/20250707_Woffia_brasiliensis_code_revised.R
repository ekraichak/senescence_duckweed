
##### Wolffia Project
##### Athita Senayai, Suzanne Chmilar, and Robert Laird
##### Note: Parts of script modified from Caloric restriction project (part 2)



##### PART I. LOAD PACKAGES AND SET WORKING DIRECTORY #####

##### 1.   Load packages #####

library(reshape)
library(geepack)
library(gee)


##### 2.   Set working directory #####

# Set working directory
setwd('/Users/athita/Downloads/Wolffia brasiliensis')



##### PART II. DATA ANALYSIS: SURVIVAL AND REPRODUCTION #####



##### 3.   Load daily survival and reproduction data, and perform preliminary processing #####

# Load daily reproduction data
wolffia.full <- read.csv('Wolffia_daily_surv_and_repro_ver_5.csv', header = T) 

# Exclude appropriate fronds
wolf <- wolffia.full[wolffia.full$exclude == "No",]

# Change light treatment names to be less cryptic 
wolf$light.treatment[wolf$light.treatment == "4"] <- "1/4"
wolf$light.treatment <- factor(wolf$light.treatment, c("1", "1/4"))

# Extract reproduction data 
wolf.repro <- subset(wolf, select = Mar.24.2023:Jun.01.2023) 

# Double check dates of birth
wolf$date.birth.alt <- apply(wolf.repro, 1, function (x) names(wolf.repro)[min(which(x == 0))])
wolf$date.birth.alt == wolf$date.birth # all entries match
wolf <- wolf[ , 1:(ncol(wolf) - 1)] # drop unnecessary date.birth.alt column

# Double check dates of last reproduction
wolf$date.last.repro.alt <- apply(wolf.repro, 1, function (x) names(wolf.repro)[max(which(x > 0))])
wolf$date.last.repro.alt == wolf$date.last.repro # all entries match 
wolf <- wolf[ , 1:(ncol(wolf) - 1)] # drop unnecessary date.last.repro.alt column 

# Calculate dates of first reproduction
wolf$date.first.repro <- apply(wolf.repro, 1, function (x) names(wolf.repro)[min(which(x > 0))])

# Convert date values to R's date format
wolf$date.birth <- as.Date(wolf$date.birth, format = "%b.%d.%Y")
wolf$date.last.repro <- as.Date(wolf$date.last.repro, format = "%b.%d.%Y")
wolf$date.first.repro <- as.Date(wolf$date.first.repro, format = "%b.%d.%Y")

# Calculate lifespan from birth
wolf$lifespan <- as.numeric(wolf$date.last.repro - wolf$date.birth + 1) # first day is Day 1, not Day 0

# Calculate lifetime reproduction
wolf$total.offspring <- apply(wolf.repro, 1, function (x) sum(x, na.rm = T))

# Make new data frame with aligned reproduction data relative to dates of birth
wolf.repro.aligned <- matrix(NA, nrow = nrow(wolf.repro), ncol = max(wolf$lifespan))
birth.col <- as.numeric(wolf$date.birth - min(wolf$date.birth) + 1) # column index for date of birth       
death.col <- as.numeric(wolf$date.last.repro - min(wolf$date.birth) + 1) # column index for date of death
for (i in 1:nrow(wolf.repro)) {
  wolf.repro.aligned[i, 1:wolf$lifespan[i]] <- as.numeric(wolf.repro[i, birth.col[i]:death.col[i]])
}

# Calculate intrinsic rate of increase of individuals (r)
for (i in 1:nrow(wolf.repro.aligned)) {
  ind <- wolf.repro.aligned[i, ]                                # extract individual in row i
  L <- max(which(ind > 0))                                     # get lifespan
  lesMat <- matrix(0, nrow = L, ncol = L)                      # create empty lifespan x lifespan Leslie matrix
  lesMat[1, ] <- ind[1:L]                                      # fill first row with reproduction data
  if (L > 1) {
    lesMat[2:L, 1:(L - 1)] <- diag(rep(1, (L - 1)))            # fill sub-diagonal with survival data
  }
  eigenvalues <- eigen(lesMat)$values                          # get list of eigenvalues of Leslie matrix
  lambda <- max(Re(eigenvalues[abs(Im(eigenvalues)) < 1e-6]))  # extract lambda
  wolf$r[i] <- log(lambda)                                      # record intrinsic rate of increase (r)
  if (abs(wolf$r[i]) < 1e-6) {
    wolf$r[i] <- 0                                              # correct for floating point arithmetic rounding errors
  }
} 



##### 4.  Analyze reproductive lifespan (Mann-Whitney) #####

# Visually inspect boxplot of lifespan
tiff(file = "lifespan_1_4.tiff", width = 5, height = 4, units = "in", res = 600, compression = "lzw")

boxplot(lifespan ~ light.treatment, dat = wolf, ylim = c(0, 50),
        las = 1, whisklty = 1, col = "white", xlab = "Light intensity treatment", ylab = "Lifespan (days)")

dev.off()


### Summary statistical values in two treatments ###
by(wolf$lifespan, wolf$light.treatment, summary)

#wolf$light.treatment: 1
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#7.00   22.00   25.00   26.11   29.00   49.00 

#wolf$light.treatment: 1/4
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#11.00   25.00   27.00   28.45   30.00   49.00

### Data Dianostics divided by treatment ###
shapiro_test <- by(wolf$lifespan, wolf$light.treatment, shapiro.test)
shapiro_test # Not normally distributed (p< 0.05)

#wolf$light.treatment: 1
#Shapiro-Wilk normality test
#W = 0.91678, p-value = 7.111e-08

#wolf$light.treatment: 1/4
#W = 0.86188, p-value = 9.432e-11

treatment_groups <- split(wolf$lifespan, wolf$light.treatment)
par(mfrow=c(1,2)) # Set up a 1x2 grid for plotting histograms side by side
for (i in 1:length(treatment_groups)) {
  hist(treatment_groups[[i]], main=paste("Treatment", names(treatment_groups)[i], "Histogram"),
       xlab="Plant Lifespan", ylab="Frequency", col="lightblue", border="white")
}

variance_test <- var.test(lifespan ~ light.treatment, data = wolf)
variance_test #Equal variance  (p>0.05)
#data:  lifespan by light.treatment
#F = 1.262, num df = 157, denom df = 154, p-value = 0.1482


# Assumption of normality unable to be met for t-test. Non-parametric Mann-Whitney test used instead.
# Mann-Whitney test on lifespan
wilcox.test(lifespan ~ light.treatment, data = wolf)

# Wilcoxon rank sum test with continuity correction
# data:  lifespan by light.treatment
# W = 9030.5, p-value = 5.697e-05
# alternative hypothesis: true location shift is not equal to 0



##### 5.   Analyze lifetime reproduction (Mann-Whitney) #####

# Visually inspect boxplot of lifetime reproduction
tiff(file = "Offspring_1_4.tiff", width = 5, height = 4, units = "in", res = 600, compression = "lzw")

boxplot(total.offspring ~ light.treatment, dat = wolf, ylim = c(0, 25), las = 1, whisklty = 1, col = "white", 
        xlab = "Light intensity treatment", ylab = "Total number of offspring")

dev.off()


### Summary statistical values in two treatments ###
by(wolf$total.offspring, wolf$light.treatment, summary)

# wolf$light.treatment: 1
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#3.00   12.00   13.00   13.47   14.75   25.00 

#  wolf$light.treatment: 1/4
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#4.00   12.00   13.00   12.85   14.00   21.00 

### Data Dianostics divided by treatment ###
shapiro_test <- by(wolf$total.offspring, wolf$light.treatment, shapiro.test)
shapiro_test # Not normally distributed (p< 0.05)

#wolf$light.treatment: 1
#W = 0.85922, p-value = 5.322e-11

#wolf$light.treatment: 1/4
#W = 0.8768, p-value = 4.958e-10

treatment_groups <- split(wolf$total.offspring, wolf$light.treatment)
par(mfrow=c(1,2)) # Set up a 1x2 grid for plotting histograms side by side
for (i in 1:length(treatment_groups)) {
  hist(treatment_groups[[i]], main=paste("Treatment", names(treatment_groups)[i], "Histogram"),
       xlab="Plant Lifespan", ylab="Frequency", col="lightblue", border="white")
}

variance_test <- var.test(total.offspring ~ light.treatment, data = wolf)
variance_test #Not equal variance  (p<0.05)
#data:  total.offspring by light.treatment
#F = 2.2492, num df = 157, denom df = 154, p-value = 6.618e-07


# Assumption of normality unable to be met for t-test. Non-parametric Mann-Whitney test used instead.
# Mann-Whitney test on total.offspring
wilcox.test(total.offspring ~ light.treatment, data = wolf)

# Wilcoxon rank sum test with continuity correction
# data:  total.offspring by light.treatment
# W = 13652, p-value = 0.07358
# alternative hypothesis: true location shift is not equal to 0



##### 6.   Analyze intrinsic rate of increase of individuals (r) (t-test) #####

# Visually inspect boxplot of intrinsic rate of increase of individuals (r)

tiff(file = "Intrinsic_1_4.tiff", width = 5, height = 4, units = "in", res = 600, compression = "lzw")

boxplot(r ~ light.treatment, dat = wolf, ylim = c(0, 0.5), las = 1,  whisklty = 1, col = "white", 
        xlab = "Light intensity treatment", ylab = "Intrinsic rate of increase (r)")

dev.off()

### Summary statistical values in two treatments ###
by(wolf$r, wolf$light.treatment, summary)

# wolf$light.treatment: 1
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.1995  0.2743  0.2866  0.2928  0.3099  0.4518 

#wolf$light.treatment: 1/4
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.1709  0.2561  0.2695  0.2738  0.2909  0.4066

### Data Dianostics divided by treatment ###
shapiro_test <- by(wolf$r, wolf$light.treatment, shapiro.test)
shapiro_test # Not normally distributed (p< 0.05)

#wolf$light.treatment: 1
#W = 0.91356, p-value = 4.43e-08

#wolf$light.treatment: 1/4
#W = 0.93999, p-value = 3.779e-06

treatment_groups <- split(wolf$r, wolf$light.treatment)
par(mfrow=c(1,2)) # Set up a 1x2 grid for plotting histograms side by side
for (i in 1:length(treatment_groups)) {
  hist(treatment_groups[[i]], main=paste("Treatment", names(treatment_groups)[i], "Histogram"),
       xlab="Plant Lifespan", ylab="Frequency", col="lightblue", border="white")
}

variance_test <- var.test(r ~ light.treatment, data = wolf)
variance_test # Equal variance  (p>0.05)
#data:  r by light.treatment
#F = 1.0003, num df = 157, denom df = 154, p-value = 0.9988


# Assumption of normality unable to be met for t-test. Non-parametric Mann-Whitney test used instead.
# Mann-Whitney test on r
wilcox.test(r ~ light.treatment, data = wolf)

# 	Wilcoxon rank sum test with continuity correction

# data:  r by light.treatment
# W = 16910, p-value = 5.647e-09
# alternative hypothesis: true location shift is not equal to 0

##### Combined three figures: Lifespan, total number of offspring, intrinsic rate of increase of individuals (r)

library(ggplot2)
library(gridExtra)
library(grid)

wolf$light.treatment <- factor(wolf$light.treatment,
                               levels = c("1", "1/4"),
                               labels = c("Full light", "Dim light"))


create_plot <- function(wolf, y_var, y_lab) {
  ggplot(wolf, aes(x = light.treatment, y = {{y_var}})) +
    geom_boxplot(outlier.shape = NA) +
    stat_boxplot(geom = 'errorbar', width = 0.2) +  # Add error bars
    geom_jitter(width = 0.2, alpha = 0.5, color = "darkgrey") +
    stat_summary(fun = mean, geom = "point", shape = 16, size = 3, color = "blue") +
    labs(y = y_lab, x = NULL) +
    expand_limits(y = 0) +  # Force y-axis to include 0
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +  # Remove padding below 0
    theme_bw() +
    theme(
      axis.title.x = element_blank(),
      axis.text.x = element_text(color = "black"),
      plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")
    )
}

# Create individual plots
p1 <- create_plot(wolf, lifespan, "Lifespan (days)")
p2 <- create_plot(wolf, total.offspring, "Total number of offspring")
p3 <- create_plot(wolf, r, expression("Intrinsic rate of increase (" * italic("r") * ")"))

# Combine plots using grid.arrange and add labels using grid.text
combined_plot <- grid.arrange(
  arrangeGrob(p1, top = textGrob("(a)", x = unit(0, "npc"), just = "left", gp = gpar(fontsize = 12, face = "bold"))),
  arrangeGrob(p2, top = textGrob("(b)", x = unit(0, "npc"), just = "left", gp = gpar(fontsize = 12, face = "bold"))),
  arrangeGrob(p3, top = textGrob("(c)", x = unit(0, "npc"), just = "left", gp = gpar(fontsize = 12, face = "bold"))),
  ncol = 1, 
  bottom = textGrob("Light intensity treatment", gp = gpar(fontsize = 12)),
  top = textGrob("Wolffia Under Different Light Intensities", gp = gpar(fontsize = 12)),
  heights = c(1, 1, 1)
)

# Display the combined plot
print(combined_plot)
ggsave("Wolffia_combined_plots_senesence_20250311.png", combined_plot, width = 5, height = 10, units = "in", dpi = 600)



##### 7.   Create aligned reproduction data for each light intensity treatment #####

### IMPORTANT: If you encounter the error "no non-missing arguments to max; returning -Inf" 
### when running this section, it means there's a factor level mismatch because the 
### light.treatment levels were changed in the plotting section (around Line 270).
### 
### EASY FIX: Clear your workspace and re-run steps 1-3, then run this step (Step 7) 
### before running any plotting sections.
###
### Alternative: Update the factor levels in this section to match Line 270


# Create separate aligned reproduction data for each light intensity treatment
wolf.repro.aligned.01 <- wolf.repro.aligned[wolf$light.treatment == "1", ]
wolf.repro.aligned.04 <- wolf.repro.aligned[wolf$light.treatment == "1/4", ]

# Convert aligned reproduction data into binomial format
wolf.repro.aligned.binom.01 <- apply(wolf.repro.aligned.01, 2, function (x) replace(x, which(x > 1), 1))
wolf.repro.aligned.binom.04 <- apply(wolf.repro.aligned.04, 2, function (x) replace(x, which(x > 1), 1))



##### 8A.   Create life tables for light intensity treatment '1' #####

n.01 <- nrow(wolf.repro.aligned.01) # total starting sample size
life.dist.01 <- wolf$lifespan[wolf$light.treatment == "1"] # distribution of lifespan
maxAge.01 <- max(life.dist.01) # max lifespan
meanAge.01 <- mean(life.dist.01) # mean lifespan
lifeTab.01 <- data.frame(matrix(nrow = maxAge.01, ncol = 7))
names(lifeTab.01) <- c("age", "age.rel", "n.last.repro", "n.surv.cum", "p.surv.cum", "n.repro", "p.repro")    
lifeTab.01$age <- 1:maxAge.01 # ages
lifeTab.01$age.rel <- lifeTab.01$age/meanAge.01 # relative ages (in mean lifespans)
lifeTab.01$n.last.repro <- tabulate(life.dist.01, nbins = maxAge.01) # number of deaths at each age
lifeTab.01$n.surv.cum <- c(n.01, n.01 - cumsum(lifeTab.01$n.last.repro[1:(maxAge.01 - 1)])) # total cumulative survivorship at age
lifeTab.01$p.surv.cum <- lifeTab.01$n.surv.cum/n.01 # proportional cumulative survivorship at age
lifeTab.01$n.repro <- colSums(wolf.repro.aligned.binom.01[ , 1:maxAge.01], na.rm = TRUE) # number reproducing at each age
lifeTab.01$p.repro <- lifeTab.01$n.repro/lifeTab.01$n.surv.cum # proportion reproducing at age


##### 8B.   Create life tables for light intensity treatment '1/4' #####

n.04 <- nrow(wolf.repro.aligned.04) # total starting sample size
life.dist.04 <- wolf$lifespan[wolf$light.treatment == "1/4"] # distribution of lifespan
maxAge.04 <- max(life.dist.04) # max lifespan
meanAge.04 <- mean(life.dist.04) # mean lifespan
lifeTab.04 <- data.frame(matrix(nrow = maxAge.04, ncol = 7))
names(lifeTab.04) <- c("age", "age.rel", "n.last.repro", "n.surv.cum", "p.surv.cum", "n.repro", "p.repro")    
lifeTab.04$age <- 1:maxAge.04 # ages
lifeTab.04$age.rel <- lifeTab.04$age/meanAge.04 # relative ages (in mean lifespans)
lifeTab.04$n.last.repro <- tabulate(life.dist.04, nbins = maxAge.04) # number of deaths at each age
lifeTab.04$n.surv.cum <- c(n.04, n.04 - cumsum(lifeTab.04$n.last.repro[1:(maxAge.04 - 1)])) # total cumulative survivorship at age
lifeTab.04$p.surv.cum <- lifeTab.04$n.surv.cum/n.04 # proportional cumulative survivorship at age
lifeTab.04$n.repro <- colSums(wolf.repro.aligned.binom.04[ , 1:maxAge.04], na.rm = TRUE) # number reproducing at each age
lifeTab.04$p.repro <- lifeTab.04$n.repro/lifeTab.04$n.surv.cum # proportion reproducing at age



##### 9A.   Fit survival models for light intensity treatment '1' #####

# control for optim function to get max log-likelihood instead of min
C <- list(fnscale = -1) 

# times to death
t.01 <- life.dist.01 

# log-likelihood functions for exponential, Weibull, Gompertz, and logistic models 
loglik.exp.01  <- function(par) sum(log(par[1]*exp(-par[1]*t.01)))
loglik.weib.01 <- function(par) sum(log(par[1]^par[2]*par[2]*t.01^(par[2] - 1)*exp(-(par[1]*t.01)^par[2])))
loglik.gomp.01 <- function(par) sum(log(par[1]*exp(par[2]*t.01 - (par[1]/par[2])*(exp(par[2]*t.01) - 1) )))
loglik.log.01  <- function(par) sum(log(par[1]*exp(par[2]*t.01)*(1 + (par[1]*par[3]/par[2])*(exp(par[2]*t.01) - 1) )^(-(par[3] + 1)/par[3])))

# maximize log-likelihood functions for each of the four survival models
# use 'Brent' method for exponential;  default 'Nelder-Mead' method is not suitable for single param optimization
ml.fit.exp.01  <- optim(c(0.1), loglik.exp.01, control = C, method = "Brent", lower = 0, upper = 1)
ml.fit.weib.01 <- optim(c(0.1, 0.1), loglik.weib.01, control = C)
ml.fit.gomp.01 <- optim(c(0.001, 0.1), loglik.gomp.01, control = C)
ml.fit.log.01  <- optim(c(0.0001, 0.75, 2), loglik.log.01, control = C)

# extract parameter estimates (parameters a, b, and c, if applicable)
exp.a.01  <- ml.fit.exp.01$par[1]
weib.a.01 <- ml.fit.weib.01$par[1]; weib.b.01 <- ml.fit.weib.01$par[2]
gomp.a.01 <- ml.fit.gomp.01$par[1]; gomp.b.01 <- ml.fit.gomp.01$par[2]
log.a.01  <- ml.fit.log.01$par[1];   log.b.01 <- ml.fit.log.01$par[2]; log.c.01 <- ml.fit.log.01$par[3]

# specify best-fit curves
exp.surv.curve.01 <- function(t.01) exp(-exp.a.01*t.01)
weib.surv.curve.01 <- function(t.01) exp(-(weib.a.01*t.01)^weib.b.01)
gomp.surv.curve.01 <- function(t.01) exp(-(gomp.a.01/gomp.b.01)*(exp(gomp.b.01*t.01) - 1))
log.surv.curve.01 <- function(t.01) (1 + log.c.01*log.a.01/log.b.01*(exp(log.b.01*t.01) - 1))^(-1/log.c.01)

# visually assess best-fit curves

tiff(file = "Surviving_1.tiff", width = 5, height = 4, units = "in", res = 600, compression = "lzw")

plot(p.surv.cum ~ age, data = lifeTab.01, log = "y", las = 1, pch = 16, xlab = "Age (days)", ylab = "Proportion surviving")
curve(exp.surv.curve.01, add = T, lwd = 2, col = 1)
curve(weib.surv.curve.01, add = T, lwd = 2, col = 2)
curve(gomp.surv.curve.01, add = T, lwd = 2, col = 3)
curve(log.surv.curve.01, add = T, lwd = 2, col = 4)
legend("bottomleft", c("Exponential", "Weibull", "Gompertz", "Logistic"), col = 1:4, bty = "n", lwd = 2)

dev.off()

# compare models with AICc
aicTab.01 <- data.frame(matrix(ncol = 8, nrow = 4))
names(aicTab.01) <- c("k", "log.lik", "deviance", "aic", "aicc", "delta.aicc", "model.lik", "aicc.weight")
row.names(aicTab.01) <- c("Exponential", "Weibull", "Gompertz", "Logistic")
aicTab.01$k <- c(1, 2, 2, 3)
aicTab.01$log.lik <- c(ml.fit.exp.01$value, ml.fit.weib.01$value, ml.fit.gomp.01$value,  ml.fit.log.01$value)
aicTab.01$deviance <- -2*aicTab.01$log.lik
aicTab.01$aic <- aicTab.01$deviance + 2*aicTab.01$k
aicTab.01$aicc <- aicTab.01$aic + 2*aicTab.01$k*(aicTab.01$k + 1)/(n.01 - aicTab.01$k - 1)
aicTab.01$delta.aicc <- aicTab.01$aicc - min(aicTab.01$aicc)
aicTab.01$model.lik <- exp(0.5*(min(aicTab.01$aicc) - aicTab.01$aicc))
aicTab.01$aicc.weight <- aicTab.01$model.lik/sum(aicTab.01$model.lik)
aicTab.01

#             k   log.lik deviance      aic     aicc delta.aicc    model.lik  aicc.weight
# Exponential 1 -673.4700 1346.940 1348.940 1348.966  314.46339 5.189702e-69 5.189702e-69
# Weibull     2 -537.3478 1074.696 1078.696 1078.773   44.27064 2.436417e-10 2.436417e-10
# Gompertz    2 -561.6462 1123.292 1127.292 1127.370   92.86756 6.824424e-21 6.824424e-21
# Logistic    3 -514.1733 1028.347 1034.347 1034.502    0.00000 1.000000e+00 1.000000e+00
#
# Conclusion: Logistic is best



##### 9B.   Fit survival models for light intensity treatment '1/4' #####

# control for optim function to get max log-likelihood instead of min
C <- list(fnscale = -1) 

# times to death
t.04 <- life.dist.04 

# log-likelihood functions for exponential, Weibull, Gompertz, and logistic models 
loglik.exp.04  <- function(par) sum(log(par[1]*exp(-par[1]*t.04)))
loglik.weib.04 <- function(par) sum(log(par[1]^par[2]*par[2]*t.04^(par[2] - 1)*exp(-(par[1]*t.04)^par[2])))
loglik.gomp.04 <- function(par) sum(log(par[1]*exp(par[2]*t.04 - (par[1]/par[2])*(exp(par[2]*t.04) - 1) )))
loglik.log.04  <- function(par) sum(log(par[1]*exp(par[2]*t.04)*(1 + (par[1]*par[3]/par[2])*(exp(par[2]*t.04) - 1) )^(-(par[3] + 1)/par[3])))

# maximize log-likelihood functions for each of the four survival models
# use 'Brent' method for exponential;  default 'Nelder-Mead' method is not suitable for single param optimization
ml.fit.exp.04  <- optim(c(0.1), loglik.exp.04, control = C, method = "Brent", lower = 0, upper = 1)
ml.fit.weib.04 <- optim(c(0.1, 0.1), loglik.weib.04, control = C)
ml.fit.gomp.04 <- optim(c(0.0001, 0.1), loglik.gomp.04, control = C)
ml.fit.log.04  <- optim(c(0.0001, 0.9, 2), loglik.log.04, control = C)

# extract parameter estimates (parameters a, b, and c, if applicable)
exp.a.04  <- ml.fit.exp.04$par[1]
weib.a.04 <- ml.fit.weib.04$par[1]; weib.b.04 <- ml.fit.weib.04$par[2]
gomp.a.04 <- ml.fit.gomp.04$par[1]; gomp.b.04 <- ml.fit.gomp.04$par[2]
log.a.04  <- ml.fit.log.04$par[1];   log.b.04 <- ml.fit.log.04$par[2]; log.c.04 <- ml.fit.log.04$par[3]

# specify best-fit curves
exp.surv.curve.04 <- function(t.04) exp(-exp.a.04*t.04)
weib.surv.curve.04 <- function(t.04) exp(-(weib.a.04*t.04)^weib.b.04)
gomp.surv.curve.04 <- function(t.04) exp(-(gomp.a.04/gomp.b.04)*(exp(gomp.b.04*t.04) - 1))
log.surv.curve.04 <- function(t.04) (1 + log.c.04*log.a.04/log.b.04*(exp(log.b.04*t.04) - 1))^(-1/log.c.04)

# visually assess best-fit curves

tiff(file = "Surviving_1_4.tiff", width = 5, height = 4, units = "in", res = 600, compression = "lzw")

plot(p.surv.cum ~ age, data = lifeTab.04, log = "y", las = 1, pch = 16, xlab = "Age (days)", ylab = "Proportion surviving")
curve(exp.surv.curve.04, add = T, lwd = 2, col = 1)
curve(weib.surv.curve.04, add = T, lwd = 2, col = 2)
curve(gomp.surv.curve.04, add = T, lwd = 2, col = 3)
curve(log.surv.curve.04, add = T, lwd = 2, col = 4)
legend("bottomleft", c("Exponential", "Weibull", "Gompertz", "Logistic"), col = 1:4, bty = "n", lwd = 2)

dev.off()

# compare models with AICc
aicTab.04 <- data.frame(matrix(ncol = 8, nrow = 4))
names(aicTab.04) <- c("k", "log.lik", "deviance", "aic", "aicc", "delta.aicc", "model.lik", "aicc.weight")
row.names(aicTab.04) <- c("Exponential", "Weibull", "Gompertz", "Logistic")
aicTab.04$k <- c(1, 2, 2, 3)
aicTab.04$log.lik <- c(ml.fit.exp.04$value, ml.fit.weib.04$value, ml.fit.gomp.04$value,  ml.fit.log.04$value)
aicTab.04$deviance <- -2*aicTab.04$log.lik
aicTab.04$aic <- aicTab.04$deviance + 2*aicTab.04$k
aicTab.04$aicc <- aicTab.04$aic + 2*aicTab.04$k*(aicTab.04$k + 1)/(n.04 - aicTab.04$k - 1)
aicTab.04$delta.aicc <- aicTab.04$aicc - min(aicTab.04$aicc)
aicTab.04$model.lik <- exp(0.5*(min(aicTab.04$aicc) - aicTab.04$aicc))
aicTab.04$aicc.weight <- aicTab.04$model.lik/sum(aicTab.04$model.lik)
aicTab.04

#            k   log.lik  deviance       aic      aicc delta.aicc    model.lik  aicc.weight
# Exponential 1 -673.9718 1347.9435 1349.9435 1349.9696  392.52418 5.813763e-86 5.813763e-86
# Weibull     2 -514.1901 1028.3803 1032.3803 1032.4592   75.01373 5.140136e-17 5.140136e-17
# Gompertz    2 -538.1365 1076.2731 1080.2731 1080.3520  122.90656 2.047291e-27 2.047291e-27
# Logistic    3 -475.6433  951.2865  957.2865  957.4455    0.00000 1.000000e+00 1.000000e+00
#
# Conclusion: Logistic is best



##### 10.  Model probability of reproduction as a function of age using GEEs #####

# Make aligned binomial reproduction data frames that include frond ID
wolf.repro.aligned.binom.id.01 <- data.frame(id = wolf$id[wolf$light.treatment == "1"], wolf.repro.aligned.binom.01)
wolf.repro.aligned.binom.id.04 <- data.frame(id = wolf$id[wolf$light.treatment == "1/4"], wolf.repro.aligned.binom.04)

# Convert reproduction data to flat form using melt.data.frame function in package 'reshape'
wolf.repro.aligned.binom.flat.01 <- melt.data.frame(wolf.repro.aligned.binom.id.01, id.vars = "id", variable_name = "age", na.rm = TRUE)
wolf.repro.aligned.binom.flat.04 <- melt.data.frame(wolf.repro.aligned.binom.id.04, id.vars = "id", variable_name = "age", na.rm = TRUE)

# Convert age to numeric ("+ 3" accounts for exclusion of first three days)
wolf.repro.aligned.binom.flat.01$age <- as.numeric(wolf.repro.aligned.binom.flat.01$age)
wolf.repro.aligned.binom.flat.04$age <- as.numeric(wolf.repro.aligned.binom.flat.04$age)

# Sort flat data frames by focal identity
wolf.repro.aligned.binom.flat.01 <-     wolf.repro.aligned.binom.flat.01[    order(wolf.repro.aligned.binom.flat.01$id), ] 
wolf.repro.aligned.binom.flat.04 <-     wolf.repro.aligned.binom.flat.04[    order(wolf.repro.aligned.binom.flat.04$id), ] 

# Convert id to factor so GEE analysis works
wolf.repro.aligned.binom.flat.01$id <-     as.factor(wolf.repro.aligned.binom.flat.01$id)
wolf.repro.aligned.binom.flat.04$id <-     as.factor(wolf.repro.aligned.binom.flat.04$id)

# Perform GEEs using three correlation structures: independence, AR1 (first order autoregressive), and exchangeable 
# For 01 (full light) treatment 
gee1.01 <- geeglm(value ~ age, data=wolf.repro.aligned.binom.flat.01, family=binomial("logit"), id=id, corstr="independence", scale.fix=T)
gee2.01 <- geeglm(value ~ age, data=wolf.repro.aligned.binom.flat.01, family=binomial("logit"), id=id, corstr="ar1", scale.fix=T)
gee3.01 <- geeglm(value ~ age, data=wolf.repro.aligned.binom.flat.01, family=binomial("logit"), id=id, corstr="exchangeable", scale.fix=T)
summary(gee1.01)
summary(gee2.01)
summary(gee3.01)

# For 04 (restricted quarter light) treatment 
gee1.04 <- geeglm(value ~ age, data=wolf.repro.aligned.binom.flat.04, family=binomial("logit"), id=id, corstr="independence", scale.fix=T)
gee2.04 <- geeglm(value ~ age, data=wolf.repro.aligned.binom.flat.04, family=binomial("logit"), id=id, corstr="ar1", scale.fix=T)
gee3.04 <- geeglm(value ~ age, data=wolf.repro.aligned.binom.flat.04, family=binomial("logit"), id=id, corstr="exchangeable", scale.fix=T)
summary(gee1.04)
summary(gee2.04)
summary(gee3.04)

# check for positive definite working correlation matrices
eigen(gee1.01$geese$vbeta.naiv)$values
eigen(gee2.01$geese$vbeta.naiv)$values
eigen(gee3.01$geese$vbeta.naiv)$values

eigen(gee1.04$geese$vbeta.naiv)$values
eigen(gee2.04$geese$vbeta.naiv)$values
eigen(gee3.04$geese$vbeta.naiv)$values # negative value
# fit gee3 using an alternate function ('gee' function in package 'gee')
# still fails to achieve positive definite working correlation matrix
#gee3b.04 <- gee(value ~ age, data=wolf.repro.aligned.binom.flat.04, family=binomial("logit"), id=id, corstr="exchangeable", scale.fix=T)


# calculate RJ values for each GEE to determine best correlation structure
# described in Wang and Carey 2004 (J. American Statistical Association 99:845-853), and Shults et al. 2009 (Statistics in Medicine 28:2338-2355)
Q.gee1.01 <- solve(gee1.01$geese$vbeta.naiv) %*% gee1.01$geese$vbeta
Q.gee2.01 <- solve(gee2.01$geese$vbeta.naiv) %*% gee2.01$geese$vbeta
Q.gee3.01 <- solve(gee3.01$geese$vbeta.naiv) %*% gee3.01$geese$vbeta

Q.gee1.04 <- solve(gee1.04$geese$vbeta.naiv) %*% gee1.04$geese$vbeta
Q.gee2.04 <- solve(gee2.04$geese$vbeta.naiv) %*% gee2.04$geese$vbeta
Q.gee3.04 <- solve(gee3.04$geese$vbeta.naiv) %*% gee3.04$geese$vbeta

RJ1.gee1.01 <- sum(diag(Q.gee1.01)) / 2; RJ1.gee2.01 <- sum(diag(Q.gee2.01)) / 2; RJ1.gee3.01 <- sum(diag(Q.gee3.01)) / 2
RJ1.gee1.04 <- sum(diag(Q.gee1.04)) / 2; RJ1.gee2.04 <- sum(diag(Q.gee2.04)) / 2; RJ1.gee3.04 <- sum(diag(Q.gee3.04)) / 2

RJ2.gee1.01 <- sum(diag(Q.gee1.01 %*% Q.gee1.01)) / 2; RJ2.gee2.01 <- sum(diag(Q.gee2.01 %*% Q.gee2.01)) / 2; RJ2.gee3.01 <- sum(diag(Q.gee3.01 %*% Q.gee3.01)) / 2
RJ2.gee1.04 <- sum(diag(Q.gee1.04 %*% Q.gee1.04)) / 2; RJ2.gee2.04 <- sum(diag(Q.gee2.04 %*% Q.gee2.04)) / 2; RJ2.gee3.04 <- sum(diag(Q.gee3.04 %*% Q.gee3.04)) / 2

RJ3.gee1.01 <- sum((eigen(Q.gee1.01)$values-1)^2) / 2; RJ3.gee2.01 <- sum((eigen(Q.gee2.01)$values-1)^2) / 2; RJ3.gee3.01 <- sum((eigen(Q.gee3.01)$values-1)^2) / 2
RJ3.gee1.04 <- sum((eigen(Q.gee1.04)$values-1)^2) / 2; RJ3.gee2.04 <- sum((eigen(Q.gee2.04)$values-1)^2) / 2; RJ3.gee3.04 <- sum((eigen(Q.gee3.04)$values-1)^2) / 2

# alternate method of calculating RJ3
RJ3alt.gee1.01 <- RJ2.gee1.01 - 2*RJ1.gee1.01 + 1; RJ3alt.gee2.01 <- RJ2.gee2.01 - 2*RJ1.gee2.01 + 1; RJ3alt.gee3.01 <- RJ2.gee3.01 - 2*RJ1.gee3.01 + 1
RJ3alt.gee1.04 <- RJ2.gee1.04 - 2*RJ1.gee1.04 + 1; RJ3alt.gee2.04 <- RJ2.gee2.04 - 2*RJ1.gee2.04 + 1; RJ3alt.gee3.04 <- RJ2.gee3.04 - 2*RJ1.gee3.04 + 1

# summary of RJ values
rjTab.01 <- cbind(c(RJ1.gee1.01, RJ2.gee1.01, RJ3.gee1.01, RJ3alt.gee1.01), c(RJ1.gee2.01, RJ2.gee2.01, RJ3.gee2.01, RJ3alt.gee2.01), c(RJ1.gee3.01, RJ2.gee3.01, RJ3.gee3.01, RJ3alt.gee3.01))
colnames(rjTab.01) <- c("gee1 (independence)", "gee2 (AR1)", "gee3 (exchangeable)")
rownames(rjTab.01) <- c("RJ1", "RJ2", "RJ3", "RJ3 Alt")
rjTab.01 # 01 - ambiguous - RJ1 is closest to 1 for gee2/AR1, but RJ2 is closest to 1 for gee1/independence, and RJ3 is closest to 0 for gee1/independence as well
#        gee1 (independence) gee2 (AR1) gee3 (exchangeable)
#RJ1                   0.710      1.145                 289
#RJ2                   0.813      2.112              166640
#RJ3                   0.392      0.821              166064
#RJ3 Alt               0.392      0.821              166064

rjTab.04 <- cbind(c(RJ1.gee1.04, RJ2.gee1.04, RJ3.gee1.04, RJ3alt.gee1.04), c(RJ1.gee2.04, RJ2.gee2.04, RJ3.gee2.04, RJ3alt.gee2.04))
colnames(rjTab.04) <- c("gee1 (independence)", "gee2 (AR1)")
rownames(rjTab.04) <- c("RJ1", "RJ2", "RJ3", "RJ3 Alt")
rjTab.04 # 04 - gee1/independence is best (RJ1&2 closest to 1, and RJ3 closest to 0)
#        gee1 (independence) gee2 (AR1)
#RJ1                   0.747       2.00
#RJ2                   0.927       6.51
#RJ3                   0.434       3.50
#RJ3 Alt               0.434       3.50

##### 10.1 Excluding the first 2 days of reproduction #####
##### Repeat GEEs excluding first 2 days #####

# Make data frame excluding first 2 days (Names NOT changed)
wolf.repro.aligned.binom.id.01 <- data.frame(id = wolf$id[wolf$light.treatment == "1"], wolf.repro.aligned.binom.01[,c(-1,-2)])
wolf.repro.aligned.binom.id.04 <- data.frame(id = wolf$id[wolf$light.treatment == "1/4"], wolf.repro.aligned.binom.04[,c(-1,-2)])

# Convert reproduction data to flat form using melt.data.frame function in package 'reshape'
wolf.repro.aligned.binom.flat.01 <- melt.data.frame(wolf.repro.aligned.binom.id.01, id.vars = "id", variable_name = "age", na.rm = TRUE)
wolf.repro.aligned.binom.flat.04 <- melt.data.frame(wolf.repro.aligned.binom.id.04, id.vars = "id", variable_name = "age", na.rm = TRUE)

# Convert age to numeric # Exclusion was done here: ("+ 2" accounts for exclusion of first two days)
wolf.repro.aligned.binom.flat.01$age <- as.numeric(wolf.repro.aligned.binom.flat.01$age)+2
wolf.repro.aligned.binom.flat.04$age <- as.numeric(wolf.repro.aligned.binom.flat.04$age)+2

# Sort flat data frames by focal identity
wolf.repro.aligned.binom.flat.01 <-     wolf.repro.aligned.binom.flat.01[    order(wolf.repro.aligned.binom.flat.01$id), ] 
wolf.repro.aligned.binom.flat.04 <-     wolf.repro.aligned.binom.flat.04[    order(wolf.repro.aligned.binom.flat.04$id), ] 

# Convert id to factor so GEE analysis works
wolf.repro.aligned.binom.flat.01$id <-     as.factor(wolf.repro.aligned.binom.flat.01$id)
wolf.repro.aligned.binom.flat.04$id <-     as.factor(wolf.repro.aligned.binom.flat.04$id)

# Perform GEEs using three correlation structures: independence, AR1 (first order autoregressive), and exchangeable 
# For 01 (full light) treatment 
gee1.01 <- geeglm(value ~ age, data=wolf.repro.aligned.binom.flat.01, family=binomial("logit"), id=id, corstr="independence", scale.fix=T)
gee2.01 <- geeglm(value ~ age, data=wolf.repro.aligned.binom.flat.01, family=binomial("logit"), id=id, corstr="ar1", scale.fix=T)
gee3.01 <- geeglm(value ~ age, data=wolf.repro.aligned.binom.flat.01, family=binomial("logit"), id=id, corstr="exchangeable", scale.fix=T)
summary(gee1.01)
summary(gee2.01)
summary(gee3.01)

#> summary(gee1.01)
#Call:
# geeglm(formula = value ~ age, family = binomial("logit"), 
#        data = wolf.repro.aligned.binom.flat.01, id = id, corstr = "independence", 
#        scale.fix = T)

#Coefficients:
#  Estimate  Std.err Wald Pr(>|W|)    
#(Intercept)  0.89110  0.05895  229   <2e-16 ***
#  age         -0.04681  0.00461  103   <2e-16 ***
#  Signif. codes:  0 â***â 0.001 â**â 0.01 â*â 0.05 â.â 0.1 â â 1
#Correlation structure = independence 
#Scale is fixed.
#Number of clusters:   158  Maximum cluster size: 47 

# For 04 (restricted quarter light) treatment 
gee1.04 <- geeglm(value ~ age, data=wolf.repro.aligned.binom.flat.04, family=binomial("logit"), id=id, corstr="independence", scale.fix=T)
gee2.04 <- geeglm(value ~ age, data=wolf.repro.aligned.binom.flat.04, family=binomial("logit"), id=id, corstr="ar1", scale.fix=T)
gee3.04 <- geeglm(value ~ age, data=wolf.repro.aligned.binom.flat.04, family=binomial("logit"), id=id, corstr="exchangeable", scale.fix=T)
summary(gee1.04)
summary(gee2.04)
summary(gee3.04)

#> summary(gee1.04)
#Call:
#  geeglm(formula = value ~ age, family = binomial("logit"), 
#        data = wolf.repro.aligned.binom.flat.04, id = id, corstr = "independence", 
#        scale.fix = T)

#Coefficients:
# Estimate  Std.err Wald Pr(>|W|)    
#(Intercept)  0.52162  0.04309  147   <2e-16 ***
#  age         -0.03667  0.00315  135   <2e-16 ***
# Signif. codes:  0 â***â 0.001 â**â 0.01 â*â 0.05 â.â 0.1 â â 1

#Correlation structure = independence 
#Scale is fixed.
#Number of clusters:   155  Maximum cluster size: 47


# check for positive definite working correlation matrices
eigen(gee1.01$geese$vbeta.naiv)$values
eigen(gee2.01$geese$vbeta.naiv)$values
eigen(gee3.01$geese$vbeta.naiv)$values

eigen(gee1.04$geese$vbeta.naiv)$values
eigen(gee2.04$geese$vbeta.naiv)$values
eigen(gee3.04$geese$vbeta.naiv)$values #negative value

# calculate RJ values for each GEE to determine best correlation structure
# described in Wang and Carey 2004 (J. American Statistical Association 99:845-853), and Shults et al. 2009 (Statistics in Medicine 28:2338-2355)
Q.gee1.01 <- solve(gee1.01$geese$vbeta.naiv) %*% gee1.01$geese$vbeta
Q.gee2.01 <- solve(gee2.01$geese$vbeta.naiv) %*% gee2.01$geese$vbeta
Q.gee3.01 <- solve(gee3.01$geese$vbeta.naiv) %*% gee3.01$geese$vbeta

Q.gee1.04 <- solve(gee1.04$geese$vbeta.naiv) %*% gee1.04$geese$vbeta
Q.gee2.04 <- solve(gee2.04$geese$vbeta.naiv) %*% gee2.04$geese$vbeta
Q.gee3.04 <- solve(gee3.04$geese$vbeta.naiv) %*% gee3.04$geese$vbeta

RJ1.gee1.01 <- sum(diag(Q.gee1.01)) / 2; RJ1.gee2.01 <- sum(diag(Q.gee2.01)) / 2; RJ1.gee3.01 <- sum(diag(Q.gee3.01)) / 2
RJ1.gee1.04 <- sum(diag(Q.gee1.04)) / 2; RJ1.gee2.04 <- sum(diag(Q.gee2.04)) / 2; RJ1.gee3.04 <- sum(diag(Q.gee3.04)) / 2

RJ2.gee1.01 <- sum(diag(Q.gee1.01 %*% Q.gee1.01)) / 2; RJ2.gee2.01 <- sum(diag(Q.gee2.01 %*% Q.gee2.01)) / 2; RJ2.gee3.01 <- sum(diag(Q.gee3.01 %*% Q.gee3.01)) / 2
RJ2.gee1.04 <- sum(diag(Q.gee1.04 %*% Q.gee1.04)) / 2; RJ2.gee2.04 <- sum(diag(Q.gee2.04 %*% Q.gee2.04)) / 2; RJ2.gee3.04 <- sum(diag(Q.gee3.04 %*% Q.gee3.04)) / 2

RJ3.gee1.01 <- sum((eigen(Q.gee1.01)$values-1)^2) / 2; RJ3.gee2.01 <- sum((eigen(Q.gee2.01)$values-1)^2) / 2; RJ3.gee3.01 <- sum((eigen(Q.gee3.01)$values-1)^2) / 2
RJ3.gee1.04 <- sum((eigen(Q.gee1.04)$values-1)^2) / 2; RJ3.gee2.04 <- sum((eigen(Q.gee2.04)$values-1)^2) / 2; RJ3.gee3.04 <- sum((eigen(Q.gee3.04)$values-1)^2) / 2

# alternate method of calculating RJ3
RJ3alt.gee1.01 <- RJ2.gee1.01 - 2*RJ1.gee1.01 + 1; RJ3alt.gee2.01 <- RJ2.gee2.01 - 2*RJ1.gee2.01 + 1; RJ3alt.gee3.01 <- RJ2.gee3.01 - 2*RJ1.gee3.01 + 1
RJ3alt.gee1.04 <- RJ2.gee1.04 - 2*RJ1.gee1.04 + 1; RJ3alt.gee2.04 <- RJ2.gee2.04 - 2*RJ1.gee2.04 + 1; RJ3alt.gee3.04 <- RJ2.gee3.04 - 2*RJ1.gee3.04 + 1

# summary of RJ values
rjTab.01 <- cbind(c(RJ1.gee1.01, RJ2.gee1.01, RJ3.gee1.01, RJ3alt.gee1.01), c(RJ1.gee2.01, RJ2.gee2.01, RJ3.gee2.01, RJ3alt.gee2.01), c(RJ1.gee3.01, RJ2.gee3.01, RJ3.gee3.01, RJ3alt.gee3.01))
colnames(rjTab.01) <- c("gee1 (independence)", "gee2 (AR1)", "gee3 (exchangeable)")
rownames(rjTab.01) <- c("RJ1", "RJ2", "RJ3", "RJ3 Alt")
rjTab.01 # 01 - gee1 is best (closest to 1 for RJ1 & 2, closest to 0 for RJ3)
#         gee1 (independence) gee2 (AR1) gee3 (exchangeable)
#RJ1                   0.911       1.91                4.45
#RJ2                   1.329       6.07               36.13
#RJ3                   0.507       3.25               28.23
#RJ3 Alt               0.507       3.25               28.23

rjTab.04 <- cbind(c(RJ1.gee1.04, RJ2.gee1.04, RJ3.gee1.04, RJ3alt.gee1.04), c(RJ1.gee2.04, RJ2.gee2.04, RJ3.gee2.04, RJ3alt.gee2.04))
colnames(rjTab.04) <- c("gee1 (independence)", "gee2 (AR1)")
rownames(rjTab.04) <- c("RJ1", "RJ2", "RJ3", "RJ3 Alt")
rjTab.04 # 04 - gee1 (independence) is best (closest to 1 for RJ1 &2, closest to 0 for RJ3)
#         gee1 (independence) gee2 (AR1) 
#RJ1                   0.532       1.72               
#RJ2                   0.415       4.45               
#RJ3                   0.351       2.00               
#RJ3 Alt               0.351       2.00                 


##### 11. Plot reproduction versus age with fitted GEE curves #####

tiff(file = "Reproduction_1.tiff", width = 5, height = 4, units = "in", res = 600, compression = "lzw")

int.gee1.01 <- coef(gee1.01)[1]; age.gee1.01 <- coef(gee1.01)[2]
int.gee2.01 <- coef(gee2.01)[1]; age.gee2.01 <- coef(gee2.01)[2]
int.gee3.01 <- coef(gee3.01)[1]; age.gee3.01 <- coef(gee3.01)[2]
rep.curve.gee1.01 <- function(t.01) exp(int.gee1.01+age.gee1.01*t.01)/(1+exp(int.gee1.01+age.gee1.01*t.01))
rep.curve.gee2.01 <- function(t.01) exp(int.gee2.01+age.gee2.01*t.01)/(1+exp(int.gee2.01+age.gee2.01*t.01))
rep.curve.gee3.01 <- function(t.01) exp(int.gee3.01+age.gee3.01*t.01)/(1+exp(int.gee3.01+age.gee3.01*t.01))
plot(p.repro ~ age, data=lifeTab.01, main = "Reproduction full light (01)")
curve(rep.curve.gee1.01, add=TRUE, lwd=2, lty=1, col=1)
curve(rep.curve.gee2.01, add=TRUE, lwd=2, lty=1, col=2)
curve(rep.curve.gee3.01, add=TRUE, lwd=2, lty=1, col=3)
legend("top", c("Independence", "AR1", "Exchangeable"), col=1:3, bty="n", lwd=2)

dev.off()

tiff(file = "Reproduction_1_4.tiff", width = 5, height = 4, units = "in", res = 600, compression = "lzw")

int.gee1.04 <- coef(gee1.04)[1]; age.gee1.04 <- coef(gee1.04)[2]
int.gee2.04 <- coef(gee2.04)[1]; age.gee2.04 <- coef(gee2.04)[2]
rep.curve.gee1.04 <- function(t.04) exp(int.gee1.04+age.gee1.04*t.04)/(1+exp(int.gee1.04+age.gee1.04*t.04))
rep.curve.gee2.04 <- function(t.04) exp(int.gee2.04+age.gee2.04*t.04)/(1+exp(int.gee2.04+age.gee2.04*t.04))
plot(p.repro ~ age, data=lifeTab.04, main = "Reproduction restricted quarter light (04)")
curve(rep.curve.gee1.04, add=TRUE, lwd=2, lty=1, col=1)
curve(rep.curve.gee2.04, add=TRUE, lwd=2, lty=1, col=2)
legend("top", c("Independence", "AR1"), col=1:2, bty="n", lwd=2)

dev.off()


##### Combined survivorship and reproduction of two treatments #####
#survivorship logistic model both
#reproduction independence both

library(patchwork)
library(ggplot2)
theme_as <- function() {
  theme_bw(base_size = 15) +
    theme(
      legend.position = "none",
      axis.text = element_text(color = "black", size = 10),
      axis.title = element_text(color = "black", size = 14),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      plot.title = element_text(hjust = 0.5, size = 14),
      plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")
    )
}


# Survivorship Plot for Full Light
p.01 <- ggplot(lifeTab.01, aes(x = age, y = p.surv.cum)) +
  geom_point(shape = 1, size = 2, stroke = 1, color = "orange") +  # Changed to orange for full light
  stat_function(fun = log.surv.curve.01, aes(color = "Logistic model"), linewidth = 1.0) +
  scale_y_log10() +
  scale_color_manual(values = c("Logistic model" = "black")) +  # Changed to black for better visibility
  labs(x = "Age (days)", 
       y = "Proportion Surviving", 
       color = "Model",
       title = "Full Light") +
  expand_limits(y = 0) +
  theme_as() 


# Survivorship Plot for Dim Light
p.04 <- ggplot(lifeTab.04, aes(x = age, y = p.surv.cum)) +
  geom_point(shape = 1, size = 2, stroke = 1, color = "dark blue") +  # Changed to dark blue for dim light
  stat_function(fun = log.surv.curve.04, aes(color = "Logistic model"), linewidth = 1.0) +
  scale_y_log10() +
  scale_color_manual(values = c("Logistic model" = "black")) +  # Changed to black for better visibility
  labs(x = "Age (days)", 
       y = "Proportion Surviving", 
       color = "Model",
       title = "Dim Light") +
  expand_limits(y = 0) +
  theme_as()+
  theme(axis.title.y = element_blank())  # Remove y-axis label

# Reproduction Plot for Full Light
p.repro01 <- ggplot(lifeTab.01, aes(x = age, y = p.repro)) +
  geom_point(shape = 1, size = 2, stroke = 1, color = "orange") +
  stat_function(fun = rep.curve.gee1.01, color = "black", linewidth = 1.0, xlim = c(3, 50)) + #starts at Day3
  labs(title = "Full Light", x = "Age (days)", y = "Proportion Reproducing") +
  theme_as() +
  scale_y_continuous(limits = c(0, 1))

# Reproduction Plot for Dim Light
p.repro04 <- ggplot(lifeTab.04, aes(x = age, y = p.repro)) +
  geom_point(shape = 1, size = 2, stroke = 1, color = "dark blue") +
  stat_function(fun = rep.curve.gee1.04, color = "black", linewidth = 1.0, xlim = c(3, 50)) +  #starts at Day3
  labs(title = "Dim Light", x = "Age (days)", y = "Proportion Reproducing") +  
  theme_as() +
  scale_y_continuous(limits = c(0, 1)) +
  theme(axis.title.y = element_blank())

# Combine the plots vertically
combined_plot_sur_repro <- (p.01 + p.04) / (p.repro01 + p.repro04) +
  plot_layout(heights = c(1, 1)) +
  plot_annotation(tag_levels = list(c("(a)", "", "(b)", ""))) &
  theme(plot.tag = element_text(size = 14, face = "bold"))

# Display the combined plot
print(combined_plot_sur_repro)
# Save the combined plot
ggsave("Wolf_Combined_Sur_Repro_Curves_excludedfirst2_test_20250612.png", combined_plot_sur_repro, width = 12, height = 8, units = "in", dpi = 600)



##### the first age at which the cumulative proportion of survivors drops below 50% #####
age_50_percent_01 <- lifeTab.01$age[which.min(abs(lifeTab.01$p.surv.cum - 0.5))]
age_50_percent_04 <- lifeTab.04$age[which.min(abs(lifeTab.04$p.surv.cum - 0.5))]

surv_prop_01 <- lifeTab.01$p.surv.cum[lifeTab.01$age == age_50_percent_01]
surv_prop_04 <- lifeTab.04$p.surv.cum[lifeTab.04$age == age_50_percent_04]

print(paste("Survival proportion at age", age_50_percent_01, "(full light):", round(surv_prop_01, 3)))
#"Survival proportion at age  (full light): 0.475
print(paste("Survival proportion at age", age_50_percent_04, "(dim light):", round(surv_prop_04, 3)))
#"Survival proportion at age  (dim light): 0.555

### Full-light see  the probabilities at specific ages
prob_age1.01 <- rep.curve.gee1.01(3)
print(prob_age1.01)
#0.679
prob_agemax.01 <- rep.curve.gee1.01(49)
print(prob_agemax.01)
#0.197

### Dim-light see  the probabilities at specific ages
prob_age1.04 <- rep.curve.gee1.04(3)
print(prob_age1.04)
#0.601
prob_agemax.04 <- rep.curve.gee1.04(49)
print(prob_agemax.04)
#0.218


##### 12. Komologorov-Smirnov test on residual log lifespan data --Stroustrop et al. (2016) #####

# Get the lifespan distribution of each treatment
life.dist.2.01 <- wolf$lifespan[wolf$light.treatment == "1"] # distribution of lifespan
life.dist.2.04 <- wolf$lifespan[wolf$light.treatment == "1/4"]

# Get the distributions of log lifespan
log.life.dist.2.01 <- log(life.dist.2.01)
log.life.dist.2.04 <- log(life.dist.2.04)

# Get the residuals of mean log lifespan
resid.log.life.dist.01 <- log.life.dist.2.01 - mean(log.life.dist.2.01)
resid.log.life.dist.04 <- log.life.dist.2.04 - mean(log.life.dist.2.04)

# Perform Komolgorov-Smirnov test

set.seed(12345)
suppressWarnings(ks.test(resid.log.life.dist.01, resid.log.life.dist.04)) # Suppress warnings; we know about ties

#Two-sample Kolmogorov-Smirnov test
#data:  resid.log.life.dist.01 and resid.log.life.dist.04
#D = 0.13434, p-value = 0.1187
#alternative hypothesis: two-sided


##### To visualize this step; run step 1-3, 8A, 8B, 12.

##### 13.   Load required packages and verify existing data #####

# Load packages 
library(ggplot2)
library(gridExtra)
library(grid)

##### 14.   Use existing survival data from life tables #####

# Use existing life tables (lifeTab.01 and lifeTab.04) that were already calculated
# These contain: age, p.surv.cum (cumulative survival proportion)

# Create survival data for absolute age using existing life tables
survival_absolute <- rbind(
  data.frame(
    time = lifeTab.01$age,
    survival = lifeTab.01$p.surv.cum,
    treatment = "Full light"
  ),
  data.frame(
    time = lifeTab.04$age,
    survival = lifeTab.04$p.surv.cum,
    treatment = "Dim light"
  )
)

##### 15.   Calculate survival for residual log lifespans #####

# Use existing residual log lifespan data from section 12.1
# resid.log.life.dist.01 and resid.log.life.dist.04 are already calculated

# Function to calculate survival curve from raw lifespan data
calculate_survival_from_lifespans <- function(lifespans) {
  lifespans <- lifespans[!is.na(lifespans)]  # Remove NA values
  lifespans <- sort(lifespans)  # Sort in ascending order
  n_total <- length(lifespans)
  unique_times <- sort(unique(lifespans))
  
  survival_data <- data.frame(
    time = unique_times,
    survival = sapply(unique_times, function(t) {
      n_surviving <- sum(lifespans >= t)  # Number still alive at time t
      return(n_surviving / n_total)  # Proportion surviving
    })
  )
  
  return(survival_data)
}

# Calculate survival curves for residual log lifespans
surv_resid_01 <- calculate_survival_from_lifespans(resid.log.life.dist.01)
surv_resid_04 <- calculate_survival_from_lifespans(resid.log.life.dist.04)

# Combine residual survival data
survival_residual <- rbind(
  data.frame(surv_resid_01, treatment = "Full light"),
  data.frame(surv_resid_04, treatment = "Dim light")
)

##### 16.   Create temporal scaling plots using existing data #####

# Function to create survival plots
create_survival_plot <- function(survival_data, x_label, y_label = NULL, show_legend = TRUE) {
  p <- ggplot(survival_data, aes(x = time, y = survival, color = treatment)) +
    geom_step(linewidth = 1, direction = "hv") +  # Step function for survival curves
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2), expand = c(0, 0)) +
    scale_color_manual(
      values = c("Full light" = "orange", "Dim light" = "dark blue"),
      breaks = c("Full light", "Dim light")
    ) +
    labs(x = x_label, y = y_label) +
    theme_bw() +
    theme(
      legend.position = if(show_legend) c(0.80, 0.85) else "none",
      legend.background = element_rect(fill = "white", color = "black"),
      legend.title = element_blank(),
      legend.text = element_text(size = 9),
      axis.text = element_text(color = "black", size = 10),
      axis.title = element_text(color = "black", size = 14),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
    ) +
    coord_cartesian(expand = FALSE)
  
  return(p)
}

# Create survival plot for absolute age (using existing life table data)
p1 <- create_survival_plot(
  survival_absolute, 
  x_label = "Absolute age (days)", 
  y_label = "Proportion surviving",
  show_legend = TRUE
)

# Create survival plot for residual log lifespans
p2 <- create_survival_plot(
  survival_residual, 
  x_label = "Residual log lifespans", 
  y_label = "",
  show_legend = TRUE
)

# Combine plots into single figure
combined_plot <- grid.arrange(
  arrangeGrob(p1, top = textGrob("(a)", x = unit(0.02, "npc"), just = "left", 
                                 gp = gpar(fontsize = 14, face = "bold"))),
  arrangeGrob(p2, top = textGrob("(b)", x = unit(0.02, "npc"), just = "left", 
                                 gp = gpar(fontsize = 14, face = "bold"))),
  ncol = 2
)

# Display the combined plot
print(combined_plot)

##### 17.   Verify Kolmogorov-Smirnov test results #####

# Use the existing residual log lifespan data from section 12.1
# Set seed for reproducible results 
set.seed(12345)

# Perform Kolmogorov-Smirnov test on residual log lifespans
ks_test_result <- suppressWarnings(ks.test(resid.log.life.dist.01, resid.log.life.dist.04))
print(ks_test_result)

# Display key statistics
cat("\nKolmogorov-Smirnov Test Results:\n")
cat("D statistic:", round(ks_test_result$statistic, 4), "\n")
cat("p-value:", round(ks_test_result$p.value, 4), "\n")

##### 18.   Save temporal scaling figure #####

# Save the combined temporal scaling plot
ggsave("temporal_scaling_wolffia_residual_1.png", combined_plot, 
       width = 10, height = 5, units = "in", dpi = 600)

##### 19.   Summary of data sources used #####

cat("\nData Sources Summary:\n")
cat("- Absolute age survival: lifeTab.01$age, lifeTab.01$p.surv.cum (Full light)\n")
cat("                        lifeTab.04$age, lifeTab.04$p.surv.cum (Dim light)\n")
cat("- Residual log lifespans: resid.log.life.dist.01 (Full light)\n")
cat("                         resid.log.life.dist.04 (Dim light)\n")
cat("- Sample sizes: n.01 =", n.01, "(Full light), n.04 =", n.04, "(Dim light)\n")

