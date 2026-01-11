##### Lemna aequinoctialis Project
##### Athita Senayai
##### Note: Parts of script modified from Wolffia project

##### PART I. LOAD PACKAGES AND SET WORKING DIRECTORY #####


##### 1.   Load packages #####

library(reshape)
library(geepack)
library(gee)


##### 2.   Set working directory #####

# Set working directory
setwd('/Users/athita/Downloads/OneDrive_2026-01-10/Lemna aequinoctialis')
##### PART II. DATA ANALYSIS: SURVIVAL AND REPRODUCTION #####


##### 3.   Load daily survival and reproduction data, and perform preliminary processing #####

# Load daily reproduction data
Lemna.full <- read.csv("20240312_Lemna_daily_surv_and_reproduction.csv", header = T) 

# Exclude appropriate fronds
Lemna <- Lemna.full[Lemna.full$exclude == "No",]

# Change light treatment names to be less cryptic 
Lemna$light.treatment[Lemna$light.treatment == "3"] <- "1/3"
Lemna$light.treatment <- factor(Lemna$light.treatment, c("1", "1/3"))

# Extract reproduction data 
Lemna.repro <- subset(Lemna, select = Oct.27.2023:Feb.09.2024) 

# Double check dates of birth
Lemna$date.birth.alt <- apply(Lemna.repro, 1, function (x) names(Lemna.repro)[min(which(x == 0))])
Lemna$date.birth.alt == Lemna$date.birth # all entries match
Lemna <- Lemna[ , 1:(ncol(Lemna) - 1)] # drop unnecessary date.birth.alt column

# Double check dates of last reproduction
Lemna$date.last.repro.alt <- apply(Lemna.repro, 1, function (x) names(Lemna.repro)[max(which(x > 0))])
Lemna$date.last.repro.alt == Lemna$date.last.repro # all entries match 
Lemna <- Lemna[ , 1:(ncol(Lemna) - 1)] # drop unnecessary date.last.repro.alt column 

# Calculate dates of first reproduction
Lemna$date.first.repro <- apply(Lemna.repro, 1, function (x) names(Lemna.repro)[min(which(x > 0))])

# Convert date values to R's date format
Lemna$date.birth <- as.Date(Lemna$date.birth, format = "%b.%d.%Y")
Lemna$date.last.repro <- as.Date(Lemna$date.last.repro, format = "%b.%d.%Y")
Lemna$date.first.repro <- as.Date(Lemna$date.first.repro, format = "%b.%d.%Y")

# Calculate lifespan from birth
Lemna$lifespan <- as.numeric(Lemna$date.last.repro - Lemna$date.birth + 1) # first day is Day 1, not Day 0

# Calculate lifetime reproduction
Lemna$total.offspring <- apply(Lemna.repro, 1, function (x) sum(x, na.rm = T))

# Make new data frame with aligned reproduction data relative to dates of birth
Lemna.repro.aligned <- matrix(NA, nrow = nrow(Lemna.repro), ncol = max(Lemna$lifespan))
birth.col <- as.numeric(Lemna$date.birth - min(Lemna$date.birth) + 1) # column index for date of birth       
death.col <- as.numeric(Lemna$date.last.repro - min(Lemna$date.birth) + 1) # column index for date of death
for (i in 1:nrow(Lemna.repro)) {
  Lemna.repro.aligned[i, 1:Lemna$lifespan[i]] <- as.numeric(Lemna.repro[i, birth.col[i]:death.col[i]])
}

# Calculate intrinsic rate of increase of individuals (r)
for (i in 1:nrow(Lemna.repro.aligned)) {
  ind <- Lemna.repro.aligned[i, ]                                # extract individual in row i
  L <- max(which(ind > 0))                                     # get lifespan
  lesMat <- matrix(0, nrow = L, ncol = L)               # create empty lifespan x lifespan Leslie matrix
  lesMat[1, ] <- ind[1:L]                                      # fill first row with reproduction data
  if (L > 1) {
    lesMat[2:L, 1:(L - 1)] <- diag(rep(1, (L - 1)))            # fill sub-diagonal with survival data
  }
  eigenvalues <- eigen(lesMat)$values                          # get list of eigenvalues of Leslie matrix
  lambda <- max(Re(eigenvalues[abs(Im(eigenvalues)) < 1e-6]))  # extract lambda
  Lemna$r[i] <- log(lambda)                                      # record intrinsic rate of increase (r)
  if (abs(Lemna$r[i]) < 1e-6) {
    Lemna$r[i] <- 0                               # correct for floating point arithmetic rounding errors
  }
} 

##### 4.  Analyze reproductive lifespan (Mann-Whitney) #####

# Visually inspect boxplot of lifespan
tiff(file = "Lemna_lifespan_1_3.tiff", width = 5, height = 4, units = "in", res = 600, compression = "lzw")

boxplot(lifespan ~ light.treatment, dat = Lemna, ylim = c(0, 70),
        las = 1, whisklty = 1, col = "white", xlab = "Light intensity treatment", ylab = "Lifespan (days)")

dev.off()

### Summary statistical values in two treatments ###
by(Lemna$lifespan, Lemna$light.treatment, summary)

#Lemna$light.treatment: 1
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#8.00   15.00   24.00   24.74   32.00   62.00

#Lemna$light.treatment: 1/3
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#11.00   25.75   33.00   33.67   39.25   87.00 


### Data Dianostics divided by treatment ###
shapiro_test <- by(Lemna$lifespan, Lemna$light.treatment, shapiro.test)
shapiro_test # Not normally distributed (p< 0.05)

#Lemna$light.treatment: 1
#W = 0.95627, p-value = 4.204e-05

#Lemna$light.treatment: 1/3
#W = 0.94825, p-value = 7.931e-06

treatment_groups <- split(Lemna$lifespan, Lemna$light.treatment)
par(mfrow=c(1,2)) # Set up a 1x2 grid for plotting histograms side by side
for (i in 1:length(treatment_groups)) {
  hist(treatment_groups[[i]], main=paste("Treatment", names(treatment_groups)[i], "Histogram"),
       xlab="Plant Lifespan", ylab="Frequency", col="lightblue", border="white")
}

variance_test <- var.test(lifespan ~ light.treatment, data = Lemna)
variance_test #Equal variance  (p>0.05)

#data:  lifespan by light.treatment
#F = 0.92147, num df = 167, denom df = 167, p-value = 0.5978

# Assumption of normality unable to be met for t-test. Non-parametric Mann-Whitney test used instead.
# Mann-Whitney test on lifespan
wilcox.test(lifespan ~ light.treatment, data = Lemna) #differrence between the groups

#Wilcoxon rank sum test with continuity correction
#data:  lifespan by light.treatment
#W = 8152.5, p-value = 2.139e-11
#alternative hypothesis: true location shift is not equal to 0



##### 5.   Analyze lifetime reproduction (Mann-Whitney) #####

# Visually inspect boxplot of lifetime reproduction
tiff(file = "Lemna_Offspring_1_3.tiff", width = 5, height = 4, units = "in", res = 600, compression = "lzw")

boxplot(total.offspring ~ light.treatment, dat = Lemna, ylim = c(0, 30), las = 1, whisklty = 1, col = "white", xlab = "Light intensity treatment", ylab = "Total number of offspring")

dev.off()

### Summary statistical values in two treatments ###
by(Lemna$total.offspring, Lemna$light.treatment, summary)

#Lemna$light.treatment: 1
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#3.00    8.00   10.00   10.37   13.00   23.00 

#Lemna$light.treatment: 1/3
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#3.00    8.00   10.00   10.43   12.00   24.00 

### Data Dianostics divided by treatment ###
shapiro_test <- by(Lemna$total.offspring, Lemna$light.treatment, shapiro.test)
shapiro_test # Not normally distributed (p< 0.05)

#Lemna$light.treatment: 1
#W = 0.96835, p-value = 0.0006982

#Lemna$light.treatment: 1/3
#W = 0.91199, p-value = 1.625e-08


treatment_groups <- split(Lemna$total.offspring, Lemna$light.treatment)
par(mfrow=c(1,2)) # Set up a 1x2 grid for plotting histograms side by side
for (i in 1:length(treatment_groups)) {
  hist(treatment_groups[[i]], main=paste("Treatment", names(treatment_groups)[i], "Histogram"),
       xlab="Plant Lifespan", ylab="Frequency", col="lightblue", border="white")
}


variance_test <- var.test(total.offspring ~ light.treatment, data = Lemna)
variance_test #Equal variance  (p>0.05)

#data:  total.offspring by light.treatment
#F = 1.108, num df = 167, denom df = 167, p-value = 0.5084

# Assumption of normality unable to be met for t-test. Non-parametric Mann-Whitney test used instead.
# Mann-Whitney test on total.offspring

wilcox.test(total.offspring ~ light.treatment, data = Lemna) #no difference between the group
#Wilcoxon rank sum test with continuity correction
#data:  total.offspring by light.treatment
#W = 14196, p-value = 0.9254
#alternative hypothesis: true location shift is not equal to 0




##### 6.   Analyze intrinsic rate of increase of individuals (r) (t-test) #####

# Visually inspect boxplot of intrinsic rate of increase of individuals (r)

tiff(file = "Intrinsic_1_3.tiff", width = 5, height = 4, units = "in", res = 600, compression = "lzw")

boxplot(r ~ light.treatment, dat = Lemna, ylim = c(0, 0.6), las = 1,  whisklty = 1, col = "white", 
        xlab = "Light intensity treatment", ylab = "Intrinsic rate of increase (r)")

dev.off()


### Summary statistical values in two treatments ###
by(Lemna$r, Lemna$light.treatment, summary)

#Lemna$light.treatment: 1
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.1589  0.3147  0.3531  0.3493  0.3869  0.5159 

#Lemna$light.treatment: 1/3
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.1086  0.2745  0.3197  0.3113  0.3570  0.4389 


### Data Dianostics divided by treatment ###
shapiro_test <- by(Lemna$r, Lemna$light.treatment, shapiro.test)
shapiro_test # Not normally distributed (p< 0.05)

#Lemna$light.treatment: 1
#W = 0.98804, p-value = 0.1646 # Normally distributed (p> 0.05)

#Lemna$light.treatment: 1/3
#W = 0.9698, p-value = 0.001008  # Not normally distributed (p< 0.05)

treatment_groups <- split(Lemna$r, Lemna$light.treatment)
par(mfrow=c(1,2)) # Set up a 1x2 grid for plotting histograms side by side
for (i in 1:length(treatment_groups)) {
  hist(treatment_groups[[i]], main=paste("Treatment", names(treatment_groups)[i], "Histogram"),
       xlab="Plant Lifespan", ylab="Frequency", col="lightblue", border="white")
}

variance_test <- var.test(r ~ light.treatment, data = Lemna)
variance_test # Equal variance  (p>0.05)

#data:  r by light.treatment
#F = 0.79511, num df = 167, denom df = 167, p-value = 0.1395


# Assumption of normality unable to be met for t-test. Non-parametric Mann-Whitney test used instead.
# Mann-Whitney test on r
wilcox.test(r ~ light.treatment, data = Lemna)

#Wilcoxon rank sum test with continuity correction

#data:  r by light.treatment
#W = 18839, p-value = 1e-07
#alternative hypothesis: true location shift is not equal to 0


##### Combined three figures: Lifespan, total number of offspring, intrinsic rate of increase of individuals (r)

library(ggplot2)
library(gridExtra)
library(grid)

Lemna$light.treatment <- factor(Lemna$light.treatment,
                                levels = c("1", "1/3"),
                                labels = c("Full light", "Dim light"))


create_plot <- function(Lemna, y_var, y_lab) {
  ggplot(Lemna, aes(x = light.treatment, y = {{y_var}})) +
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
p1 <- create_plot(Lemna, lifespan, "Lifespan (days)")
p2 <- create_plot(Lemna, total.offspring, "Total number of offspring")
p3 <- create_plot(Lemna, r, expression("Intrinsic rate of increase (" * italic("r") * ")"))

# Combine plots using grid.arrange and add labels using grid.text
combined_plot <- grid.arrange(
  arrangeGrob(p1, top = textGrob("(a)", x = unit(0, "npc"), just = "left", gp = gpar(fontsize = 12, face = "bold"))),
  arrangeGrob(p2, top = textGrob("(b)", x = unit(0, "npc"), just = "left", gp = gpar(fontsize = 12, face = "bold"))),
  arrangeGrob(p3, top = textGrob("(c)", x = unit(0, "npc"), just = "left", gp = gpar(fontsize = 12, face = "bold"))),
  ncol = 1, 
  bottom = textGrob("Light intensity treatment", gp = gpar(fontsize = 12)),
  top = textGrob("Lemna Under Different Light Intensities", gp = gpar(fontsize = 12)),
  heights = c(1, 1, 1)
)

# Display the combined plot
print(combined_plot)
ggsave("Lemna_combined_plots_senesence_20250314.png", combined_plot, width = 5, height = 10, units = "in", dpi = 600)



##### 7.   Create aligned reproduction data for each light intensity treatment #####

### IMPORTANT: If you encounter the error "no non-missing arguments to max; returning -Inf" 
### when running this section, it means there's a factor level mismatch because the 
### light.treatment levels were changed in the plotting section (around Line 264).
### 
### EASY FIX: Clear your workspace and re-run steps 1-3, then run this step (Step 7) 
### before running any plotting sections.
###
### Alternative: Update the factor levels in this section to match Line 264


# Create separate aligned reproduction data for each light intensity treatment
Lemna.repro.aligned.01 <- Lemna.repro.aligned[Lemna$light.treatment == "1", ]
Lemna.repro.aligned.03 <- Lemna.repro.aligned[Lemna$light.treatment == "1/3", ]

# Convert aligned reproduction data into binomial format
Lemna.repro.aligned.binom.01 <- apply(Lemna.repro.aligned.01, 2, function (x) replace(x, which(x > 1), 1))
Lemna.repro.aligned.binom.03 <- apply(Lemna.repro.aligned.03, 2, function (x) replace(x, which(x > 1), 1))



##### 8A.   Create life tables for light intensity treatment '1' #####

n.01 <- nrow(Lemna.repro.aligned.01) # total starting sample size
life.dist.01 <- Lemna$lifespan[Lemna$light.treatment == "1"] # distribution of lifespan
maxAge.01 <- max(life.dist.01) # max lifespan
meanAge.01 <- mean(life.dist.01) # mean lifespan
lifeTab.01 <- data.frame(matrix(nrow = maxAge.01, ncol = 7))
names(lifeTab.01) <- c("age", "age.rel", "n.last.repro", "n.surv.cum", "p.surv.cum", "n.repro", "p.repro")    
lifeTab.01$age <- 1:maxAge.01 # ages
lifeTab.01$age.rel <- lifeTab.01$age/meanAge.01 # relative ages (in mean lifespans)
lifeTab.01$n.last.repro <- tabulate(life.dist.01, nbins = maxAge.01) # number of deaths at each age
lifeTab.01$n.surv.cum <- c(n.01, n.01 - cumsum(lifeTab.01$n.last.repro[1:(maxAge.01 - 1)])) # total cumulative survivorship at age
lifeTab.01$p.surv.cum <- lifeTab.01$n.surv.cum/n.01 # proportional cumulative survivorship at age
lifeTab.01$n.repro <- colSums(Lemna.repro.aligned.binom.01[ , 1:maxAge.01], na.rm = TRUE) # number reproducing at each age
lifeTab.01$p.repro <- lifeTab.01$n.repro/lifeTab.01$n.surv.cum # proportion reproducing at age


##### 8B.   Create life tables for light intensity treatment '1/3' #####

n.03 <- nrow(Lemna.repro.aligned.03) # total starting sample size
life.dist.03 <- Lemna$lifespan[Lemna$light.treatment == "1/3"] # distribution of lifespan
maxAge.03 <- max(life.dist.03) # max lifespan
meanAge.03 <- mean(life.dist.03) # mean lifespan
lifeTab.03 <- data.frame(matrix(nrow = maxAge.03, ncol = 7))
names(lifeTab.03) <- c("age", "age.rel", "n.last.repro", "n.surv.cum", "p.surv.cum", "n.repro", "p.repro")    
lifeTab.03$age <- 1:maxAge.03 # ages
lifeTab.03$age.rel <- lifeTab.03$age/meanAge.03 # relative ages (in mean lifespans)
lifeTab.03$n.last.repro <- tabulate(life.dist.03, nbins = maxAge.03) # number of deaths at each age
lifeTab.03$n.surv.cum <- c(n.03, n.03 - cumsum(lifeTab.03$n.last.repro[1:(maxAge.03 - 1)])) # total cumulative survivorship at age
lifeTab.03$p.surv.cum <- lifeTab.03$n.surv.cum/n.03 # proportional cumulative survivorship at age
lifeTab.03$n.repro <- colSums(Lemna.repro.aligned.binom.03[ , 1:maxAge.03], na.rm = TRUE) # number reproducing at each age
lifeTab.03$p.repro <- lifeTab.03$n.repro/lifeTab.03$n.surv.cum # proportion reproducing at age

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
# Exponential 1 -707.0423 1414.085 1416.085 1416.109  137.34366 1.500386e-30 1.497388e-30
# Weibull     2 -637.3461 1274.692 1278.692 1278.765    0.00000 1.000000e+00 9.980019e-01
# Gompertz    2 -651.0560 1302.112 1306.112 1306.185   27.41981 1.111382e-06 1.109161e-06
# Logistic    3 -642.5235 1285.047 1291.047 1291.193   12.42826 2.000957e-03 1.996958e-03
#
# Conclusion: Weibull is the best


##### 9B.   Fit survival models for light intensity treatment '1/3' #####

# control for optim function to get max log-likelihood instead of min
C <- list(fnscale = -1) 

# times to death
t.03 <- life.dist.03 

# log-likelihood functions for exponential, Weibull, Gompertz, and logistic models 
loglik.exp.03  <- function(par) sum(log(par[1]*exp(-par[1]*t.03)))
loglik.weib.03 <- function(par) sum(log(par[1]^par[2]*par[2]*t.03^(par[2] - 1)*exp(-(par[1]*t.03)^par[2])))
loglik.gomp.03 <- function(par) sum(log(par[1]*exp(par[2]*t.03 - (par[1]/par[2])*(exp(par[2]*t.03) - 1) )))
loglik.log.03  <- function(par) sum(log(par[1]*exp(par[2]*t.03)*(1 + (par[1]*par[3]/par[2])*(exp(par[2]*t.03) - 1) )^(-(par[3] + 1)/par[3])))

# maximize log-likelihood functions for each of the four survival models
# use 'Brent' method for exponential;  default 'Nelder-Mead' method is not suitable for single param optimization
ml.fit.exp.03  <- optim(c(0.1), loglik.exp.03, control = C, method = "Brent", lower = 0, upper = 1)
ml.fit.weib.03 <- optim(c(0.1, 0.1), loglik.weib.03, control = C)
ml.fit.gomp.03 <- optim(c(0.0001, 0.1), loglik.gomp.03, control = C)
ml.fit.log.03  <- optim(c(0.0001, 0.9, 2), loglik.log.03, control = C)

# extract parameter estimates (parameters a, b, and c, if applicable)
exp.a.03  <- ml.fit.exp.03$par[1]
weib.a.03 <- ml.fit.weib.03$par[1]; weib.b.03 <- ml.fit.weib.03$par[2]
gomp.a.03 <- ml.fit.gomp.03$par[1]; gomp.b.03 <- ml.fit.gomp.03$par[2]
log.a.03  <- ml.fit.log.03$par[1];   log.b.03 <- ml.fit.log.03$par[2]; log.c.03 <- ml.fit.log.03$par[3]

# specify best-fit curves
exp.surv.curve.03 <- function(t.03) exp(-exp.a.03*t.03)
weib.surv.curve.03 <- function(t.03) exp(-(weib.a.03*t.03)^weib.b.03)
gomp.surv.curve.03 <- function(t.03) exp(-(gomp.a.03/gomp.b.03)*(exp(gomp.b.03*t.03) - 1))
log.surv.curve.03 <- function(t.03) (1 + log.c.03*log.a.03/log.b.03*(exp(log.b.03*t.03) - 1))^(-1/log.c.03)

# visually assess best-fit curves

tiff(file = "Surviving_1_3.tiff", width = 5, height = 4, units = "in", res = 600, compression = "lzw")

plot(p.surv.cum ~ age, data = lifeTab.03, log = "y", las = 1, pch = 16, xlab = "Age (days)", ylab = "Proportion surviving")
curve(exp.surv.curve.03, add = T, lwd = 2, col = 1)
curve(weib.surv.curve.03, add = T, lwd = 2, col = 2)
curve(gomp.surv.curve.03, add = T, lwd = 2, col = 3)
curve(log.surv.curve.03, add = T, lwd = 2, col = 4)
legend("bottomleft", c("Exponential", "Weibull", "Gompertz", "Logistic"), col = 1:4, bty = "n", lwd = 2)

dev.off()

# compare models with AICc
aicTab.03 <- data.frame(matrix(ncol = 8, nrow = 4))
names(aicTab.03) <- c("k", "log.lik", "deviance", "aic", "aicc", "delta.aicc", "model.lik", "aicc.weight")
row.names(aicTab.03) <- c("Exponential", "Weibull", "Gompertz", "Logistic")
aicTab.03$k <- c(1, 2, 2, 3)
aicTab.03$log.lik <- c(ml.fit.exp.03$value, ml.fit.weib.03$value, ml.fit.gomp.03$value,  ml.fit.log.03$value)
aicTab.03$deviance <- -2*aicTab.03$log.lik
aicTab.03$aic <- aicTab.03$deviance + 2*aicTab.03$k
aicTab.03$aicc <- aicTab.03$aic + 2*aicTab.03$k*(aicTab.03$k + 1)/(n.03 - aicTab.03$k - 1)
aicTab.03$delta.aicc <- aicTab.03$aicc - min(aicTab.03$aicc)
aicTab.03$model.lik <- exp(0.5*(min(aicTab.03$aicc) - aicTab.03$aicc))
aicTab.03$aicc.weight <- aicTab.03$model.lik/sum(aicTab.03$model.lik)
aicTab.03

#            k   log.lik  deviance       aic      aicc delta.aicc    model.lik  aicc.weight
# Exponential 1 -758.7734 1517.547 1519.547 1519.571  223.06501 3.647945e-49 3.625157e-49
# Weibull     2 -651.2860 1302.572 1306.572 1306.645   10.13885 6.286023e-03 6.246755e-03
# Gompertz    2 -678.6427 1357.285 1361.285 1361.358   64.85219 8.270402e-15 8.218739e-15
# Logistic    3 -645.1798 1290.360 1296.360 1296.506    0.00000 1.000000e+00 9.937532e-01
#
# Conclusion: Logistic is the best



##### 10.  Model probability of reproduction as a function of age using GEEs #####

# Make aligned binomial reproduction data frames that include frond ID
Lemna.repro.aligned.binom.id.01 <- data.frame(id = Lemna$id[Lemna$light.treatment == "1"], Lemna.repro.aligned.binom.01)
Lemna.repro.aligned.binom.id.03 <- data.frame(id = Lemna$id[Lemna$light.treatment == "1/3"], Lemna.repro.aligned.binom.03)

# Convert reproduction data to flat form using melt.data.frame function in package 'reshape'
Lemna.repro.aligned.binom.flat.01 <- melt.data.frame(Lemna.repro.aligned.binom.id.01, id.vars = "id", variable_name = "age", na.rm = TRUE)
Lemna.repro.aligned.binom.flat.03 <- melt.data.frame(Lemna.repro.aligned.binom.id.03, id.vars = "id", variable_name = "age", na.rm = TRUE)

# Convert age to numeric
Lemna.repro.aligned.binom.flat.01$age <- as.numeric(Lemna.repro.aligned.binom.flat.01$age)
Lemna.repro.aligned.binom.flat.03$age <- as.numeric(Lemna.repro.aligned.binom.flat.03$age)

# Sort flat data frames by focal identity
Lemna.repro.aligned.binom.flat.01 <-     Lemna.repro.aligned.binom.flat.01[    order(Lemna.repro.aligned.binom.flat.01$id), ] 
Lemna.repro.aligned.binom.flat.03 <-     Lemna.repro.aligned.binom.flat.03[    order(Lemna.repro.aligned.binom.flat.03$id), ] 

# Convert id to factor so GEE analysis works
Lemna.repro.aligned.binom.flat.01$id <-     as.factor(Lemna.repro.aligned.binom.flat.01$id)
Lemna.repro.aligned.binom.flat.03$id <-     as.factor(Lemna.repro.aligned.binom.flat.03$id)

# Perform GEEs using three correlation structures: independence, AR1 (first order autoregressive), and exchangeable 
# For 01 (full light) treatment 
gee1.01 <- geeglm(value ~ age, data=Lemna.repro.aligned.binom.flat.01, family=binomial("logit"), id=id, corstr="independence", scale.fix=T)
gee2.01 <- geeglm(value ~ age, data=Lemna.repro.aligned.binom.flat.01, family=binomial("logit"), id=id, corstr="ar1", scale.fix=T)
gee3.01 <- geeglm(value ~ age, data=Lemna.repro.aligned.binom.flat.01, family=binomial("logit"), id=id, corstr="exchangeable", scale.fix=T)
summary(gee1.01)
summary(gee2.01)
summary(gee3.01)

#Exchageable is the best model 
#Coefficients:
#           Estimate Std.err  Wald Pr(>|W|)    
#(Intercept)  -0.0368  0.0497  0.55     0.46    
#age          -0.0343  0.0036 90.48   <2e-16 ***
#Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

#Correlation structure = exchangeable 
#Scale is fixed.

#Link = identity 

#Estimated Correlation Parameters:
#      Estimate Std.err
#alpha   0.0161 0.00611
#Number of clusters:   168  Maximum cluster size: 62 


# For 03 (restricted quarter light) treatment 
gee1.03 <- geeglm(value ~ age, data=Lemna.repro.aligned.binom.flat.03, family=binomial("logit"), id=id, corstr="independence", scale.fix=T)
gee2.03 <- geeglm(value ~ age, data=Lemna.repro.aligned.binom.flat.03, family=binomial("logit"), id=id, corstr="ar1", scale.fix=T)
gee3.03 <- geeglm(value ~ age, data=Lemna.repro.aligned.binom.flat.03, family=binomial("logit"), id=id, corstr="exchangeable", scale.fix=T)
summary(gee1.03)
summary(gee2.03)
summary(gee3.03)

#Exchageable is the best model 

#Coefficients:
#           Estimate  Std.err  Wald Pr(>|W|)    
#(Intercept) -0.49298  0.04373 127.1   <2e-16 ***
#age         -0.02525  0.00278  82.5   <2e-16 ***
#Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

#Correlation structure = exchangeable 
#Scale is fixed.

# Link = identity 

#Estimated Correlation Parameters:
#      Estimate Std.err
#alpha  0.00357 0.00361
#Number of clusters:   168  Maximum cluster size: 87

# check for positive definite working correlation matrices
eigen(gee1.01$geese$vbeta.naiv)$values
eigen(gee2.01$geese$vbeta.naiv)$values
eigen(gee3.01$geese$vbeta.naiv)$values

eigen(gee1.03$geese$vbeta.naiv)$values
eigen(gee2.03$geese$vbeta.naiv)$values
eigen(gee3.03$geese$vbeta.naiv)$values 

# calculate RJ values for each GEE to determine best correlation structure
# described in Wang and Carey 2004 (J. American Statistical Association 99:845-853), and Shults et al. 2009 (Statistics in Medicine 28:2338-2355)
Q.gee1.01 <- solve(gee1.01$geese$vbeta.naiv) %*% gee1.01$geese$vbeta
Q.gee2.01 <- solve(gee2.01$geese$vbeta.naiv) %*% gee2.01$geese$vbeta
Q.gee3.01 <- solve(gee3.01$geese$vbeta.naiv) %*% gee3.01$geese$vbeta

Q.gee1.03 <- solve(gee1.03$geese$vbeta.naiv) %*% gee1.03$geese$vbeta
Q.gee2.03 <- solve(gee2.03$geese$vbeta.naiv) %*% gee2.03$geese$vbeta
Q.gee3.03 <- solve(gee3.03$geese$vbeta.naiv) %*% gee3.03$geese$vbeta

RJ1.gee1.01 <- sum(diag(Q.gee1.01)) / 2; RJ1.gee2.01 <- sum(diag(Q.gee2.01)) / 2; RJ1.gee3.01 <- sum(diag(Q.gee3.01)) / 2
RJ1.gee1.03 <- sum(diag(Q.gee1.03)) / 2; RJ1.gee2.03 <- sum(diag(Q.gee2.03)) / 2; RJ1.gee3.03 <- sum(diag(Q.gee3.03)) / 2

RJ2.gee1.01 <- sum(diag(Q.gee1.01 %*% Q.gee1.01)) / 2; RJ2.gee2.01 <- sum(diag(Q.gee2.01 %*% Q.gee2.01)) / 2; RJ2.gee3.01 <- sum(diag(Q.gee3.01 %*% Q.gee3.01)) / 2
RJ2.gee1.03 <- sum(diag(Q.gee1.03 %*% Q.gee1.03)) / 2; RJ2.gee2.03 <- sum(diag(Q.gee2.03 %*% Q.gee2.03)) / 2; RJ2.gee3.03 <- sum(diag(Q.gee3.03 %*% Q.gee3.03)) / 2

RJ3.gee1.01 <- sum((eigen(Q.gee1.01)$values-1)^2) / 2; RJ3.gee2.01 <- sum((eigen(Q.gee2.01)$values-1)^2) / 2; RJ3.gee3.01 <- sum((eigen(Q.gee3.01)$values-1)^2) / 2
RJ3.gee1.03 <- sum((eigen(Q.gee1.03)$values-1)^2) / 2; RJ3.gee2.03 <- sum((eigen(Q.gee2.03)$values-1)^2) / 2; RJ3.gee3.03 <- sum((eigen(Q.gee3.03)$values-1)^2) / 2

# alternate method of calculating RJ3
RJ3alt.gee1.01 <- RJ2.gee1.01 - 2*RJ1.gee1.01 + 1; RJ3alt.gee2.01 <- RJ2.gee2.01 - 2*RJ1.gee2.01 + 1; RJ3alt.gee3.01 <- RJ2.gee3.01 - 2*RJ1.gee3.01 + 1
RJ3alt.gee1.03 <- RJ2.gee1.03 - 2*RJ1.gee1.03 + 1; RJ3alt.gee2.03 <- RJ2.gee2.03 - 2*RJ1.gee2.03 + 1; RJ3alt.gee3.03 <- RJ2.gee3.03 - 2*RJ1.gee3.03 + 1

# summary of RJ values
rjTab.01 <- cbind(c(RJ1.gee1.01, RJ2.gee1.01, RJ3.gee1.01, RJ3alt.gee1.01), c(RJ1.gee2.01, RJ2.gee2.01, RJ3.gee2.01, RJ3alt.gee2.01), c(RJ1.gee3.01, RJ2.gee3.01, RJ3.gee3.01, RJ3alt.gee3.01))
colnames(rjTab.01) <- c("gee1 (independence)", "gee2 (AR1)", "gee3 (exchangeable)")
rownames(rjTab.01) <- c("RJ1", "RJ2", "RJ3", "RJ3 Alt")
rjTab.01 

#           gee1 (independence) gee2 (AR1) gee3  (exchangeable)
#RJ1                   1.373      1.559               1.033
#RJ2                   2.216      2.855               1.206
#RJ3                   0.469      0.736               0.141
#RJ3 Alt               0.469      0.736               0.141

rjTab.03 <- cbind(c(RJ1.gee1.03, RJ2.gee1.03, RJ3.gee1.03, RJ3alt.gee1.03), c(RJ1.gee2.03, RJ2.gee2.03, RJ3.gee2.03, RJ3alt.gee2.03), c(RJ1.gee3.03, RJ2.gee3.03, RJ3.gee3.03, RJ3alt.gee3.03))
colnames(rjTab.03) <- c("gee1 (independence)", "gee2 (AR1)", "gee3 (exchangeable)")
rownames(rjTab.03) <- c("RJ1", "RJ2", "RJ3", "RJ3 Alt")
rjTab.03 

#         gee1 (independence) gee2 (AR1) gee3 (exchangeable)
#RJ1                   1.163      1.369               1.056
#RJ2                   1.664      2.316               1.350
#RJ3                   0.338      0.579               0.239
#RJ3 Alt               0.338      0.579               0.239


##### 10.1 Excluding the first 2 days of reproduction #####
##### Repeat GEEs excluding first 2 days #####

# Make data frame excluding first 2 days (Names NOT changed)
Lemna.repro.aligned.binom.id.01 <- data.frame(id = Lemna$id[Lemna$light.treatment == "1"], Lemna.repro.aligned.binom.01[,c(-1,-2)])
Lemna.repro.aligned.binom.id.03 <- data.frame(id = Lemna$id[Lemna$light.treatment == "1/3"], Lemna.repro.aligned.binom.03[,c(-1,-2)])

# Convert reproduction data to flat form using melt.data.frame function in package 'reshape'
Lemna.repro.aligned.binom.flat.01 <- melt.data.frame(Lemna.repro.aligned.binom.id.01, id.vars = "id", variable_name = "age", na.rm = TRUE)
Lemna.repro.aligned.binom.flat.03 <- melt.data.frame(Lemna.repro.aligned.binom.id.03, id.vars = "id", variable_name = "age", na.rm = TRUE)

# Convert age to numeric # Exclusion was done here: ("+ 2" accounts for exclusion of first two days)
Lemna.repro.aligned.binom.flat.01$age <- as.numeric(Lemna.repro.aligned.binom.flat.01$age)+2
Lemna.repro.aligned.binom.flat.03$age <- as.numeric(Lemna.repro.aligned.binom.flat.03$age)+2

# Sort flat data frames by focal identity
Lemna.repro.aligned.binom.flat.01 <-     Lemna.repro.aligned.binom.flat.01[    order(Lemna.repro.aligned.binom.flat.01$id), ] 
Lemna.repro.aligned.binom.flat.03 <-     Lemna.repro.aligned.binom.flat.03[    order(Lemna.repro.aligned.binom.flat.03$id), ] 

# Convert id to factor so GEE analysis works
Lemna.repro.aligned.binom.flat.01$id <-     as.factor(Lemna.repro.aligned.binom.flat.01$id)
Lemna.repro.aligned.binom.flat.03$id <-     as.factor(Lemna.repro.aligned.binom.flat.03$id)

# Perform GEEs using three correlation structures: independence, AR1 (first order autoregressive), and exchangeable 

# For 01 (full light) treatment 
gee1.01 <- geeglm(value ~ age, data=Lemna.repro.aligned.binom.flat.01, family=binomial("logit"), id=id, corstr="independence", scale.fix=T)
gee2.01 <- geeglm(value ~ age, data=Lemna.repro.aligned.binom.flat.01, family=binomial("logit"), id=id, corstr="ar1", scale.fix=T)
gee3.01 <- geeglm(value ~ age, data=Lemna.repro.aligned.binom.flat.01, family=binomial("logit"), id=id, corstr="exchangeable", scale.fix=T)
summary(gee1.01)
summary(gee2.01)
summary(gee3.01)

#> summary(gee3.01)
#Call:
#  geeglm(formula = value ~ age, family = binomial("logit"), 
#        data = Lemna.repro.aligned.binom.flat.01, id = id, corstr = "exchangeable", 
#         scale.fix = T)

#Coefficients:
# Estimate  Std.err  Wald Pr(>|W|)    
#(Intercept)  0.44988  0.06608  46.4  9.9e-12 ***
# age         -0.05811  0.00442 172.8  < 2e-16 ***

#  Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#Correlation structure = exchangeable 
#Scale is fixed.

#Link = identity 

#Estimated Correlation Parameters:
#  Estimate Std.err
#alpha   0.0157 0.00625
#Number of clusters:   168  Maximum cluster size: 60

# For 03 (restricted quarter light) treatment 
gee1.03 <- geeglm(value ~ age, data=Lemna.repro.aligned.binom.flat.03, family=binomial("logit"), id=id, corstr="independence", scale.fix=T)
gee2.03 <- geeglm(value ~ age, data=Lemna.repro.aligned.binom.flat.03, family=binomial("logit"), id=id, corstr="ar1", scale.fix=T)
gee3.03 <- geeglm(value ~ age, data=Lemna.repro.aligned.binom.flat.03, family=binomial("logit"), id=id, corstr="exchangeable", scale.fix=T)
summary(gee1.03)
summary(gee2.03)
summary(gee3.03)

#> summary(gee3.03)
#Call:
#  geeglm(formula = value ~ age, family = binomial("logit"), 
#         data = Lemna.repro.aligned.binom.flat.03, id = id, corstr = "exchangeable", 
#         scale.fix = T)

#Coefficients:
#  Estimate  Std.err  Wald Pr(>|W|)    
#(Intercept) -0.21268  0.05453  15.2  9.6e-05 ***
#  age      -0.03621  0.00333 118.1  < 2e-16 ***

#  Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#Correlation structure = exchangeable 
#Scale is fixed.
#Link = identity 
#Estimated Correlation Parameters:
#  Estimate Std.err
#alpha  0.00481 0.00395
#Number of clusters:   168  Maximum cluster size: 85 

# check for positive definite working correlation matrices
eigen(gee1.01$geese$vbeta.naiv)$values
eigen(gee2.01$geese$vbeta.naiv)$values
eigen(gee3.01$geese$vbeta.naiv)$values

eigen(gee1.03$geese$vbeta.naiv)$values
eigen(gee2.03$geese$vbeta.naiv)$values
eigen(gee3.03$geese$vbeta.naiv)$values

# calculate RJ values for each GEE to determine best correlation structure
# described in Wang and Carey 2004 (J. American Statistical Association 99:845-853), and Shults et al. 2009 (Statistics in Medicine 28:2338-2355)
Q.gee1.01 <- solve(gee1.01$geese$vbeta.naiv) %*% gee1.01$geese$vbeta
Q.gee2.01 <- solve(gee2.01$geese$vbeta.naiv) %*% gee2.01$geese$vbeta
Q.gee3.01 <- solve(gee3.01$geese$vbeta.naiv) %*% gee3.01$geese$vbeta

Q.gee1.03 <- solve(gee1.03$geese$vbeta.naiv) %*% gee1.03$geese$vbeta
Q.gee2.03 <- solve(gee2.03$geese$vbeta.naiv) %*% gee2.03$geese$vbeta
Q.gee3.03 <- solve(gee3.03$geese$vbeta.naiv) %*% gee3.03$geese$vbeta

RJ1.gee1.01 <- sum(diag(Q.gee1.01)) / 2; RJ1.gee2.01 <- sum(diag(Q.gee2.01)) / 2; RJ1.gee3.01 <- sum(diag(Q.gee3.01)) / 2
RJ1.gee1.03 <- sum(diag(Q.gee1.03)) / 2; RJ1.gee2.03 <- sum(diag(Q.gee2.03)) / 2; RJ1.gee3.03 <- sum(diag(Q.gee3.03)) / 2

RJ2.gee1.01 <- sum(diag(Q.gee1.01 %*% Q.gee1.01)) / 2; RJ2.gee2.01 <- sum(diag(Q.gee2.01 %*% Q.gee2.01)) / 2; RJ2.gee3.01 <- sum(diag(Q.gee3.01 %*% Q.gee3.01)) / 2
RJ2.gee1.03 <- sum(diag(Q.gee1.03 %*% Q.gee1.03)) / 2; RJ2.gee2.03 <- sum(diag(Q.gee2.03 %*% Q.gee2.03)) / 2; RJ2.gee3.03 <- sum(diag(Q.gee3.03 %*% Q.gee3.03)) / 2

RJ3.gee1.01 <- sum((eigen(Q.gee1.01)$values-1)^2) / 2; RJ3.gee2.01 <- sum((eigen(Q.gee2.01)$values-1)^2) / 2; RJ3.gee3.01 <- sum((eigen(Q.gee3.01)$values-1)^2) / 2
RJ3.gee1.03 <- sum((eigen(Q.gee1.03)$values-1)^2) / 2; RJ3.gee2.03 <- sum((eigen(Q.gee2.03)$values-1)^2) / 2; RJ3.gee3.03 <- sum((eigen(Q.gee3.03)$values-1)^2) / 2

# alternate method of calculating RJ3
RJ3alt.gee1.01 <- RJ2.gee1.01 - 2*RJ1.gee1.01 + 1; RJ3alt.gee2.01 <- RJ2.gee2.01 - 2*RJ1.gee2.01 + 1; RJ3alt.gee3.01 <- RJ2.gee3.01 - 2*RJ1.gee3.01 + 1
RJ3alt.gee1.03 <- RJ2.gee1.03 - 2*RJ1.gee1.03 + 1; RJ3alt.gee2.03 <- RJ2.gee2.03 - 2*RJ1.gee2.03 + 1; RJ3alt.gee3.03 <- RJ2.gee3.03 - 2*RJ1.gee3.03 + 1

# summary of RJ values
rjTab.01 <- cbind(c(RJ1.gee1.01, RJ2.gee1.01, RJ3.gee1.01, RJ3alt.gee1.01), c(RJ1.gee2.01, RJ2.gee2.01, RJ3.gee2.01, RJ3alt.gee2.01), c(RJ1.gee3.01, RJ2.gee3.01, RJ3.gee3.01, RJ3alt.gee3.01))
colnames(rjTab.01) <- c("gee1 (independence)", "gee2 (AR1)", "gee3 (exchangeable)")
rownames(rjTab.01) <- c("RJ1", "RJ2", "RJ3", "RJ3 Alt")
rjTab.01 
#gee1              (independence) gee2 (AR1) gee3 (exchangeable)
#RJ1                   1.405      1.670              1.0947
#RJ2                   2.179      3.098              1.2753
#RJ3                   0.368      0.758              0.0859
#RJ3 Alt               0.368      0.758              0.0859
#Exchangeble is the best

rjTab.03 <- cbind(c(RJ1.gee1.03, RJ2.gee1.03, RJ3.gee1.03, RJ3alt.gee1.03), c(RJ1.gee2.03, RJ2.gee2.03, RJ3.gee2.03, RJ3alt.gee2.03), c(RJ1.gee3.03, RJ2.gee3.03, RJ3.gee3.03, RJ3alt.gee3.03))
colnames(rjTab.03) <- c("gee1 (independence)", "gee2 (AR1)", "gee3 (exchangeable)")
rownames(rjTab.03) <- c("RJ1", "RJ2", "RJ3", "RJ3 Alt")
rjTab.03
#gee1              (independence) gee2 (AR1) gee3 (exchangeable)
#RJ1                    1.26      1.510               1.128
#RJ2                    1.94      2.796               1.526
#RJ3                    0.41      0.777               0.269
#RJ3 Alt                0.41      0.777               0.269
#exchangeable is the best


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

tiff(file = "Reproduction_1_3.tiff", width = 5, height = 4, units = "in", res = 600, compression = "lzw")

int.gee1.03 <- coef(gee1.03)[1]; age.gee1.03 <- coef(gee1.03)[2]
int.gee2.03 <- coef(gee2.03)[1]; age.gee2.03 <- coef(gee2.03)[2]
int.gee3.03 <- coef(gee3.03)[1]; age.gee3.03 <- coef(gee3.03)[2]
rep.curve.gee1.03 <- function(t.03) exp(int.gee1.03+age.gee1.03*t.03)/(1+exp(int.gee1.03+age.gee1.03*t.03))
rep.curve.gee2.03 <- function(t.03) exp(int.gee2.03+age.gee2.03*t.03)/(1+exp(int.gee2.03+age.gee2.03*t.03))
rep.curve.gee3.03 <- function(t.03) exp(int.gee3.03+age.gee3.03*t.03)/(1+exp(int.gee3.03+age.gee3.03*t.03))
plot(p.repro ~ age, data=lifeTab.03, main = "Reproduction dim light (03)")
curve(rep.curve.gee1.03, add=TRUE, lwd=2, lty=1, col=1)
curve(rep.curve.gee2.03, add=TRUE, lwd=2, lty=1, col=2)
curve(rep.curve.gee3.03, add=TRUE, lwd=2, lty=1, col=3)
legend("top", c("Independence", "AR1", "Exchangeable"), col=1:3, bty="n", lwd=2)


dev.off()

##### Combined survivorship and reproduction of two treatments #####
#survivorship: Weibull model --full-light, logistic model --dim-light 
#reproduction exchangeable both

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
  geom_point(shape = 25, size = 1.5, stroke = 1, color = "orange") +  # Changed to yellow for full light
  stat_function(fun = weib.surv.curve.01, aes(color = "Weibull model"), linewidth = 1.0) +
  scale_y_log10() +
  scale_color_manual(values = c("Weibull model" = "black")) +  # Changed to orange for better visibility
  labs(x = "Age (days)", 
       y = "Proportion Surviving", 
       color = "Model",
       title = "Full Light") +
  expand_limits(y = 0) +
  theme_as() 


# Survivorship Plot for Dim Light
p.03 <- ggplot(lifeTab.03, aes(x = age, y = p.surv.cum)) +
  geom_point(shape = 25, size = 1.5, stroke = 1, color = "dark blue") +  # Changed to dark blue for dim light
  stat_function(fun = log.surv.curve.03, aes(color = "Logistic model"), linewidth = 1.0) +
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
  geom_point(shape = 25, size = 1.5, stroke = 1, color = "orange") +
  stat_function(fun = rep.curve.gee1.01, color = "black", linewidth = 1.0, xlim = c(3, 62)) + #starts at Day3
  labs(title = "Full Light", x = "Age (days)", y = "Proportion Reproducing") +
  theme_as() +
  scale_y_continuous(limits = c(0, 1))

# Reproduction Plot for Dim Light
p.repro03 <- ggplot(lifeTab.03, aes(x = age, y = p.repro)) +
  geom_point(shape = 25, size = 1.5, stroke = 1, color = "dark blue") +
  stat_function(fun = rep.curve.gee1.03, color = "black", linewidth = 1.0, xlim = c(3, 87)) +  #starts at Day3
  labs(title = "Dim Light", x = "Age (days)", y = "Proportion Reproducing") +  
  theme_as() +
  scale_y_continuous(limits = c(0, 1)) +
  theme(axis.title.y = element_blank())


# Combine the plots vertically
combined_plot_sur_repro <- (p.01 + p.03) / (p.repro01 + p.repro03) +
  plot_layout(heights = c(1, 1)) +
  plot_annotation(tag_levels = list(c("(a)", "", "(b)", ""))) &
  theme(plot.tag = element_text(size = 14, face = "bold"))

# Display the combined plot
print(combined_plot_sur_repro)
# Save the combined plot
ggsave("Lemna_Combined_Sur_Repro_Curves_excludedfirst2_test20250612.png", combined_plot_sur_repro, width = 12, height = 8, units = "in", dpi = 600)



##### the first age at which the cumulative proportion of survivors drops below 50% #####
age_50_percent_01 <- lifeTab.01$age[which.min(abs(lifeTab.01$p.surv.cum - 0.5))]
age_50_percent_03 <- lifeTab.03$age[which.min(abs(lifeTab.03$p.surv.cum - 0.5))]

surv_prop_01 <- lifeTab.01$p.surv.cum[lifeTab.01$age == age_50_percent_01]
surv_prop_03 <- lifeTab.03$p.surv.cum[lifeTab.03$age == age_50_percent_03]

print(paste("Survival proportion at age", age_50_percent_01, "(full light):", round(surv_prop_01, 3)))
#"Survival proportion at age 24 (full light): 0.53
print(paste("Survival proportion at age", age_50_percent_03, "(dim light):", round(surv_prop_03, 3)))
#"Survival proportion at age 33 (dim light): 0.506"

### Full-light see  the probabilities at specific ages
prob_age1.01 <- rep.curve.gee3.01(3)
print(prob_age1.01)
#0.568
prob_agemax.01 <- rep.curve.gee3.01(62)
print(prob_agemax.01)
#0.041

### Dim-light see  the probabilities at specific ages
prob_age1.03 <- rep.curve.gee3.03(3)
print(prob_age1.03)
#0.42
prob_agemax.03 <- rep.curve.gee3.03(87)
print(prob_agemax.03)
#0.0335 


##### 12. Komologorov-Smirnov test on residual log lifespan data --Stroustrop et al. (2016) #####

# Get the lifespan distribution of each treatment
life.dist.2.01 <- Lemna$lifespan[Lemna$light.treatment == "1"] # distribution of lifespan
life.dist.2.03 <- Lemna$lifespan[Lemna$light.treatment == "1/3"]


# Get the distributions of log lifespan
log.life.dist.2.01 <- log(life.dist.2.01)
log.life.dist.2.03 <- log(life.dist.2.03)

# Get the residuals of mean log lifespan
resid.log.life.dist.01 <- log.life.dist.2.01 - mean(log.life.dist.2.01)
resid.log.life.dist.03 <- log.life.dist.2.03 - mean(log.life.dist.2.03)

# Perform Komolgorov-Smirnov test

set.seed(12345)
suppressWarnings(ks.test(resid.log.life.dist.01, resid.log.life.dist.03)) # Suppress warnings; we know about ties

#data:  resid.log.life.dist.01 and resid.log.life.dist.03
#D = 0.14286, p-value = 0.06486

##### To visualize this step; run step 1-3, 8A, 8B, 12. 

##### 13.   Load required packages and verify existing data #####

# Load packages 
library(ggplot2)
library(gridExtra)
library(grid)

##### 14.   Use existing survival data from life tables #####

# Use existing life tables (lifeTab.01 and lifeTab.03) that were already calculated
# These contain: age, p.surv.cum (cumulative survival proportion)

# Create survival data for absolute age using existing life tables
survival_absolute <- rbind(
  data.frame(
    time = lifeTab.01$age,
    survival = lifeTab.01$p.surv.cum,
    treatment = "Full light"
  ),
  data.frame(
    time = lifeTab.03$age,
    survival = lifeTab.03$p.surv.cum,
    treatment = "Dim light"
  )
)

##### 15.   Calculate survival for residual log lifespans #####

# Use existing residual log lifespan data from section 12.1
# resid.log.life.dist.01 and resid.log.life.dist.03 are already calculated

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
surv_resid_03 <- calculate_survival_from_lifespans(resid.log.life.dist.03)

# Combine residual survival data
survival_residual <- rbind(
  data.frame(surv_resid_01, treatment = "Full light"),
  data.frame(surv_resid_03, treatment = "Dim light")
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

# Use the existing residual log lifespan data from section 12.
# Set seed for reproducible results 
set.seed(12345)

# Perform Kolmogorov-Smirnov test on residual log lifespans
ks_test_result <- suppressWarnings(ks.test(resid.log.life.dist.01, resid.log.life.dist.03))
print(ks_test_result)

# Display key statistics
cat("\nKolmogorov-Smirnov Test Results:\n")
cat("D statistic:", round(ks_test_result$statistic, 4), "\n")
cat("p-value:", round(ks_test_result$p.value, 4), "\n")

##### 18.   Save temporal scaling figure #####

# Save the combined temporal scaling plot
ggsave("temporal_scaling_lemna_residual.png", combined_plot, 
       width = 10, height = 5, units = "in", dpi = 600)

##### 19.   Summary of data sources used #####

cat("\nData Sources Summary:\n")
cat("- Absolute age survival: lifeTab.01$age, lifeTab.01$p.surv.cum (Full light)\n")
cat("                        lifeTab.03$age, lifeTab.03$p.surv.cum (Dim light)\n")
cat("- Residual log lifespans: resid.log.life.dist.01 (Full light)\n")
cat("                         resid.log.life.dist.03 (Dim light)\n")
cat("- Sample sizes: n.01 =", n.01, "(Full light), n.03 =", n.03, "(Dim light)\n")

