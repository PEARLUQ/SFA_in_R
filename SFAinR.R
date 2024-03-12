rm(list=ls())
graphics.off()

# Data Process for the Queensland hospital data for illustrations------------

# Read data
data <- read.csv("QLD.csv")
names(data)[names(data)=="HOSID"] <- "id"
names(data)[names(data)=="Yeardummy"] <- "Year"

# Convert to panel data
library(plm)
paneldata<- pdata.frame(data, c("id","Year"))
# Define formulas for different models
form = lAggout ~ lBEDS + lAgglabours + lSUPP
formz = lAggout ~ lBEDS + lAgglabours + lSUPP | TEACH + Remote + Small
formt = lAggout ~ lBEDS + lAgglabours + lSUPP + factor(id)
t = data$Year
t2 = t^2
formc = lAggout ~ lBEDS + lAgglabours + lSUPP + t + t2

# Define factor variables and formulas for SVKZ
attach (data)
fTEACH = as.factor(TEACH)
fRemote = as.factor(Remote)
fSmall = as.factor(Small)
forms = lAggout ~ lBEDS + lAgglabours + lSUPP + fTEACH + fRemote + fSmall
forms3 = ehat3 ~ lBEDS + lAgglabours + lSUPP + fTEACH + fRemote + fSmall

# Define output variable
Y = lAggout

################################################################# 
### The following code can be applied with the synthetic data ###
###     or the Queensland hospital data processed above.      ###
################################################################# 

# Code for SFMs ------------
library(frontier)
########################## ALS77 ##########################
attach(data)
# Model estimation
als77 <- sfa (form, data = data, ineffDecrease = T, truncNorm = F, timeEffect = F)
summary(als77)
# Individual inefficiency by JLMS
  # Coefficients
  als77coef <- coef(als77, which = "mle", extraPar = T)
  # epsilon from fitted values
  fals77 <- fitted(als77, asInData = T)
  ei = Y - fals77
  # E(ui|ei)
  us2 = (als77coef[["sigmaU"]])^2
  vs2 = (als77coef[["sigmaV"]])^2
  sigmastar = sqrt((vs2*us2)/(vs2+us2))
  ustari = (-us2*ei)/(vs2+us2)
  uals77 = ((sigmastar*dnorm(ustari/sigmastar))/(pnorm(ustari/sigmastar)))+ustari
  inals77 = 1 - exp(-uals77) # Inefficiency estimation of individual unit
  summary(inals77)

########################## SS84 ##########################
attach(paneldata)
# Fixed effects for ui
ss84 <- plm(form, data = paneldata,
            model = "within", index = c("id","Year"), effect = "individual")
summary(ss84)
  # ai and uihat
  ai = as.numeric(unname(fixef(ss84)))
  uss84 = max(ai) - ai
  inss84 = (uss84 - min(uss84))/(max(uss84)-min(uss84)) # Relative inefficiency of individual unit
  summary(inss84)

########################## PL81 ##########################
attach(paneldata)
# Fixed effects
pl81f <- plm(form, data = paneldata,
            model = "random", index = c("id","Year"), effect = "individual")
summary(pl81f)
  # Extract initial values
  sigmau2 = as.numeric(ercomp(pl81f)[["sigma2"]][2])
  sigmav2 = as.numeric(ercomp(pl81f)[["sigma2"]][1])
  sigmasq = sigmau2+sigmav2
  gamma = sigmau2/sigmasq

  init_para = c(as.numeric(pl81f$coefficients), sigmasq, gamma)
# MLE with initial values
pl81 <- sfa (form, data = paneldata, ineffDecrease = T, truncNorm = F, timeEffect =F,
             startVal = init_para)
summary(pl81)

  # Estimate individual inefficiency
  inpl81 <- 1-efficiencies (pl81, asInData = T, logDepVar = T, minusU = T)
  summary(inpl81)

########################## CSS90 ##########################
attach(paneldata)
# Fixed effects for ui
c90ori <- plm(formc, data = paneldata,
                model = "within", index = c("HOSID","Yeardummy"), effect = "twoways")
summary(c90ori)
ait = lAggout - lBEDS*as.numeric(c90ori[["coefficients"]])[1] - lAgglabours*as.numeric(c90ori[["coefficients"]])[2] - lSUPP*as.numeric(c90ori[["coefficients"]])[3]

  # Calculate ajt
  ti = c(1, 2, 3, 4)
  ti2 = ti^2
  aith = numeric()
  
  for (i in 1:length(unique(id))){
    aith = c(aith, as.numeric(unname(
      lm(ait[(4*(i-1)+1):(4*(i-1)+4)] ~ ti2 + ti)$fitted.values)))
  }
  
  c90 = cbind(aith, t)
  # athat for each year
  ajtyear = aggregate(aith ~ t, data = c90, max)
  # Compute uithat
  c90 = as.data.frame(merge(c90, ajtyear, by = "t"))
  colnames(c90) <- c("t","ait","ajt")
  uc90 = c90$ajt - c90$ait
  summary(uc90)
  inc90 = (uc90 - min(uc90))/(max(uc90)-min(uc90)) # Relative inefficiency of individual unit
  summary(inc90)
  
########################## BC92 ##########################
attach(paneldata)
# Model estimation
bc92 <- sfa (form, data = paneldata, ineffDecrease = T, truncNorm = T, timeEffect = T)
summary(bc92, effic = F, logDepVar = T, effMinusU = T)
  
  # Estimate individual inefficiency
  inbc92 <- 1-efficiencies (bc92, asInData = T, logDepVar = T, minusU = T)
  summary(inbc92)
  
########################## G05 ##########################
attach(paneldata)
# Model estimation
g05 <- sfa (formt, data = paneldata, ineffDecrease = T, truncNorm = F, timeEffect = T)
summary(g05, effic = F, logDepVar = T, effMinusU = T)

  # Estimate individual inefficiency
  ing05 <- 1-efficiencies (g05, asInData = T, logDepVar = T, minusU = T)
  summary(ing05)

########################## KLH14 ##########################
attach(paneldata)
# Fixed effects for alphai and epsilonit
klh14ori <- plm(form, data = paneldata,
              model = "within", index = c("id","Year"), effect = "individual")
summary(klh14ori)
# alphai
alphai = as.numeric(unname(fixef(klh14ori)))
# epsilonit
epsiloni = as.numeric(klh14ori[["residuals"]])
summary(alphai)
summary(epsiloni)
# Constant for sfa estimation
constant = 1
# Persistent inefficiency
klh14per = as.data.frame(cbind(alphai,constant))
klh14p <- sfa(alphai~constant-1, data = klh14per, ineffDecrease = T, truncNorm = F, timeEffect = F)
summary(klh14p, effic = F, logDepVar = T, effMinusU = T)

  # Estimate individual persistent inefficiency
  inklh14p <- 1-efficiencies (klh14p, asInData = T, logDepVar = T, minusU = T)
  summary(inklh14p)
  
# Transitory inefficiency
klh14tran = as.data.frame(cbind(epsiloni,constant))
klh14t <- sfa(epsiloni~constant-1, data = klh14tran, ineffDecrease = T, truncNorm = F, timeEffect = T)
summary(klh14t, effic = F, logDepVar = T, effMinusU = T)

  # Estimate individual transitory inefficiency
  inklh14t <- 1-efficiencies (klh14t, asInData = T, logDepVar = T, minusU = T)
  summary(inklh14t)  

########################## BC95 ##########################
attach(paneldata)
# Model estimation
bc95 <- sfa(formz, data = paneldata, ineffDecrease = T, truncNorm = T, timeEffect = F)
summary(bc95, effic = F, logDepVar = T, effMinusU = T)
  
  # Estimate individual inefficiency
  inbc95 <- 1-efficiencies (bc95, asInData = T, logDepVar = T, minusU = T)
  summary(inbc95)
########################## SVKZ ##########################
library(np)

# Bandwidths selection
bws.r1 <- npregbw(forms, regtype="ll", data=data, bwmethod = "cv.ls", ckertype = "epanechnikov")
# Estimate conditional mean and extract fitted values and residuals
r1.est <- npreg(bws=bws.r1, gradients=TRUE)  
r1hat   <- fitted(r1.est)
e1hat  <- residuals(r1.est)
# Generate for r3
ehat3 <- e1hat^3

# Bandwidths selection for r3 (skewness measures)
bws.r3 <- npregbw(forms3, data=data, regtype="ll", bwmethod = "cv.ls", ckertype = "epanechnikov")
# Estimate conditional mean and extract fitted values and residuals
r3.est  <- npreg(bws=bws.r3, gradients=TRUE)
r3hat  <- fitted(r3.est)
e3hat  <- residuals(r3.est)

  # Estimate individual inefficiency
  sigu3.hat <- sqrt(pi/2)*(pi/(pi-4))*r3hat 
  sigu.hat <- apply(cbind(0,sigu3.hat),1,FUN=max)^(1/3)
  muhat <- sqrt(2/pi)*sigu.hat
  summary(muhat)
  insvkz = 1 - exp(-muhat) # Inefficiency estimation of individual unit
  summary(insvkz)

# Statistics of results ------------
  library(qpcR)
  # Combine inefficiency vectors of different length
  ineffs <- qpcR:::cbind.na(inals77, inss84, inpl81, inc90, inbc92, ing05, inklh14t, inklh14p,inbc95, insvkz)
  # Apply statistical analysis to each model
  stat = list("mean" = apply(ineffs, 2, mean, na.rm = T),
              "sd" = apply(ineffs, 2, sd, na.rm = T),
              "min" = apply(ineffs, 2, min, na.rm = T),
              "Q1" = apply(ineffs, 2, quantile, probs=0.25, na.rm = T),
              "Median" = apply(ineffs, 2, quantile, probs=0.5, na.rm = T),
              "Q3" = apply(ineffs, 2, quantile, probs=0.75, na.rm = T),
              "max" = apply(ineffs, 2, max, na.rm = T))
  
  write.csv(stat, file = 'Stats of inefficiency.csv')
  
  # Correlation between estimates and confidence interval of mean inefficiency
  cor(ineffs, method = c("spearman"))
  library(tidyverse)
  ci <- function(x) {
    t.test(x, conf.level = 0.95)$conf.int
  }
  apply(ineffs, 2, ci)

# Rank of individual unit ------------
  # Ranks of panel data setting
  Rankp = cbind(data$id, data$Year, 
                rank(inals77),
                as.integer(as.factor(rank(inpl81, ties.method = "average"))), #Tune for duplicate values
                as.integer(as.factor(rank(inc90, ties.method = "average"))),
                rank(inbc92),
                rank(inbc95),
                rank(ing05),
                rank(inklh14t),
                as.integer(as.factor(rank(insvkz, ties.method = "average"))))
  
  # Rename for Rankp
  Name <- c("id","Year","ALS77","PL81","CSS90","BC92","BC95","G05","KLH14-transitory","SVKZ")
  colnames(Rankp) <- Name
  
  # Ranks of cross-sectional data setting
  Rankc = cbind(unique(data$id), #Extract id for individual unit
                rank(inss84),
                rank(inklh14p))
  
  # Rename for Rankc
  Name <- c("id","SS84", "KLH14-persistent")
  colnames(Rankc) <- Name
  
  cor(Rankp, method = c("spearman"))
  cor(Rankc, method = c("spearman"))
  
# Kernel densities of estimations ------------
  require("ggplot2")
  .df <- na.omit(data.frame(x = inals77))
  .nbins <- pretty(range(.df$x), n = nclass.FD(.df$x), min.n = 1)
  .als<-ggplot(data = .df, aes(x = x, y = ..density..)) +
    # Epanechnikov kernel and CV bandwidth
    geom_density(
      kernel = "gaussian",
      bw = "ucv",
      alpha = 0.5,
      #Here for single group: color and fill without aes()
      color = "brown1", fill = "brown1",
      #shut the legend
      show.legend = FALSE
    ) +
    scale_y_continuous(expand = c(0.01, 0)) +
    xlab("ALS77") +
    ylab("Density") +
    
    RcmdrPlugin.KMggplot2::theme_simple(base_size = 14, base_family = "sans")
  print(.als)
  rm(.df, .nbins)
  
  require("ggplot2")
  .df <- na.omit(data.frame(x = inss84))
  .nbins <- pretty(range(.df$x), n = nclass.FD(.df$x), min.n = 1)
  .ss84<-ggplot(data = .df, aes(x = x, y = ..density..)) +
    # Epanechnikov kernel and CV bandwidth
    geom_density(
      kernel = "gaussian",
      bw = "ucv",
      alpha = 0.5,
      #Here for single group: color and fill without aes()
      color = "darkgreen", fill = "darkgreen",
      #shut the legend
      show.legend = FALSE
    ) +
    scale_y_continuous(expand = c(0.01, 0)) +
    xlab("SS84") +
    ylab("Density") +
    
    RcmdrPlugin.KMggplot2::theme_simple(base_size = 14, base_family = "sans")
  print(.ss84)
  rm(.df, .nbins)
  
  require("ggplot2")
  .df <- na.omit(data.frame(x = inpl81))
  .nbins <- pretty(range(.df$x), n = nclass.FD(.df$x), min.n = 1)
  .pl81<-ggplot(data = .df, aes(x = x, y = ..density..)) +
    # Epanechnikov kernel and CV bandwidth
    geom_density(
      kernel = "gaussian",
      bw = "ucv",
      alpha = 0.5,
      #Here for single group: color and fill without aes()
      color = "dodgerblue", fill = "dodgerblue",
      #shut the legend
      show.legend = FALSE
    ) +
    scale_y_continuous(expand = c(0.01, 0)) +
    xlab("PL81") +
    ylab("Density") +
    
    RcmdrPlugin.KMggplot2::theme_simple(base_size = 14, base_family = "sans")
  print(.pl81)
  rm(.df, .nbins)
  
  require("ggplot2")
  .df <- na.omit(data.frame(x = inc90))
  .nbins <- pretty(range(.df$x), n = nclass.FD(.df$x), min.n = 1)
  .c90<-ggplot(data = .df, aes(x = x, y = ..density..)) +
    # Epanechnikov kernel and CV bandwidth
    geom_density(
      kernel = "gaussian",
      bw = "ucv",
      alpha = 0.5,
      #Here for single group: color and fill without aes()
      color = "gold", fill = "gold",
      #shut the legend
      show.legend = FALSE
    ) +
    scale_y_continuous(expand = c(0.01, 0)) +
    xlab("C90") +
    ylab("Density") +
    
    RcmdrPlugin.KMggplot2::theme_simple(base_size = 14, base_family = "sans")
  print(.c90)
  rm(.df, .nbins)
  
  require("ggplot2")
  .df <- na.omit(data.frame(x = inbc92))
  .nbins <- pretty(range(.df$x), n = nclass.FD(.df$x), min.n = 1)
  .bc92<-ggplot(data = .df, aes(x = x, y = ..density..)) +
    # Epanechnikov kernel and CV bandwidth
    geom_density(
      kernel = "gaussian",
      bw = "ucv",
      alpha = 0.5,
      #Here for single group: color and fill without aes()
      color = "blueviolet", fill = "blueviolet",
      #shut the legend
      show.legend = FALSE
    ) +
    scale_y_continuous(expand = c(0.01, 0)) +
    xlab("BC92") +
    ylab("Density") +
    
    RcmdrPlugin.KMggplot2::theme_simple(base_size = 14, base_family = "sans")
  print(.bc92)
  rm(.df, .nbins)
  
  require("ggplot2")
  .df <- na.omit(data.frame(x = ing05))
  .nbins <- pretty(range(.df$x), n = nclass.FD(.df$x), min.n = 1)
  .g05<-ggplot(data = .df, aes(x = x, y = ..density..)) +
    # Epanechnikov kernel and CV bandwidth
    geom_density(
      kernel = "gaussian",
      bw = "ucv",
      alpha = 0.5,
      #Here for single group: color and fill without aes()
      color = "coral1", fill = "coral1",
      #shut the legend
      show.legend = FALSE
    ) +
    scale_y_continuous(expand = c(0.01, 0)) +
    xlab("G05") +
    ylab("Density") +
    
    RcmdrPlugin.KMggplot2::theme_simple(base_size = 14, base_family = "sans")
  print(.g05)
  rm(.df, .nbins)
  
  require("ggplot2")
  # Disable scientific notation
  options(scipen=1000)
  
  .df <- na.omit(data.frame(x = inklh14t))
  .nbins <- pretty(range(.df$x), n = nclass.FD(.df$x), min.n = 1)
  .klh14t<-ggplot(data = .df, aes(x = x, y = ..density..)) +
    # Epanechnikov kernel and CV bandwidth
    geom_density(
      kernel = "gaussian",
      bw = "ucv",
      alpha = 0.5,
      #Here for single group: color and fill without aes()
      color = "darkolivegreen1", fill = "darkolivegreen1",
      #shut the legend
      show.legend = FALSE
    ) +
    scale_y_continuous(expand = c(0.01, 0)) +
    xlab("KLH14-Transitory") +
    ylab("Density") +
    
    RcmdrPlugin.KMggplot2::theme_simple(base_size = 14, base_family = "sans")
  print(.klh14t)
  rm(.df, .nbins)
  
  require("ggplot2")
  .df <- na.omit(data.frame(x = inklh14p))
  .nbins <- pretty(range(.df$x), n = nclass.FD(.df$x), min.n = 1)
  .klh14p<-ggplot(data = .df, aes(x = x, y = ..density..)) +
    # Epanechnikov kernel and CV bandwidth
    geom_density(
      kernel = "gaussian",
      bw = "ucv",
      alpha = 0.5,
      #Here for single group: color and fill without aes()
      color = "royalblue", fill = "royalblue",
      #shut the legend
      show.legend = FALSE
    ) +
    scale_y_continuous(expand = c(0.01, 0)) +
    xlab("KLH14-Persistent") +
    ylab("Density") +
    
    RcmdrPlugin.KMggplot2::theme_simple(base_size = 14, base_family = "sans")
  print(.klh14p)
  rm(.df, .nbins)
  
  require("ggplot2")
  .df <- na.omit(data.frame(x = inbc95))
  .nbins <- pretty(range(.df$x), n = nclass.FD(.df$x), min.n = 1)
  .bc95<-ggplot(data = .df, aes(x = x, y = ..density..)) +
    # Epanechnikov kernel and CV bandwidth
    geom_density(
      kernel = "gaussian",
      bw = "ucv",
      alpha = 0.5,
      #Here for single group: color and fill without aes()
      color = "khaki3", fill = "khaki3",
      #shut the legend
      show.legend = FALSE
    ) +
    scale_y_continuous(expand = c(0.01, 0)) +
    xlab("BC95") +
    ylab("Density") +
    
    RcmdrPlugin.KMggplot2::theme_simple(base_size = 14, base_family = "sans")
  print(.bc95)
  rm(.df, .nbins)
  
  require("ggplot2")
  .df <- na.omit(data.frame(x = insvkz))
  .nbins <- pretty(range(.df$x), n = nclass.FD(.df$x), min.n = 1)
  .svkz<-ggplot(data = .df, aes(x = x, y = ..density..)) +
    # Epanechnikov kernel and CV bandwidth
    geom_density(
      kernel = "gaussian",
      bw = "ucv",
      alpha = 0.5,
      #Here for single group: color and fill without aes()
      color = "mediumorchid1", fill = "mediumorchid1",
      #shut the legend
      show.legend = FALSE
    ) +
    scale_y_continuous(expand = c(0.01, 0)) +
    xlab("SVKZ") +
    ylab("Density") +
    
    RcmdrPlugin.KMggplot2::theme_simple(base_size = 14, base_family = "sans")
  print(.svkz)
  rm(.df, .nbins)
  
  library("ggpubr")
  Inefficiencies <- ggarrange(.als, .ss84, .pl81, .c90, .bc92, .g05, .klh14t, .klh14p, .bc95, .svkz,
                              ncol = 2, nrow = 5)
  Inefficiencies
#----------------------END-----------------------------