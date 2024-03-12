rm(list=ls())
graphics.off()
# Generate synthetic panel data set ------------------------------------

set.seed(99)  # Set fixed seed for synthetic data

# Unit ID and Year:
id = sort(rep(1:100, times = 4))
Year = rep(1:4, times = 100)
t = Year
t2 = t^2

# Output
Y = rnorm(400, mean=-2, sd=1.6)

# Inputs
X1 = rnorm(400, mean=3.3, sd=1.3)
X2 = rnorm(400, mean=-1.6, sd=1.5)
X3 = rnorm(400, mean=13.6, sd=2)

# Environmental
#Assume 20% of units are 1 in Z1
Z1 = rep(c(1, 0), times= c(20,80))[order(runif(100))]
Z1 = rep(Z1, each = 4) # Replicate for a panel data set
Z2 = rep(c(1, 0), times= c(30,70))[order(runif(100))]
Z2 = rep(Z2, each = 4)
Z3 = rep(c(1, 0), times= c(80,20))[order(runif(100))]
Z3 = rep(Z3, each = 4)

# Composed data set
data = as.data.frame(cbind(id, Year, Y, X1, X2, X3, Z1, Z2, Z3, t, t2))
# Convert to panel data
library(plm)
paneldata<- pdata.frame(data, c("id","Year"))

# Define formulas for different models
form = Y ~ X1 + X2 + X3 # For ALS77/SS84/PL81/BC92/KLH14
formz = Y ~ X1 + X2 + X3 | Z1 + Z2 + Z3 # For BC95
formt = Y~ X1 + X2 + X3 + factor(id) # For G05
formc = Y~ X1 + X2 + X3 + t + t2 # For CSS90

fZ1 = as.factor(Z1)
fZ2 = as.factor(Z2)
fZ3 = as.factor(Z3)

forms = Y ~ X1 + X2 + X3 + fZ1 + fZ2 + fZ3 # For SVKZ
forms3 = ehat3 ~ X1 + X2 + X3 + fZ1 + fZ2 + fZ3 #For SVKZ

# The synthetic data set can be then applied with the SFM code as illustrated in Appendix B.
