# Data upload
# Variables information available in the relative document
dati = read.csv("FL_hazardCN.csv")
library(cluster)
library(CCA)
library(ggfortify)
View(dati)
str(dati)
head(dati)
pairs(dati)

# Linear transformation (Y=Ax+b) aiming to obtain two synthetic indicators:
# Fuel_hazard
# Propagation_hazard
A <- matrix(c(1/6, 1/3, 1/2, 0, 1, 0 ,0, 0,
              0, 0, 0, 1, 0, 0, 1.5, -1), byrow=T, nrow=2)
b <- matrix(c(-20, 0), ncol=1)
n <- dim(dati)[1]
Y <- as.matrix(dati) %*% t(A) + rep(1, n) %*% t(b)

# High hazard for observations with fuel hazard higher than 500:
length(which(Y[,1]>500))

# Summary statistics about indicators
colMeans(Y)
var(Y)
cor(Y)

# Principal Components Analysis (PCA)
# Performed on standardized data due to different scales of variables
pc_s <- prcomp(as.matrix(dati), scale. = TRUE)
summary(pc_s)
plot(pc_s, type="l", main="Screeplot") # see plot 1
cor(dati, pc_s$x)
# due to the criteria of the screeplot, of the value of the eigenvalues and by
# analyzing the correlation between the variables and the PCs, it's assumed that
# two PCs would be sufficient, however an appleasant amoutn of the explained variance
# (around 80%) is reached only by including four PCs
biplot(pc_s) # see plot 2
# Within the first CP all the variables, except RH and wind, act in the same direction as more
# the dry the soil, the higher the temperature, the faster the wind, the less humid the environment and the more it spreads
# quickly a fire. Hence it can be considered a measure of the speed of propagation of a
# fire. The second, third and fourth CP are more influenced by the variables respectively RH,
# wind and rain. By putting the information together, if the principle of parsimony is adopted, three can be kept
# components, if you prefer to have a larger explained variance you can keep four.
# From the graph it can be seen that there are two anomalous observations that differ from the others: number 380 and
# number 500. Then, we proceed with a new principal component analysis eliminating these observations:
pc_s2 <- prcomp(as.matrix(dati[-c(380,500),]), scale. = TRUE)
summary(pc_s2)
# Importance of components:
# PC1 PC2 PC3 PC4 PC5 PC6 PC7 PC8
# Standard deviation 1.6974 1.1728 1.1072 0.9300 0.85166 0.63754 0.5571
# Proportion of Variance 0.3602 0.1719 0.1532 0.1081 0.09067 0.05081 0.0388
# Cumulative Proportion 0.3602 0.5321 0.6853 0.7934 0.88409 0.93489 0.9737
# Standard deviation 0.45877
# Proportion of Variance 0.02631
# Cumulative Proportion 1.00000
plot(pc_s2, type="l", main="Screeplot") # see plot 3
biplot(pc_s2)  # see plot 4

# Even by looking at the previous correlations, there are no major changes compared to the previous case
# and the three-component solution could only be supported a little more. The results don't change anymore
# much because it must be kept in mind that the anomalous observations are two out of more than 500 statistical units,
# and they would have to have values far removed from those of the other units for them to affect the matrix of
# covariance enough to change the results of the principal component analysis.
# In the previous biplot, the observations relating to the days with the greatest risk or danger of fires are placed
# further to the right and lower than the others, i.e. in correspondence with high dryness and speed values 22
# fire propagation and low humidity and rainfall values.
# In the biplot in which the two anomalous observations were still present, these observations were always located at
# right but up (and not down towards higher values of RH and rain).

# Factor Analysis (EFA)
# checking for variables' normality
shapiro.test(dati$FFMC)
shapiro.test(dati$DMC)
shapiro.test(dati$DC)
# just by looking at the first three, the normality hypothesis is refused at all levels
# of statistical significance of alpha, therefore we prefer to conduct an EFA by using the
# method of the principal components instead of the MLE method which assumes normality for the variables.
eis = eigen(cor(dati))
loads = eis$vectors %*% diag(sqrt(eis$values))
eis$values
cumsum(eis$values/8) 

# reaching almost 80% of explained variance in the fourth factor.
# comunality and specifity estimation:
com4 = diag(loads[,1:4]%*%t(loads[,1:4]))
spec4 = 1 - com4
loads[,1:4]
# first factors appears highly described by the first 5 variables, exhibiting
# high negative correlation. [...]

# Canonical Correlations Analysis (CCA)
# Analyzing the correlations between the two group of variables:
# - Variables made by the Canadian fire prevention service
# - Atmospherical factors
X1 <- as.matrix(dati[,c(1:4)])
X2 <- as.matrix(dati[,-c(1:4)])
rho <- matcor(X1,X2)

# Heatmap
img.matcor(rho, type = 2) # see plot 5
cc.dati <- cc(X1,X2)
cc.dati$cor # 0.66790976 0.37712731 0.24674008 0.03043379 Canonical correlations

# Assessing ideal number of canonical dimensions
# TRV Test
hat.rho <- cc.dati$cor
p <- dim(X1)[2]
q <- dim(X2)[2]
n <- nrow(dati)
(w_b = -(n-1-0.5*(p+q+1))*sum(log(1-hat.rhoˆ2)))
w_b> qchisq(p = 0.95, df = p*q)
(w_b2 = -(n-1-0.5*(p+q+1))*sum(log(1-hat.rhoˆ2)[-1]))
w_b2 > qchisq(p = 0.95, df = (p-1)*(q-1))
(w_b3 = -(n-1-0.5*(p+q+1))*sum(log(1-hat.rhoˆ2)[-(1:2)]))
w_b3 > qchisq(p = 0.95, df = (p-2)*(q-2))
(w_b4 = -(n-1-0.5*(p+q+1))*sum(log(1-hat.rhoˆ2)[-(1:3)]))
w_b4 > qchisq(p = 0.95, df = (p-3)*(q-3))
# keeping three canonical dimension to describe the correlation 
# phenomena between the two groups

# Interpretation and canonical correlations plot
cc.dati$scores$corr.X.xscores
cc.dati$scores$corr.Y.yscores
plt.cc(cc.dati, var.label = TRUE, type="b") # see plot 6
# First of all, it is noted that uˆ1 is strongly negatively correlated with DMC and DC (but also with ISI e                                                                       FFMS) 
# while vˆ1 is strongly negatively correlated only with temp. Since uˆ1 and vˆ1 are positive
# correlated, the following reasoning is obtained: if the dryness of the fuel increases (just as it increases
# the index of how quickly a fire spreads) then uˆ1 decreases, but if uˆ1 decreases then also
# vˆ1 decreases, and if vˆ1 decreases then the temperature increases. In summary, a high (resp. low)
# temperature increases (resp. decreases) both the dryness of the fuel of any layer and the
# rate of spread of a fire.
# As regards the second canonical dimension, it is noted that uˆ2 is strongly positively correlated
# with FFMC (but also with ISI) while vˆ2 is negatively correlated with RH. Since again uˆ2 and vˆ2 are
# positively correlated, with the same reasoning as before we obtain that a high (resp. low) humidity
# leads to a decrease (resp. increase) of the dryness of the fuel, and thus to a lower (resp.
# greater) rapidity of spread of a fire.

# HIERARCHICAL AGGLOMERATIVE CLUSTERING
# Assessing euclidean distance with Ward's linkage
# Working on the standardized dataset due to the different scales of variables
dati_st = scale(dati, scale = TRUE, center = TRUE)
dm <- dist(dati_st, method="euclidean")
hcw <- hclust(dm, method="ward.D2")

# dendrogram
plot(hcw)  # see plot 7

# Looking at the dendrogram based on standardized data, there would appear to be two groups: that
# small one on the left and the larger one on the right, obtained for example by cutting the dendogram at height 3 40.
# All other distances seem relatively small, so there would appear to be no possibility of
# detect more than two groups.
# Now let's consider a division into 7 groups and the relative scatterplots:
pairs(dati, col=(cutree(hcw, 7))) # see plot 8

# DMC vs DC
pairs(dati[,2:3], col=(cutree(hcw, 7)), pch = 19) # see plot 9

# temperature vs wind speed
pairs(dati[,5:6], col=(cutree(hcw, 7)), pch = 19) # see plot 10

# representation on the first two PCs dimensions
clusplot(dati_st, cutree(hcw, k = 7)) # see plot 11

# NON-HIERARCHICAL CLUSTERING EVALUTAION
# assessing the number of cluster to consider:
km1 = kmeans(dati_st, centers = 1)
km2 = kmeans(dati_st, centers = 2)
km3 = kmeans(dati_st, centers = 3)
km4 = kmeans(dati_st, centers = 4)
km5 = kmeans(dati_st, centers = 5)
km6 = kmeans(dati_st, centers = 6)
explained.var = c(km1$betweenss/km1$totss,
                  km2$betweenss/km2$totss,
                  km3$betweenss/km3$totss,
                  km4$betweenss/km4$totss,
                  km5$betweenss/km5$totss,
                  km6$betweens/km6$totss)
plot(c(1:6), -explained.var, type = "l", lwd = 2) # see plot 12
# Assuming 2 or 4 clusters, further evaluations to decide:

# Representation of the cluster (2 and 4) on the first two PCs dimension:
autoplot(km2, data = dati_st, frame = TRUE) # see plot 13
autoplot(km4, data = dati_st, frame = TRUE) # see plot 14
# two groups appear to be better identifying with less overlapping

# Silhouette Analysis
# Silhouette analysis for two-cluster identification assessment:
d = dist(dati)
plot(silhouette(km2$cluster, d)) # see plot 15