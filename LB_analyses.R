#____________________________________________________________________________________________
#
# Entomological Society of America Conference
#
# Statistics Workshop
#
# Using museum records to infer changes in insect communities over time
#
# Changes in lady beetle communities over time in Ohio using museum records
#
# Perry, KI, CA Bahlai, TJ Assal, CB Riley, KJ Turo, L Taylor, J Radl,
# YA Delgado de la flor, FS Sivakoff, and MM Gardiner. 2024. Landscape change
# and alien invasions drive shifts in native lady beetle communities over a
# century, Ecological Applications, 34(7): e3024.
# https://doi-org.proxy.lib.ohio-state.edu/10.1002/eap.3024
#
# Bahlai, C and KI Perry. 2024. BahlaiLab/Ohio_ladybeetles: Ecological
# applications manuscript final code (V1.1). Zenodo. https://doi.org/10.5281/zenodo.11263088 
#
# Kayla I Perry & Christie A Bahlai
# 9 November 2025
#____________________________________________________________________________________________

# Part 1: Partition beta-diversity into turnover and nestedness components ----

# load in the data
lb.b <- read.csv("data/ladybeetle_betadiversity.csv")
str(lb.b)

# install.packages("reshape2")
library(reshape2)

# install.packages("betapart")
library(betapart)

# change data from long format to wide format
lb.b.matrix <- dcast(lb.b, Decade ~ Name, length)
str(lb.b.matrix)

# remove 1890 due to few collections
# restricting analyses to 1900-2018
lb.b.matrix <- lb.b.matrix[-1,]
str(lb.b.matrix)

# remove single collection of Mulsantina luteodorsa because it hasn't been
# verified by LB experts
lb.b.matrix <- lb.b.matrix[,-23]
str(lb.b.matrix)
as.data.frame(lb.b.matrix)

# change decade variable to factor
lb.b.matrix$Decade <- as.factor(lb.b.matrix$Decade)
levels(lb.b.matrix$Decade)

# change data set to presence/absence (1/0)
# to calculate occurrence based beta-diversity metric
lb.b.matrix[lb.b.matrix>0]<-1
str(lb.b.matrix)

# add a column for site (i.e. Ohio) so can make row names later
# will look at changes in capture over time across the entire state
lb.b.matrix$site <- c("OH", "OH", "OH", "OH", "OH", "OH", "OH", "OH", "OH", "OH", "OH", "OH")
str(lb.b.matrix)

# Descriptive patterns for changes in total lady beetle beta diversity
# across decades in Ohio

# pull out data from each decade
# remove factor columns, change site to row names

time1 <- lb.b.matrix[which(lb.b.matrix$Decade == "1900"),]
rownames(time1) <- time1[,29]
time1 <- time1[,-1]
time1 <- time1[,-28]
str(time1)

time2 <- lb.b.matrix[which(lb.b.matrix$Decade == "1910"),]
rownames(time2) <- time2[,29]
time2 <- time2[,-1]
time2 <- time2[,-28]
str(time2)

time3 <- lb.b.matrix[which(lb.b.matrix$Decade == "1920"),]
rownames(time3) <- time3[,29]
time3 <- time3[,-1]
time3 <- time3[,-28]
str(time3)

time4 <- lb.b.matrix[which(lb.b.matrix$Decade == "1930"),]
rownames(time4) <- time4[,29]
time4 <- time4[,-1]
time4 <- time4[,-28]
str(time4)

time5 <- lb.b.matrix[which(lb.b.matrix$Decade == "1940"),]
rownames(time5) <- time5[,29]
time5 <- time5[,-1]
time5 <- time5[,-28]
str(time5)

time6 <- lb.b.matrix[which(lb.b.matrix$Decade == "1950"),]
rownames(time6) <- time6[,29]
time6 <- time6[,-1]
time6 <- time6[,-28]
str(time6)

time7 <- lb.b.matrix[which(lb.b.matrix$Decade == "1960"),]
rownames(time7) <- time7[,29]
time7 <- time7[,-1]
time7 <- time7[,-28]
str(time7)

time8 <- lb.b.matrix[which(lb.b.matrix$Decade == "1970"),]
rownames(time8) <- time8[,29]
time8 <- time8[,-1]
time8 <- time8[,-28]
str(time8)

time9 <- lb.b.matrix[which(lb.b.matrix$Decade == "1980"),]
rownames(time9) <- time9[,29]
time9 <- time9[,-1]
time9 <- time9[,-28]
str(time9)

time10 <- lb.b.matrix[which(lb.b.matrix$Decade == "1990"),]
rownames(time10) <- time10[,29]
time10 <- time10[,-1]
time10 <- time10[,-28]
str(time10)

time11 <- lb.b.matrix[which(lb.b.matrix$Decade == "2000"),]
rownames(time11) <- time11[,29]
time11 <- time11[,-1]
time11 <- time11[,-28]
str(time11)

time12 <- lb.b.matrix[which(lb.b.matrix$Decade == "2010"),]
rownames(time12) <- time12[,29]
time12 <- time12[,-1]
time12 <- time12[,-28]
str(time12)

# create beta part object for each decade
# betapart.core function computes the basic quantities needed for computing the
# multiple-site beta diversity measures and pairwise dissimilarity matrices

t1 <- betapart.core(time1)
t2 <- betapart.core(time2)
t3 <- betapart.core(time3)
t4 <- betapart.core(time4)
t5 <- betapart.core(time5)
t6 <- betapart.core(time6)
t7 <- betapart.core(time7)
t8 <- betapart.core(time8)
t9 <- betapart.core(time9)
t10 <- betapart.core(time10)
t11 <- betapart.core(time11)
t12 <- betapart.core(time12)

# assess temporal change in community composition among the decades for all of Ohio
# beta.temp function computes the dissimilarity for each locality between time 1 and time 2,
# including the total beta-diversity and the turnover and nestedness components

t1.2 <- beta.temp(t1, t2, index.family = "sorensen")
t1.2

t2.3 <- beta.temp(t2, t3, index.family = "sorensen")
t2.3

t3.4 <- beta.temp(t3, t4, index.family = "sorensen")
t3.4

t4.5 <- beta.temp(t4, t5, index.family = "sorensen")
t4.5

t5.6 <- beta.temp(t5, t6, index.family = "sorensen")
t5.6

t6.7 <- beta.temp(t6, t7, index.family = "sorensen")
t6.7

t7.8 <- beta.temp(t7, t8, index.family = "sorensen")
t7.8

t8.9 <- beta.temp(t8, t9, index.family = "sorensen")
t8.9

t9.10 <- beta.temp(t9, t10, index.family = "sorensen")
t9.10

t10.11 <- beta.temp(t10, t11, index.family = "sorensen")
t10.11

t11.12 <- beta.temp(t11, t12, index.family = "sorensen")
t11.12

lb.beta <- rbind(t1.2, t2.3, t3.4, t4.5, t5.6, t6.7, t7.8, t8.9, t9.10, t10.11, t11.12)
str(lb.beta)

# add decade back in to the data set

lb.beta$decade <- c("1900-1910", "1910-1920", "1920-1930", "1930-1940", "1940-1950", "1950-1960", "1960-1970", "1970-1980", "1980-1990", "1990-2000", "2000-2010")
str(lb.beta)
lb.beta$decade <- as.factor(lb.beta$decade)

# let's graph the results
# first need to change from wide format to long format
# and remove total beta-diversity

lb.beta.long <- melt(lb.beta, id.vars = c("decade"))
colnames(lb.beta.long) <- c("decade", "beta.var", "beta.div")

beta.sim <- lb.beta.long[which(lb.beta.long$beta.var == "beta.sim"),]
beta.sne <- lb.beta.long[which(lb.beta.long$beta.var == "beta.sne"),]

beta.final <- rbind(beta.sim, beta.sne)
beta.final <- droplevels(beta.final)
str(beta.final)

library(ggplot2)

ggplot(data = beta.final, aes(x = decade, y = beta.div, fill = beta.var)) +
  geom_bar(stat = "identity", position = "stack") +
  coord_flip() +
  theme_classic() +
  ylab("Total Beta-diversity (Bsor)") +
  xlab("Decade") +
  theme(legend.title = element_blank()) +
  scale_fill_manual(values = c("gray", "black"), labels = c("Turnover", "Nestedness"))


#________________________________________________________________________________________
# Part 2: Nonmetric multidimensional scaling ----

# load in the data
lb_raw <- read.csv("data/ladybeetle_data.csv")
str(lb_raw)

LULC <- read.csv("data/landscape_data.csv")
str(LULC)

LULC.raw.wide <- dcast(LULC, Year + County ~ Class, value.var = "Percentage")

lb_raw$Name <- paste(lb_raw$Genus, lb_raw$Species, sep = ".")
lb.wide.yeargroups <- dcast(lb_raw, YearGroup + County ~ Name, length)

# NMDS analysis will not run if any sites (i.e., rows) have 0 occurrences across all species
# So, need to remove county/year combinations with insufficient observations

cutpoint <- rowSums(lb.wide.yeargroups[3:30])
lb.wide.yeargroups <- lb.wide.yeargroups[which(cutpoint>2), ]

# create a dataset with just native lady beetle species (i.e., remove exotic species)
lb.wide.yeargroups.natives <- lb.wide.yeargroups[,-c(11,14,17,23,29)]

# Same issue, need to remove county/year combinations with insufficient observations
cutpoint1 <- rowSums(lb.wide.yeargroups.natives[3:25])
lb.wide.yeargroups.natives <- lb.wide.yeargroups.natives[which(cutpoint1>2), ]

# install.packages("vegan")
library(vegan)

# landscape
landscape.ord<-metaMDS(LULC.raw.wide[3:6])
landscape.ord # stress is quality of fit

plot(landscape.ord, disp='sites', type="n")
points(landscape.ord, display="sites", select=which(LULC.raw.wide$Year=="1938"), pch=19, col="pink")
points(landscape.ord, display="sites", select=which(LULC.raw.wide$Year=="1970"), pch=19, col="orange")
points(landscape.ord, display="sites", select=which(LULC.raw.wide$Year=="1992"), pch=19, col="green")
points(landscape.ord, display="sites", select=which(LULC.raw.wide$Year=="2016"), pch=19, col="blue")
text(landscape.ord, display="species", col="red", pch=1)
ordiellipse(landscape.ord, LULC.raw.wide$Year, draw="polygon", col=c("pink", "orange", "green", "blue"), kind="se", conf=0.999, label=TRUE)

# all lady beetles
lb.ord<-metaMDS(lb.wide.yeargroups[3:30])
lb.ord

plot(lb.ord, disp='sites', type="n")
points(lb.ord, display="sites", select=which(lb.wide.yeargroups$YearGroup=="1938"), pch=19, col="pink")
points(lb.ord, display="sites", select=which(lb.wide.yeargroups$YearGroup=="1970"), pch=19, col="orange")
points(lb.ord, display="sites", select=which(lb.wide.yeargroups$YearGroup=="1992"), pch=19, col="green")
points(lb.ord, display="sites", select=which(lb.wide.yeargroups$YearGroup=="2016"), pch=19, col="blue")
ordiellipse(lb.ord, lb.wide.yeargroups$YearGroup, draw="polygon", col=c("pink", "orange", "green", "blue"), kind="se", conf=0.999, label=TRUE)

# native lady beetles
lb.ord.n<-metaMDS(lb.wide.yeargroups.natives[3:25])
lb.ord.n

plot(lb.ord.n, disp='sites', type="n")
points(lb.ord.n, display="sites", select=which(lb.wide.yeargroups.natives$YearGroup=="1938"), pch=19, col="pink")
points(lb.ord.n, display="sites", select=which(lb.wide.yeargroups.natives$YearGroup=="1970"), pch=19, col="orange")
points(lb.ord.n, display="sites", select=which(lb.wide.yeargroups.natives$YearGroup=="1992"), pch=19, col="green")
points(lb.ord.n, display="sites", select=which(lb.wide.yeargroups.natives$YearGroup=="2016"), pch=19, col="blue")
ordiellipse(lb.ord.n, lb.wide.yeargroups.natives$YearGroup, draw="polygon", col=c("pink", "orange", "green", "blue"), kind="se", conf=0.999, label=TRUE)


library(ggordiplots)

# all lady beetles
# another approach to graph the data

lbnmds <- gg_ordiplot(lb.ord, groups = lb.wide.yeargroups$YearGroup, kind = "se", conf = 0.99, pt.size = -1, plot = F)

gglbnmds <- lbnmds$plot +  
  geom_point(data=lbnmds$df_ord, aes(x=x,y=y, color=Group, shape=Group))+
  scale_colour_manual(values=c("1938" = "darkred", "1970" = "darkorange", "1992"= "yellow3", "2016"="darkgreen"), 
                      labels=c("Before 1938", "1939-1970", "1971-1992", "1993-2016"), title("Year"))+
  scale_fill_manual(values=c("1938" = "darkred", "1970" = "darkorange", "1992"= "yellow3", "2016"="darkgreen"), 
                    labels=c("1938", "1970", "1992", "2016"), title("Year"))+
  scale_shape_manual(values=c("1938" = 1, "1970" = 2, "1992"= 3, "2016"=4), 
                     labels=c("Before 1938", "1939-1970", "1971-1992", "1993-2016"), title("Year"))+
  geom_polygon(data = lbnmds$df_ellipse, aes(x = x, y = y,  fill=Group), show.legend = FALSE, color="black", alpha =0.25)+
  geom_path(data = lbnmds$df_mean.ord, aes(x = x, y = y), 
            show.legend = FALSE, color="black")+
  geom_point(data = lbnmds$df_mean.ord, aes(x = x, y = y), 
             show.legend = FALSE, color="black", pch=16)+
  theme_classic()+
  theme(aspect.ratio = 0.9)+
  labs(x="NMDS1", y="NMDS2")

gglbnmds

#________________________________________________________________________________________
# Part 3: Responses of native lady beetles using generalized additive models (GAMs) ----

library(mgcv) 
library(visreg)
library(ggplot2)


lb_all <- read.csv("data/ladybeetle_gams.csv")
str(lb_all)

# Evaluate relative importance of temporal, spatial, landscape, and community factors on
# the captures of native lady beetle species using negative binomial generalized additive models (GAMs)
# Use Hippodamia covergens as an example

# first, how many captures are we working with?
sum(lb_all$Hippodamia.convergens)


# iterative process for model selection - use AIC as selection criterion
# the Akaike Information Criterion (AIC) is a comparative metric of model performance,
# where lower scores indicate improved performance/fit compared to other considered models.

# GAMs capture non-linear relationships by fitting a smooth function through the data
# smooth terms are specified using s(x), where x is the predictor variable
# multiple smooth terms to select: https://stat.ethz.ch/R-manual/R-devel/library/mgcv/html/smooth.terms.html
# s() indicates the smooth function, of which the default is thin ploate regression spline
# The smooth functions have several parameters that can be set to change their behavior.


hcon.gam1 <- gam(Hippodamia.convergens ~ offset(log(1+Totalcount-Hippodamia.convergens)) +
                 s(Decade, sp=0.5, k=4)+
                 s(lon, sp=0.5)+
                 s(lat, sp=0.5)+
                 s(log(1+Totalinvasive), sp=0.5, k=4)+
                 s(Agriculture, sp=0.5)+
                 s(Forest, sp=0.5)+
                 s(Developed, sp=0.5), 
               data=lb_all, family="nb")
summary(hcon.gam1)
# Effective degrees of freedon (edf) is related to the smoothing parameter, with higher edf
# indicating the term is more non-linear
# k sets the number of base functions used to create a smooth function
# higher k = more smoothing, but be careful to not overfit!

AIC(hcon.gam1)

# visreg produces a plot showing the predicted relationship between a predictor and the response,
# while holding other variables constant
# plot includes fitted curve, confidience band, and partial residuals

visreg(hcon.gam1, "Decade",  ylab="Captures")
visreg(hcon.gam1, "lon",  ylab="Captures")
visreg(hcon.gam1, "lat",  ylab="Captures")
visreg(hcon.gam1, "Totalinvasive",  ylab="Captures")
visreg(hcon.gam1, "Agriculture",  ylab="Captures")
visreg(hcon.gam1, "Forest", ylab="Captures")
visreg(hcon.gam1, "Developed",  ylab="Captures")


# replace total invasives with proportion of invasives

hcon.gam2 <- gam(Hippodamia.convergens ~ offset(log(1+Totalcount-Hippodamia.convergens)) +
                   s(Decade, sp=0.5, k=4)+
                   s(lon, sp=0.5)+
                   s(lat, sp=0.5)+
                   s(Propinvasive, sp=0.5)+
                   s(Agriculture, sp=0.5)+
                   s(Forest, sp=0.5)+
                   s(Developed, sp=0.5),  
                 data=lb_all, family="nb")
summary(hcon.gam2)
AIC(hcon.gam2)

visreg(hcon.gam2, "Decade",  ylab="Captures")
visreg(hcon.gam2, "lon",  ylab="Captures")
visreg(hcon.gam2, "lat",  ylab="Captures")
visreg(hcon.gam2, "Propinvasive",  ylab="Captures")
visreg(hcon.gam2, "Agriculture",  ylab="Captures")
visreg(hcon.gam2, "Forest", ylab="Captures")
visreg(hcon.gam2, "Developed",  ylab="Captures")

# replace proportion invasive with the two major invasive species

hcon.gam3 <- gam(Hippodamia.convergens ~ offset(log(1+Totalcount-Hippodamia.convergens)) +
                 s(Decade, sp=0.5, k=4)+
                 s(lon, sp=0.5)+
                 s(lat, sp=0.5)+
                 s(log(1+Coccinella.septempunctata), sp=0.5, k=4)+
                 s(log(1+Harmonia.axyridis), sp=0.5, k=4)+
                 s(Agriculture, sp=0.5)+
                 s(Forest, sp=0.5)+
                 s(Developed, sp=0.5),  
               data=lb_all, family="nb")
summary(hcon.gam3)
AIC(hcon.gam3)

visreg(hcon.gam3, "Decade",  ylab="Captures")
visreg(hcon.gam3, "lon",  ylab="Captures")
visreg(hcon.gam3, "lat",  ylab="Captures")
visreg(hcon.gam3, "Coccinella.septempunctata",  ylab="Captures")
visreg(hcon.gam3, "Harmonia.axyridis",  ylab="Captures")
visreg(hcon.gam3, "Agriculture",  ylab="Captures")
visreg(hcon.gam3, "Forest", ylab="Captures")
visreg(hcon.gam3, "Developed",  ylab="Captures")


# Used AIC model selection to assess landscape parameters (not shown)

# final model: 

hcon.gam4 <- gam(Hippodamia.convergens ~ offset(log(1+Totalcount-Hippodamia.convergens)) +
                 s(Decade, sp=0.5, k=4)+
                 s(lat, sp=0.5)+
                 s(log(1+Coccinella.septempunctata), sp=0.5, k=4)+
                 s(log(1+Harmonia.axyridis), sp=0.5, k=4)+
                 s(Developed, sp=0.5),  
               data=lb_all, family="nb")
summary(hcon.gam4)
AIC(hcon.gam4)

# for all the 'final' models as chosen by AIC, check the concurvity to identify variables where concurvity is 
# likely to be causing an issue in the final fit. Eliminate variables with overall estimated concurvity >0.8 
# Concurvity occurs when some smooth term in a model could be approximated by one or more of the other smooth terms in the model.
# Causes similar problems as colinearity

concurvity(hcon.gam4)
#https://stat.ethz.ch/R-manual/R-devel/library/mgcv/html/concurvity.html

# Decade
hcon.decade <- visreg(hcon.gam3, "Decade",  ylab="Residual captures",
                    xlab=expression(paste("Year")),
                    gg=T,
                    line=list(col="gray19"),
                    fill=list(col="gray57", fill="gray57"),
                    points=list(size=1, pch=21, fill="gray19", col="black"))+
  theme_classic()
hcon.decade

# Latitude
hcon.lat<-visreg(hcon.gam3, "lat",  ylab="Residual captures",
                 xlab=expression(paste("Latitude")),
                 gg=T,
                 line=list(col="gray28"),
                 fill=list(col="gray", fill="gray"),
                 points=list(size=1, pch=21, fill="gray28", col="black"))+
  theme_classic()
hcon.lat

# C. septempunctata
hcon.c7<-visreg(hcon.gam3, "Coccinella.septempunctata",  ylab="Residual captures",
                xlab=expression(paste("Captures of ", italic("Coccinella septempunctata"))),
                gg=T,
                line=list(col="darkred"),
                fill=list(col="lightcoral", fill="lightcoral"),
                points=list(size=1, pch=22, fill="darkred", col="black"))+
  theme_classic()
hcon.c7

# H. axyridis
hcon.ha<-visreg(hcon.gam3, "Harmonia.axyridis",  ylab="Residual captures",
                xlab=expression(paste("Captures of ", italic("Harmonia axyridis"))),
                gg=T,
                line=list(col="darkorange"),
                fill=list(col="burlywood1", fill="burlywood1"),
                points=list(size=1, pch=22, fill="darkorange", col="black"))+
  theme_classic()
hcon.ha

# Developed
hcon.developed<-visreg(hcon.gam3, "Developed", ylab="Residual Captures", xlab="% Developed cover",
                       gg=T,
                       line=list(col="slategray4"),
                       fill=list(col="slategray2", fill="slategray2"),
                       points=list(size=1, pch=23, fill="slategray4", col="black"))+
  theme_classic()
hcon.developed
