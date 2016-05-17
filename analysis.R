library('vtreat')
library('ggplot2')
library('tidyr')
library('devtools')

devtools::install_github('WinVector/WVPlots',build_vignettes = TRUE)
library('WVPlots')

setwd('C:/Users/EvanMi/Desktop/Stats/PCA_part1')
source('helper_functions.R')

# make data
set.seed(23525)
dTrain <- mkData(1000)
dTest <- mkData(1000)

summary(dTrain[, c("y", "x.01", "x.02",
                   "noise1.01", "noise1.02")])

goodVars <-  colnames(dTrain)[grep('^x.',colnames(dTrain))]
dTrainIdeal <- dTrain[,c('y',goodVars)]
dTestIdeal <-  dTrain[,c('y',goodVars)]

# do the PCA
dmTrainIdeal <- as.matrix(dTrainIdeal[,goodVars])
princIdeal <- prcomp(dmTrainIdeal,center = TRUE,scale. = TRUE)

# extract the principal components
rot5Ideal <- extractProjection(5,princIdeal)

# prepare the data to plot the variable loadings
rotfIdeal = as.data.frame(rot5Ideal)
rotfIdeal$varName = rownames(rotfIdeal)
rotflongIdeal = gather(rotfIdeal, "PC", "loading",
                       starts_with("PC"))
rotflongIdeal$vartype = ifelse(grepl("noise", 
                                     rotflongIdeal$varName),
                               "noise", "signal")

# plot the singular values
dotplot_identity(frame = data.frame(pc=1:length(princIdeal$sdev), 
                                    magnitude=princIdeal$sdev), 
                 xvar="pc",yvar="magnitude") +
  ggtitle("Ideal case: Magnitudes of singular values")

dotplot_identity(rotflongIdeal, "varName", "loading", "vartype") + 
  facet_wrap(~PC,nrow=1) + coord_flip() + 
  ggtitle("x scaled variable loadings, first 5 principal components") + 
  scale_color_manual(values = c("noise" = "#d95f02", "signal" = "#1b9e77"))

# signs are arbitrary on PCA, so instead of calling predict we pull out
# (and alter) the projection by hand
projectedTrainIdeal <-
  as.data.frame(dmTrainIdeal %*% extractProjection(2,princIdeal),
                stringsAsFactors = FALSE)
projectedTrainIdeal$y <- dTrain$y
ScatterHistN(projectedTrainIdeal,'PC1','PC2','y',
             "Ideal Data projected to first two principal components")

# Without scaling

vars <- setdiff(colnames(dTrain),'y')

duTrain <- as.matrix(dTrain[,vars])
prinU <- prcomp(duTrain,center = TRUE,scale. = FALSE) 

dotplot_identity(frame = data.frame(pc=1:length(prinU$sdev), 
                                    magnitude=prinU$sdev), 
                 xvar="pc",yvar="magnitude") +
  ggtitle("Unscaled case: Magnitudes of singular values")

rot5U <- extractProjection(5,prinU)
rot5U = as.data.frame(rot5U)
rot5U$varName = rownames(rot5U)
rot5U = gather(rot5U, "PC", "loading",
               starts_with("PC"))
rot5U$vartype = ifelse(grepl("noise", 
                             rot5U$varName),
                       "noise", "signal")

dotplot_identity(rot5U, "varName", "loading", "vartype") + 
  facet_wrap(~PC,nrow=1) + coord_flip() + 
  ggtitle("unscaled variable loadings, first 5 principal components") + 
  scale_color_manual(values = c("noise" = "#d95f02", "signal" = "#1b9e77"))

# get all the principal components
# not really a projection as we took all components!
projectedTrain <- as.data.frame(predict(prinU,duTrain),
                                stringsAsFactors = FALSE)
vars = colnames(projectedTrain)
projectedTrain$y <- dTrain$y

varexpr = paste(vars, collapse="+")
fmla = paste("y ~", varexpr)

model <- lm(fmla,data=projectedTrain)
summary(model)

estimate <- predict(model,newdata=projectedTrain)
trainrsq <- rsq(estimate,projectedTrain$y)

# Scaling X

dTrainNTreatedUnscaled <- dTrain
dTestNTreatedUnscaled <- dTest

# scale the data
dTrainNTreatedXscaled <- 
  as.data.frame(scale(dTrainNTreatedUnscaled[,colnames(dTrainNTreatedUnscaled)!='y'],
                      center=TRUE,scale=TRUE),stringsAsFactors = FALSE)
dTrainNTreatedXscaled$y <- dTrainNTreatedUnscaled$y
dTestNTreatedXscaled <- 
  as.data.frame(scale(dTestNTreatedUnscaled[,colnames(dTestNTreatedUnscaled)!='y'],
                      center=TRUE,scale=TRUE),stringsAsFactors = FALSE)
dTestNTreatedXscaled$y <- dTestNTreatedUnscaled$y

# get the variable ranges
ranges = vapply(dTrainNTreatedXscaled, FUN=function(col) c(min(col), max(col)), numeric(2))
rownames(ranges) = c("vmin", "vmax") 
rframe = as.data.frame(t(ranges))  # make ymin/ymax the columns
rframe$varName = rownames(rframe)
varnames = setdiff(rownames(rframe), "y")
rframe = rframe[varnames,]
rframe$vartype = ifelse(grepl("noise", rframe$varName),
                        "noise", "signal")

summary(dTrainNTreatedXscaled[, c("y", "x.01", "x.02", 
                                  "noise1.01", "noise1.02")])

barbell_plot(rframe, "varName", "vmin", "vmax", "vartype") +
  coord_flip() + ggtitle("x scaled variables: ranges") + 
  scale_color_manual(values = c("noise" = "#d95f02", "signal" = "#1b9e77"))

vars = setdiff(colnames(dTrainNTreatedXscaled), "y")

dmTrain <- as.matrix(dTrainNTreatedXscaled[,vars])
dmTest <- as.matrix(dTestNTreatedXscaled[,vars])
princ <- prcomp(dmTrain,center = TRUE,scale. = TRUE) 
dotplot_identity(frame = data.frame(pc=1:length(princ$sdev), 
                                    magnitude=princ$sdev), 
                 xvar="pc",yvar="magnitude") +
  ggtitle("x scaled variables: Magnitudes of singular values")

# Just looking at the first five PCs

rot5 <- extractProjection(5,princ)
rotf = as.data.frame(rot5)
rotf$varName = rownames(rotf)
rotflong = gather(rotf, "PC", "loading", starts_with("PC"))
rotflong$vartype = ifelse(grepl("noise", rotflong$varName), 
                          "noise", "signal")

dotplot_identity(rotflong, "varName", "loading", "vartype") + 
  facet_wrap(~PC,nrow=1) + coord_flip() + 
  ggtitle("x scaled variable loadings, first 5 principal components") + 
  scale_color_manual(values = c("noise" = "#d95f02", "signal" = "#1b9e77"))

# Using the first 20 PCs to model

# get all the principal components
# not really a projection as we took all components!
projectedTrain <- as.data.frame(predict(princ,dmTrain),
                                stringsAsFactors = FALSE)
projectedTrain$y <- dTrainNTreatedXscaled$y

ncomp = 20
# here we will only model with the first ncomp principal components
varexpr = paste(paste("PC", 1:ncomp, sep=''), collapse='+')
fmla = paste("y ~", varexpr)

model <- lm(fmla,data=projectedTrain)
summary(model)

projectedTrain$estimate <- predict(model,newdata=projectedTrain)
ScatterHist(projectedTrain,'estimate','y','Recovered 20 variable model versus truth (train)',
            smoothmethod='identity',annot_size=3)

trainrsq <- rsq(projectedTrain$estimate,projectedTrain$y)

projectedTest <- as.data.frame(predict(princ,dmTest),
                               stringsAsFactors = FALSE)
projectedTest$y <- dTestNTreatedXscaled$y
projectedTest$estimate <- predict(model,newdata=projectedTest)
testrsq <- rsq(projectedTest$estimate,projectedTest$y)
testrsq

proj <- extractProjection(2,princ)
# apply projection
projectedTrain <- as.data.frame(dmTrain %*% proj,
                                stringsAsFactors = FALSE)
projectedTrain$y <- dTrainNTreatedXscaled$y
# plot data sorted by principal components
ScatterHistN(projectedTrain,'PC1','PC2','y',
             "x scaled Data projected to first two principal components")
