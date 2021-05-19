# initialize before presentation
{
  setwd("C:/Users/rl02898/Documents/Multivariate Statistics/Final")
  # setwd("~/Documents/Multivariate")
  library("MASS")
  library("fastDummies")
  library("factoextra")
  library("cluster")
  
  # custom functions
  {
    pause = function(){
      if (interactive())
      {
        invisible(readline(prompt = "Press <Enter> to continue..."))
      }
      else
      {
        cat("Press <Enter> to continue...")
        invisible(readLines(file("stdin"), 1))
      }
    }
  }
}



###############################################################################
################################ DATA CLEANING ################################
###############################################################################
# first turn character variables into dummy variables
cars <- dummy_cols(Cars93, select_columns = c("Type",
                                              "AirBags",
                                              "DriveTrain",
                                              "Man.trans.avail",
                                              "Origin"))

# now eliminate redundant and extra unwanted variables
cars <- cars[,-c(1:4, 6, 9:10, 15:16, 26, 27, 40, 43)]

## small changes
cars$Luggage.room[is.na(cars$Luggage.room)] <- 0                       # change NA to 0
cars$Rear.seat.room[is.na(cars$Rear.seat.room)] <- 0                   # change NA to 0
cars$Cylinders <- as.integer(levels(cars$Cylinders))[cars$Cylinders]   # turning a factor into integer
cars$Cylinders[is.na(cars$Cylinders)] <- 0                             # change NA to 0

# fix empty spaces in variables names
colnames(cars)[23] <- "AirBags_DriverAndPassenger"
colnames(cars)[24] <- "AirBags_DriverOnly"

# let's get a look at the cleaned data
head(cars)



###############################################################################
######################## FORWARD VARIABLE SELECTION ###########################
###############################################################################
# original forward stepwise variable selection
{
  # N <- nrow(cars)                                 # number of rows in Cars93
  # k <- 10                                         # the number of times I will run forward selection
  # variable <- "MPG.highway"                       # the response variables
  # set.seed(1)
  # shuffle <- sample(N)                            # indexes to randomize training data
  # folds <- cut(1:N, breaks = k, labels = FALSE)   # creating folds for CV
  # X <- cbind(cars[,-which(colnames(cars) == variable)], cars[,which(colnames(cars) == variable)])
  # colnames(X)[ncol(X)] <- variable
  # threshold <- 2                                  # cut-off when model is not improving by .003
  # effect <- vector(mode = "list", length = k)     # an empty list that will contain the "sumSquares" of the model
  # models <- vector(mode = "list", length = k)     # an empty list that will contain the variables used
  # 
  # for (m in 1:k) {
  #   # separate Cars93 into a training set and validation set
  #   shuffledTrain <- shuffle[-which(folds == m)]
  #   shuffledValidate <- shuffle[which(folds == m)]
  #   shuffledXtrain <- X[shuffledTrain,]               # k-th fold training set
  #   shuffledXvalidate <- X[shuffledValidate,]         # k-th fold validation set
  #   p <- ncol(X) - 1
  #   sumSquarespast <- c(0, rep(NA, p-1))
  #   accepted <- c()
  #   # forward step-wise selection algorithm
  #   for (i in 1:p) {
  #     sumSquares <- rep(NA, p)
  #     if (i == 1) {
  #       leftOver <- colnames(X[,-(p+1)])
  #     }
  #     else {
  #       leftOver <- colnames(X[,-(p+1)])[-which(colnames(X) %in% accepted)]
  #     }
  #     # loop over each variable not already in the model
  #     for (j in leftOver) {
  #       if (length(accepted) == 0) {
  #         formula <- as.formula(paste(variable, "~ ", j))
  #         accepted <- c()
  #       }
  #       else {
  #         formula <- as.formula(paste(variable, "~", paste(accepted, collapse = " + "), " + ", j))
  #       }
  #       model <- lm(formula, data = shuffledXtrain)
  #       sumSquares[which(colnames(shuffledXvalidate) == j)] <- sum((as.matrix(cbind(1,shuffledXvalidate[,c(accepted, j)]))%*%model$coefficients - shuffledXvalidate[, p + 1])^2)
  #     }
  #     accepted <- c(accepted, colnames(X)[which.min(sumSquares)])   # adding the best predictor variable to the model
  #     if (abs(min(sumSquares[!is.na(sumSquares)]) - min(sumSquarespast[!is.na(sumSquarespast)])) < threshold) {   # if no improvement, break loop
  #       print(min(sumSquares[!is.na(sumSquares)]))
  #       effect[[m]] <- min(sumSquares[!is.na(sumSquares)])
  #       print(accepted)
  #       models[[m]] <- accepted
  #       break
  #     }
  #     sumSquarespast <- sumSquares
  #   }
  # }
  # 
  # # After 10 loops we have some lists of optimal variables for models
  # modelVariables <- c()
  # for (i in 1:k) {     # put variables into one list
  #   modelVariables <- c(modelVariables, models[[i]])
  # }
  # modelVariables <- c(variable, unique(modelVariables))
}

# built-in stepwise variable selection, both ways
{
  model <- lm(MPG.city + MPG.highway ~ ., data = cars)
  FSS <- stepAIC(model, direction = "both")
  modelVariables <- all.vars(FSS$terms[[3]])
  modelVariables
}



###############################################################################
########################## CROSS VALIDATION FOR PCs ###########################
###############################################################################
## Cross validation for best number of PCs
{
  r <- "MPG.highway"
  ind <- c(seq(ncol(cars))[-which(colnames(cars) %in% modelVariables)], which(colnames(cars) %in% r))
  
  CVerror <- rep(0, (ncol(cars)-length(ind)-1))
  for (i in 1:(ncol(cars)-length(ind))) {
    dist <- rep(NA, nrow(cars))
    for (p in 1:nrow(cars)) {
      X <- cars[-p,-ind]
      Y <- cars[-p,which(colnames(cars) == r)]
      PCA <- prcomp(X)
      PC <- as.matrix(X) %*% PCA$rotation[,1:i]
      lm <- lm(Y ~ PC)
      pred <- as.matrix(cars[p,])
      prediction <- c(1, pred[-ind] %*% PCA$rotation[,1:i]) %*% lm$coefficients
      dist[p] <- (cars[p, which(colnames(cars) == r)]-prediction)^2
    }
    CVerror[i] <- sum(dist)
  }
  pdf(file = "CV.pdf", height = 4)
  plot(CVerror, main = "Cross Validation Error", xlab = "Number of PCs")
  dev.off()
  plot(CVerror, main = "Cross Validation Error", xlab = "Number of PCs")
  print(which.min(CVerror))
  k <- 6
  print(length(which(colnames(cars) %in% modelVariables)))
}


###############################################################################
## shhhh, this is for detecting outliers, leave this alone
{
  testCars <- cars[,1:16]
  mu <- colMeans(testCars)
  S <- cov(testCars)
  T2 <- rep(0, nrow(testCars))
  for (i in 1:nrow(testCars)) {
    T2[i] <- as.matrix(testCars[i,] - mu)%*%solve(S)%*%t(as.matrix(testCars[i,] - mu))
  }
  crit <- qchisq(0.95, ncol(testCars))
  outliers <- Cars93[which(T2 > crit), 27]
}
###############################################################################



###############################################################################
###################### PRINCIPAL COMPONENT ANALYSIS ###########################
###############################################################################
r <- "MPG.city"       # response
r <- "MPG.highway"    # response

## Analysis
for (p in 1:nrow(cars)) {
  X <- cars[-p,-ind]
  Y <- cars[-p,which(colnames(cars) == r)]
  PCA <- prcomp(X)
  PC <- as.matrix(X) %*% PCA$rotation[,1:k]
  lm <- lm(Y ~ PC)
  pred <- t(as.matrix(cars[p,]))
  PCpred <- c(1, t(pred[-ind]) %*% PCA$rotation[,1:k])
  if (sum(outliers %in% Cars93[p,27]) > 0) {
    print(paste(as.character(Cars93[p, 27]), "***"))
  }
  else {
    print(as.character(Cars93[p, 27]))
  }
  print("prediction")
  print(as.vector(round(PCpred %*% lm$coefficients)))
  print("true")
  print(cars[p, which(colnames(cars) == r)])
  pause()
}


## Predict my car
response <- c("MPG.city", "MPG.highway")      # response

# my cars specifications
myCar <- c(
  20.8,   # price
  24,     # city MPG
  34,     # highway MPG
  4,      # cylinders
  1.8,    # engine size
  170,    # horsepower
  4800,   # RPM
  18.5,   # fuel tank capacity
  5,      # passengers
  191.6,  # length
  110.4,  # wheelbase
  72.2,   # width
  36.4,   # turn circle
  39.1,   # rear seat room
  15.9,   # luggage room
  3166,   # weight
  0,      # compact
  0,      # large
  1,      # midsize
  0,      # small
  0,      # sporty
  0,      # van
  1,      # driver and passenger airbags
  0,      # driver only airbags
  0,      # no airbags
  0,      # 4WD
  1,      # FWD
  0,      # RWD
  0,      # manual transmission available
  0       # origin
)

for (r in response) {
  X <- cars[,-ind]
  Y <- cars[,which(colnames(cars) == r)]
  PCA <- prcomp(X)
  PC <- as.matrix(X) %*% PCA$rotation[,1:k]
  lm <- lm(Y ~ PC)
  PCpred <- c(1, t(myCar[-ind]) %*% PCA$rotation[,1:k])
  print(round(PCpred %*% lm$coefficients))
}

# true MPG.city is    24
# true MPG.highway is 34
# great!



###############################################################################
###################### MULTIVARIATE LINEAR REGRESSION #########################
###############################################################################
 # make a prediction using multivariate multiple linear regression
{
  r <- c("MPG.city", "MPG.highway")    # responses
  ind <- c(seq(ncol(cars))[-which(colnames(cars) %in% modelVariables)], which(colnames(cars) == r))
  
  for (p in 1:nrow(cars)) {
    X <- cars[-p,-ind]
    Y <- cars[-p,which(colnames(cars) %in% r)]
    PCA <- prcomp(X)
    PC <- as.matrix(X) %*% PCA$rotation[,1:k]
    lm <- lm(as.matrix(Y) ~ PC)
    pred <- t(as.matrix(cars[p,]))
    PCpred <- c(1, t(pred[-ind]) %*% PCA$rotation[,1:k])
    if (sum(outliers %in% Cars93[p,27]) > 0) {
      print(paste(as.character(Cars93[p, 27]), "***"))
    }
    else {
      print(as.character(Cars93[p, 27]))
    }
    print("prediction")
    print(as.vector(round(PCpred %*% lm$coefficients)))
    print("true")
    print(as.vector(unlist(cars[p, which(colnames(cars) %in% r)])))
    pause()
  }
}


# Make a prediction for origin using logistic regression
{
  r <- "Origin_USA"
  ind <- c(seq(ncol(cars))[-which(colnames(cars) %in% modelVariables)], which(colnames(cars) == r))
  
  for (p in 1:nrow(cars)) {
    X <- cars[-p,-ind]
    Y <- cars[-p,which(colnames(cars) == r)]
    PCA <- prcomp(X)
    PC <- as.matrix(X) %*% PCA$rotation[,1:k]
    data <- as.data.frame(cbind(PC, Y))
    formula <- as.formula(paste("Y ~", paste(colnames(data[,-ncol(data)]), collapse = " + ")))
    lm <- glm(formula, data = data, family = binomial)
    pred <- as.matrix(cars[p,])
    PCpred <- as.data.frame(pred[-ind] %*% PCA$rotation[,1:k])
    if (sum(outliers %in% Cars93[p,27]) > 0) {
      print(paste(as.character(Cars93[p, 27]), "***"))
    }
    else {
      print(as.character(Cars93[p, 27]))
    }
    print("prediction")
    print(as.vector(round(predict(lm, newdata = PCpred, type = "response"))))
    print("true")
    print(cars[p, which(colnames(cars) == r)])
    pause()
  }
}



###############################################################################
##################################### TESTING #################################
###############################################################################
## Mahalanobis Test for Normality
{
  data <- cars[,-c(17:ncol(cars))]
  for (i in 1:ncol(data)) {
    data[which(data[,i] == -Inf),i] <- mean(data[-which(data[,i] == -Inf),i])
  }
  mu <- colMeans(data)
  cov <- cov(data)
  mah <- sort(mahalanobis(data, center = mu, cov = cov))
  chi <- qchisq((1:nrow(cars)-0.5)/(nrow(cars)), df = ncol(data))
  pdf(file = "mahalanobis.pdf", height = 6)
  plot(chi, mah, main = "Mahalanobis Distances Plot")
  dev.off()
  plot(chi, mah, main = "Mahalanobis Distances Plot")
}


## Hypothesis Test
## Is my car similar to the other cars
{
  testCars <- cars[,1:16]
  testPred <- myCar[1:16]
  n <- nrow(testCars)
  p <- ncol(testCars)
  xbar <- colMeans(testCars)
  S <- cov(testCars)
  mu0 <- testPred
  test <- n*t(xbar - mu0) %*% solve(S) %*% (xbar - mu0)
  crit <- (n-1)*p/(n-p)*qf(0.05, p, n-p, lower.tail = FALSE)
  print("Test Statistic:")
  print(test)
  print("Critical Value:")
  crit
}
## Reject the null, my car is very different than an average car


## simultaneous confidence intervals
## shadow's of confidence region
{
  alpha <- 0.05
  CIs1 <- matrix(0, nrow = 16, ncol = 2)
  for (i in 1:16) {
    CIs1[i,1] <- xbar[i] - sqrt((n-1)*p/(n-p)*qf(alpha, p, n-p, lower.tail = FALSE))*sqrt(S[i,i]/n)
    CIs1[i,2] <- xbar[i] + sqrt((n-1)*p/(n-p)*qf(alpha, p, n-p, lower.tail = FALSE))*sqrt(S[i,i]/n)
  }
  CIs1 <- as.data.frame(cbind(CIs1, mu0))
  colnames(CIs1) <- c("LB", "UB", "     My Car")
  rownames(CIs1) <- colnames(testCars)
  sink(file = "CIs1.txt")
  CIs1
  sink()
  CIs1
}


# Bonferroni simultaneous confidence intervals
{
  CIs2 <- matrix(0, nrow = 16, ncol = 2)
  for (i in 1:16) {
    CIs2[i,1] <- xbar[i] + qt(alpha/2/p, n-1)*sqrt(S[i,i]/n)
    CIs2[i,2] <- xbar[i] - qt(alpha/2/p, n-1)*sqrt(S[i,i]/n)
  }
  CIs2 <- as.data.frame(cbind(CIs2, mu0))
  colnames(CIs2) <- c("LB", "UB", "     My Car")
  rownames(CIs2) <- colnames(testCars)
  sink(file = "CIs2.txt")
  CIs2
  sink()
  CIs2
}


## Outlier detection
{
  T2 <- rep(0, nrow(testCars))
  for (i in 1:nrow(testCars)) {
    T2[i] <- as.matrix(testCars[i,] - mu)%*%solve(S)%*%t(as.matrix(testCars[i,] - mu))
  }
  crit <- qchisq(0.95, ncol(testCars))
  pdf(file = "outliers.pdf", height = 6)
  plot(1:length(T2), T2,
       xlab = "Car Index",
       ylab = expression("T"^2~"value"),
       ylim = c(min(T2, crit), max(T2, crit)),
       main = "Outlier Detection with Chi-Square Quantile")
  abline(h = crit)
  dev.off()
  plot(1:length(T2), T2,
       xlab = "Car Index",
       ylab = expression("T"^2~"value"),
       ylim = c(min(T2, crit), max(T2, crit)),
       main = "Outlier Detection with Chi-Square Quantile")
  abline(h = crit)
  Cars93[which(T2 > crit), 27]
}



###############################################################################
############################### CLUSTER ANALYSIS ##############################
###############################################################################
# variables for clustering
{
  N <- nrow(cars)                                 # number of rows in Cars93
  k <- 10                                         # the number of times I will run forward selection
  variable <- "MPG.highway"                       # the response variables
  set.seed(1)
  shuffle <- sample(N)                            # indexes to randomize training data
  folds <- cut(1:N, breaks = k, labels = FALSE)   # creating folds for CV
  X <- cbind(cars[,-which(colnames(cars) == variable)], cars[,which(colnames(cars) == variable)])
  colnames(X)[ncol(X)] <- variable
  threshold <- 2                                  # cut-off when model is not improving by .003
  effect <- vector(mode = "list", length = k)     # an empty list that will contain the "sumSquares" of the model
  models <- vector(mode = "list", length = k)     # an empty list that will contain the variables used
  
  for (m in 1:k) {
    # separate Cars93 into a training set and validation set
    shuffledTrain <- shuffle[-which(folds == m)]
    shuffledValidate <- shuffle[which(folds == m)]
    shuffledXtrain <- X[shuffledTrain,]               # k-th fold training set
    shuffledXvalidate <- X[shuffledValidate,]         # k-th fold validation set
    p <- ncol(X) - 1
    sumSquarespast <- c(0, rep(NA, p-1))
    accepted <- c()
    # forward step-wise selection algorithm
    for (i in 1:p) {
      sumSquares <- rep(NA, p)
      if (i == 1) {
        leftOver <- colnames(X[,-(p+1)])
      }
      else { 
        leftOver <- colnames(X[,-(p+1)])[-which(colnames(X) %in% accepted)]
      }
      # loop over each variable not already in the model
      for (j in leftOver) {
        if (length(accepted) == 0) {
          formula <- as.formula(paste(variable, "~ ", j))
          accepted <- c()
        }
        else {
          formula <- as.formula(paste(variable, "~", paste(accepted, collapse = " + "), " + ", j))
        }
        model <- lm(formula, data = shuffledXtrain)
        sumSquares[which(colnames(shuffledXvalidate) == j)] <- sum((as.matrix(cbind(1,shuffledXvalidate[,c(accepted, j)]))%*%model$coefficients - shuffledXvalidate[, p + 1])^2)
      }
      accepted <- c(accepted, colnames(X)[which.min(sumSquares)])   # adding the best predictor variable to the model
      if (abs(min(sumSquares[!is.na(sumSquares)]) - min(sumSquarespast[!is.na(sumSquarespast)])) < threshold) {   # if no improvement, break loop
        # print(min(sumSquares[!is.na(sumSquares)]))
        effect[[m]] <- min(sumSquares[!is.na(sumSquares)])
        # print(accepted)
        models[[m]] <- accepted
        break
      }
      sumSquarespast <- sumSquares
    }
  }
  
  # After 10 loops we have some lists of optimal variables for models
  modelVariables <- c()
  for (i in 1:k) {     # put variables into one list
    modelVariables <- c(modelVariables, models[[i]])
  }
  modelVariables <- c(variable, unique(modelVariables))
}

# non-hierarchial clustering
data <- cars[,-seq(ncol(cars))[-which(colnames(cars) %in% modelVariables)]]
rownames(data) <- Cars93[,27]
set.seed(4)
kmlist <- vector("list", ncol(data))
for(k in 2:(ncol(data)+1)){
  kmlist[[k-1]] <- kmeans(data, k, nstart=5000)
}

## Find the number of clusters, K
{
tot.withinss <- sapply(kmlist, function(x) x$tot.withinss)
pdf(file = "scree.pdf", height = 6)
plot(2:(ncol(data)+1),
     tot.withinss/kmlist[[1]]$totss,
     xlab="Clusters",
     ylab="SSW / SST",
     main="Scree Plot for Clustering")
dev.off()
plot(2:(ncol(data)+1),
     tot.withinss/kmlist[[1]]$totss,
     xlab="Clusters",
     ylab="SSW / SST",
     main="Scree Plot for Clustering")
}
K <- 15

## Plot the clusters
{
  pdf(file = "clusterPlot.pdf", width = 8, height = 6)
  fviz_cluster(kmlist[[K-1]], data = data)
  dev.off()
  fviz_cluster(kmlist[[K-1]], data = data)
}

# look at the clusters
sort(kmlist[[K-1]]$cluster)
# cluster 3
sort(kmlist[[K-1]]$cluster)[which(sort(kmlist[[K-1]]$cluster) == 3)]
# cluster 10
sort(kmlist[[K-1]]$cluster)[which(sort(kmlist[[K-1]]$cluster) == 10)]
# cluster 12
sort(kmlist[[K-1]]$cluster)[which(sort(kmlist[[K-1]]$cluster) == 12)]


# Hierarchial Clustering
distances <- matrix(0, nrow = nrow(cars), ncol = nrow(cars))
for (i in 1:nrow(distances)) {
  for (j in 1:ncol(distances)) {
    if (i > j) {
      distances[i,j] <- sqrt(sum((cars[i,] - cars[j,])^2))
    }
    else {
      distances[i,j] <- 0
    }
  }
}
rownames(distances) <- colnames(distances) <- Cars93[,27]
d <- as.dist(distances)
cluster <- hclust(d, method = "complete")
pdf(file = "complete.pdf", height = 6)
plot(cluster, cex = 0.3, hang = -1)
abline(h = 650)
dev.off()
plot(cluster, cex = 0.3, hang = -1)
abline(h = 650)
sort(cutree(cluster, k = 15))

