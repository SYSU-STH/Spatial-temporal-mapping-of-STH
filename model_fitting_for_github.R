library(INLA)
library(maps)
library(foreign)
library(raster)
library(maptools)
library(sp)
library(rgeos)
library(ggplot2)
library(MASS)
library(MBA)
library(fields)
library(combinat)


#### 1 read data ####
setwd('/Users/data_path/') #data path

data <- read.csv("Hookworm_total.csv") #i.e. Hookworm


#### 2 model setting ####
#2.1 setting####
validate <- c("NumPositive") #the number of positive individuals
exam <- c("NumExamine") #the number of examined individuals
data$Intercept <- 1

#2.2 mesh####
a <- shapefile("guangdong map/guangdong.shp") 
IDs <- a$ID_1
outer_border <- unionSpatialPolygons(a, IDs) 
boundaries <- fortify(outer_border)
boundaries  <- boundaries[,1:2]
colnames(boundaries) <- c("x","y")
plot(boundaries)

#coordinates of points for making mesh
coords <- cbind(data$longitude,data$latitude)
colnames(coords) <- c("x","y")

dist <- dist(as.matrix(coords), method = "euclidean", diag = FALSE, 
           upper = FALSE, p = 2)
kappa <- 3*2/max(dist)
mesh <- inla.mesh.create.helper(points = as.matrix(coords),
                                points.domain = cbind(boundaries$x,boundaries$y),
                                max.edge = c(0.2/kappa,0.4/kappa),
                                cutoff = 0.045/kappa, 
                                offset = c(0.1,0.5),
                                plot.delay = NULL)

plot(mesh)
points(boundaries,lwd = 3)
points(as.matrix(coords),pch = 20,cex = 1,col = 2)

#2.3 SPDE model####
#define SPDE model
spde <- inla.spde2.matern(mesh = mesh,alpha = 2) 
hyper = spde$f$hyper.default

#construct spatial-temporal model
field.indices <- inla.spde.make.index("field", mesh$n, 
                                     n.group = length(table(data$surveydate)))

#define projection matrix
A.est <- inla.spde.make.A(mesh,loc = as.matrix(coords), group = data$surveydate) 

#effects in the stack
predictor <- c("sex","age1","age2","age3","age4","surveydate1","surveydate2",
               "shii","shumidity","salt1","salt2","loc_id")
stack.est <- inla.stack(data = list(validate = data[,validate]),
                       A = list(A.est, 1),
                       effects = list(c(field.indices, list(Intercept = 1)),
                                      list(data[predictor])),
                       tag = "est") 


#### 3 model fitting ####
for1 <- c("-1","Intercept","f(loc_id,model = 'iid')",
        "f(field,model = spde,hyper = hyper,group = field.group, 
        control.group = list(model = 'ar1'))")
formula <- as.formula(paste("validate ~ ", paste(paste(for1, collapse = "+"), 
                                                 paste(predictor, collapse = "+"),
                                                 sep = "+")))
stack <- inla.stack(stack.est)

result <- inla(formula,data = inla.stack.data(stack,spde = spde),family = "binomial",
               Ntrials = data[,exam],control.family = list(link = "logit"),
               control.predictor = list(compute = TRUE,A = inla.stack.A(stack),
                                      quantiles = c(0.025,0.5,0.975),link = link),
               control.compute = list(dic = TRUE,config = TRUE,cpo = TRUE),
               control.inla = list(int.strategy = "eb"),verbose = TRUE)

#to obtain summaries of the spatial parameters
result.field <- inla.spde.result(result,"field",spde,do.transform = TRUE)

#eg for the range
inla.qmarginal(c(0.025,0.5,0.975),result.field$marginals.range.nominal[[1]])

#for spatial variance
inla.qmarginal(c(0.025,0.5,0.975), result.field$marginals.variance.nominal[[1]])

#for non-spatial variance
sigma2nonspatial = inla.tmarginal(function(x)1/x,
                                  result$marginals.hyper$"Precision for loc_id")
inla.qmarginal(c(0.025,0.5,0.975),sigma2nonspatial)

summary(result)
log.score <- -mean(log(result$cpo$cpo))
log.score
dic <- result$dic$dic[1]
dic

save.image("model_Hookworm.RData")


#### 4 prediction ####
#4.1 get 500 samples from posteriori distribution####
set.seed(100)
nsample = 500
samp <- inla.posterior.sample(n = 500, result, use.improved.mean = TRUE)
h <- t(sapply(samp, function(x) x$latent))

#get random field
a <- paste(c("field.",1),collapse = "")
b <- paste(c("field.",dim(A.est)[2]),collapse = "")
beg <- grep(a, rownames(samp[[1]]$latent))
end <- grep(b, rownames(samp[[1]]$latent))
samp.field <- t(h[,beg:end])

#get coefficient 
beg1 <- grep("Intercept", rownames(samp[[1]]$latent))
end1 <- grep("salt2", rownames(samp[[1]]$latent))
samp.coef <- t(h[,beg1:end1])

#get location exchangeable hyperparameter 
hy <- t(sapply(samp, function(x) x$hyperpar))
beg2 <- grep("loc_id", colnames(hy))
samp.precision <- t(hy[,beg2])

#4.2 prediction####
#predicted data (i.e. first national survey)
path <- "/Users/result/Hookworm/pre_1/" #path for result
data.whole <- read.csv("data_whole_1.csv")

data.whole$surveydate1 <- 1
data.whole$surveydate2 <- 0
data.whole$Intercept <- 1

#predict age- and gender-specific prevalence
for (agegr in 1:5) {
  for (k in 1:2) {
    data.whole$sex <- ifelse(k = 2,1,0)
    data.whole$age1 <- ifelse(agegr = 1,1,0)
    data.whole$age2 <- ifelse(agegr = 2,1,0)
    data.whole$age3 <- ifelse(agegr = 4,1,0)
    data.whole$age4 <- ifelse(agegr = 5,1,0)
    
    coords.predict <- cbind(data.whole$longitude,data.whole$latitude)
    plot(mesh)
    points(boundaries,lwd = 3)
    points(as.matrix(coords),pch = 20,cex = 1,col = 2)
    points(as.matrix(coords.predict),pch = 20,cex = 1,col = 3)
   
    #prediction of points
     A.pre = inla.spde.make.A(mesh, loc = matrix(rep(t(as.matrix(coords.predict)),3),
                                              ncol = ncol(as.matrix(coords.predict)),
                                              byrow = TRUE),
                             group = rep(1:3,each = nrow(data.whole)))
    
    pred.field <- A.pre%*%samp.field
    a1 = 0*dim(data.whole)[1]+1
    a2 = 1*dim(data.whole)[1]
    pred.field1 <- pred.field[a1:a2,]
    env <- cbind(rep(1,nrow(data.whole)),data.whole[,c("sex","age1","age2","age3","age4",
                                                     "surveydate1","surveydate2",
                                                     "shii","shumidity","salt1","salt2")])
    N <- nrow(data.whole)
    pred.loc <- matrix(NA,nrow = N,ncol = 500)
    for (j in 1:500){
     pred.loc[,j] <- rnorm(N, mean = 0, sd = sqrt(1/samp.precision[j]))
    }
    pred.linear <- as.matrix(env)%*%samp.coef+pred.field1+pred.loc
    
    #transfer to prevalence
    pred.prev <- exp(pred.linear)/(1+exp(pred.linear))
    pred.prev <- as.matrix(pred.prev)
    pred.prev <- data.frame(pred.prev)
    data.whole <- data.frame(data.whole)
    pred.samp <- cbind(data.whole$longitude,data.whole$latitude,pred.prev,data.whole[,16+k+2*(agegr-1)])
    pred.samp <- data.frame(pred.samp)
    
    write.dta(pred.samp,paste0(path,'pred.samp-age',agegr,'-sex',k,'.dta'))
    
    #summary low, median, up, sd
    library(matrixStats)
    n <- nrow(data.whole)
    pred.prev <- as.matrix(pred.prev)
    lower <- rowQuantiles(pred.prev, probs = 0.025)
    upper <- rowQuantiles(pred.prev, probs = 0.975)
    median <- rowQuantiles(pred.prev, probs = 0.5)
    stdev <- rowSds(pred.prev)
    results <- cbind(coords.predict,lower,median,upper,stdev,data.whole[,16+k+2*(agegr-1)])
    colnames(results) <- c("long","lat","lower","median","upper","stdev","pop")
    results <- data.frame(results)
    results$NumPos <- results$median*results$pop
    write.dbf(results,paste0(path,'result-age',agegr,'-sex',k,'.dbf'))
  }
}

#predict prevalence for whole people
pred.pos <- matrix(0,nrow = nrow(data.whole),ncol = nsample)
for (k in 1:2) {
  print(k)
  for(agegr in 1:5){
    print(agegr)
    dta <- read.dta(paste0(path,'pred.samp-age',agegr,'-sex',k,'.dta'))
    pred.age <- apply(dta[,3:(nsample+2)],2,FUN = function(x) x*dta[,503])
    pred.pos <- pred.pos+pred.age
  }
}

pred.pos <- data.frame(pred.pos)
pred.prev <- apply(pred.pos,2,FUN = function(x) x/data.whole$popc)
pred.prev <- data.frame(pred.prev)
data.whole <- data.frame(data.whole)
pred.samp <- cbind(data.whole$longitude,data.whole$latitude,pred.prev,data.whole$popc)
pred.samp <- data.frame(pred.samp)
write.dta(pred.samp,paste0(path,"pred.samp-all.dta"))

# summary low, median, up, sd#
pred.prev <- as.matrix(pred.prev)
lower <- rowQuantiles(pred.prev, probs = 0.025)
upper <- rowQuantiles(pred.prev, probs = 0.975)
median <- rowQuantiles(pred.prev, probs = 0.5)
stdev <- rowSds(pred.prev)
results <- cbind(data.whole$long,data.whole$lat,lower,median,upper,stdev,data.whole$popc)
colnames(results) <- c("long","lat","lower","median","upper","stdev","popc")
results <- data.frame(results)
write.dbf(results,paste0(path,"predict_result-all.dbf"))

save.image(paste(path,'gd_prevalence1.RData'))