library(bclust)

#import data
data("gaelle")

#set seed for reproducibility
 set.seed(150)

 # rename data for ease of coding
 x <- gaelle
 
 
 #add a set of label ids
 x.id <- rep(1:14, c(3, rep(4, 13)))
 
 #the following lines are where we are setting priors in this case 
 #we are using the data to calculate the maximum marginal likelihood
 
 #calculate the mean, corrected sum of squares for calculation of the loglikelihood
 meansumsq <- meancss (x, x.id)
 
 #a shortened function to calculate -loglikelihood
 #in this instance we are assuming gaussian noise in our error estimations. another estimate could be laplace if you think the 
 optimfunc <- function(phi) {
    -loglikelihood(x.mean = meansumsq$mean, x.css = meansumsq$css,
                     repno = meansumsq$repno, transformed.par = phi, var.select = FALSE) }
 
 
 
xinit.tpar <- optim(rep(0, 5), optimfunc, method = "BFGS")$par

#adding variable selection to make a different parameter set. 
optimfunc <- function(phi) {
   -loglikelihood(x.mean = meansumsq$mean, x.css = meansumsq$css,
                     repno = meansumsq$repno, transformed.par = c(xinit.tpar[1:4], phi)) }


x.tpar <- c(xinit.tpar[1:4], optim(rep(0, 2), optimfunc,
                                       method = "BFGS")$par)
 
 

#create genotype labels
x.labels <- c("ColWT", "d172", "d263", "isa2", "sex4", "dpe2", "mex1",
                 "sex3", "pgm", "sex1", "WsWT", "tpt", "RLDWT", "ke103")


#first question if we dont know what colWT is and assume its a seed 
#we found on the ground that corresponds to one of our known mutants
#what do we think it could be 

#predictive training w/o variable selection 
bdiscrim(training = x[-(1:3), ], training.id = (x.id[-(1:3)] - 1),
             training.labels = x.labels[-1], predict = x[1:3, ],
             transformed.par = xinit.tpar, var.select = FALSE)$probs * 100
#predictive training w/ variable selection 
bdiscrim(training = x[-(1:3), ], training.id = (x.id[-(1:3)] - 1),
             training.labels = x.labels [-1], predict = x[1:3, ],
             transformed.par = x.tpar)$probs * 100
#what does variable selection provide? 


#lets futher examine the data.
#first off is there evidence that colWT is highly similar to WsWT?


#make a cluster object that we will manipulate 
#if there are a lot of traits or samples this may take longer

bclust.obj <- bclust(x, rep.id = x.id, transformed.par = x.tpar,
                            labels = x.labels)

#create a tree
plot(bclust.obj)

#apply a dashed line where the higher posterior probability of clustering 
#place it 2% the height of the tree lower for visualization and interpretation
abline(h = (bclust.obj$cut-(bclust.obj$height[length(bclust.obj$height)]*.02)), lty = 2, col = "red", lwd = 3)


#how many suggested clusters based on evidence?


#which metabolites are useful to identify mutants?
#Positive bayes factors are highlighted in red!

viplot(imp(bclust.obj)$var,xlab = colnames(x), col = as.numeric(imp(bclust.obj)$var > 0) * 2)

#lets see our clustering relative to our variables

plot.new()

dptplot(bclust.obj, scale = 10, horizbar.plot = F,
               varimp = imp(bclust.obj)$var, horizbar.distance = 5, dendrogram.lwd = 2)
plot.new()
#alternative visualization
#red means high bayes factor
ditplot(bclust.obj, horizbar.plot = F,
        varimp = imp(bclust.obj)$var, horizbar.distance = 5, dendrogram.lwd = 2)


#WAIT this suggests more clusters than we previously identified what happened?

#lets check what the cut of the graph is 
bclust.obj$cut

# & check the height of each grouping & logposterior 
bclust.obj$height

bclust.obj$logposterior

#these are really close alternative criteria can be used to estimate the cut

