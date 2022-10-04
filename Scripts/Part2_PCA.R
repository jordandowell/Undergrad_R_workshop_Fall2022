
library(factoextra)
library(missMDA)
library(corrplot)
library(ggfortify)
library(BayesFactor)
library(ggpubr)
library(dplyr)
library(vtable)
library(psych)
#import data
#import data
data("gaelle")

#set seed for reproducibility
set.seed(150)

# rename data for ease of coding
x <- gaelle

View(x)


#create a vector of PlantID and Genotype

x<-cbind(x,rownames(x))
x<-cbind(x,gsub("\\..*","",as.character(x[,44])))

x<-as.data.frame(x)
x[,1:43]<-sapply(x[,1:43], as.numeric)

names(x)[44:45]<-c("PlantID","Genotype")
x[,45]<-as.factor(x[,45])


#PCA of flavenoids and phenolics


Metabolite.PCA <- prcomp(x[,1:43], center = T, scale. = T)


#produce eigenvalue plot

pdf("OUTPUT/FlavPhen_eigen.pdf")
fviz_eig(Metabolite.PCA, choice = "eigenvalue", addlabels = T)
dev.off()


#produce variable loading plots
pdf("OUTPUT/FlavPhen_VariableLoading_axes_1_2.pdf")
fviz_pca_var(
  Metabolite.PCA,
  axes = c(1, 2),
  label = "var",
  alpha.var = "cos2",
  col.var = "cos2",
  # Color by contributions to the PC
  gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
  repel = TRUE     # Avoid text overlapping
)
dev.off()

#produce plot of individuals labels by core 12
pdf("OUTPUT/FlavPhen_IndividualLoading_axes_1_2.pdf")
PCA.IND <-
  fviz_pca_ind(
    Metabolite.PCA,
    geom.ind = "point",
    pointshape = 19,
    #label = x$PlantID,
    col.ind = x$Genotype,
    mean.point = F
  )

#assess cos2

var <- get_pca_var(Metabolite.PCA)


#create plot for Cos2
pdf("OUTPUT/FlavPhen_cos2.pdf")
corrplot(
  var$cos2,
  method = "color",
  #p.mat = round(var$cos2, digits = 2),
  insig = "p-value",
  sig.level = 0.01,
  outline = T,
  cl.lim = c(0, 1),
  is.corr = FALSE
)
dev.off()




system(say "I'm Done")


bf <-
  ttestBF(formula = reformulate(
    termlabels = "Genotype",
    response = colnames(Comparison)[l]
  ),
  data = Comparison)

bf.frame <- data.frame(bf)
chains <- posterior(bf, iterations = 100000)
chains.frame <- data.frame(unclass(summary(chains)))





for (l in 1:43) {
  trait.range <-
    range(Comparison[, l])[2] - range(Comparison[, l])[1]
  
  Trait <- colnames(Comparison)[l]
  #bftest has a error but we need the loop to continue
  bf <-
    ttestBF(formula = reformulate(
      termlabels = "Genotype",
      response = colnames(Comparison)[l]
    ),
    data = Comparison)
  bf.frame <- data.frame(bf)
  chains <- posterior(bf, iterations = 100000)
  chains.frame <- data.frame(unclass(summary(chains)))
  #save results
  results <- data.frame()
  results <-
    cbind(
      bf.frame$bf,
      bf.frame$error,
      (chains.frame$statistics.Mean[2] / abs(trait.range)),
      (chains.frame$statistics.SD[2] / abs(trait.range)),
      (chains.frame$quantiles.2.5.[2] / abs(trait.range)),
      (chains.frame$quantiles.97.5.[2] / abs(trait.range))
    )
  row.names(results) <- paste(Trait,row.names(chains.frame[2,]))
  colnames(results) <-
    c("BF", "BFerror", "mean", "mean.SD", "quant.2.5", "quant.97.5")
  BF_TTest_results <- rbind(BF_TTest_results, results)
}
#add them to large results file

View(BF_TTest_results)

write.csv(ECOPHYS_results, "OUTPUT/Ecophys_violinplots_breed.csv")






#summary statistics!
describe.by()