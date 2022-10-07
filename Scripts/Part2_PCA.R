
library(factoextra)
library(missMDA)
library(corrplot)
library(ggfortify)
library(BayesFactor)
library(ggpubr)
library(dplyr)
library(vtable)
library(psych)
library(bclust)
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


#PCA of Metabolite variation


Metabolite.PCA <- prcomp(x[,1:43], center = T, scale. = T)


#produce eigenvalue plot

pdf("OUTPUT/Metabolite_eigen.pdf")
fviz_eig(Metabolite.PCA, choice = "eigenvalue", addlabels = T)
dev.off()


#produce variable loading plots
pdf("OUTPUT/Metabolite_VariableLoading_axes_1_2.pdf")
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
pdf("OUTPUT/Metabolite_IndividualLoading_axes_1_2.pdf")
  fviz_pca_ind(
    Metabolite.PCA,
    geom.ind = "point",
    pointshape = 19,
    #label = x$PlantID,
    col.ind = x$Genotype,
    mean.point = F
  )
dev.off()



#assess cos2 this is a metric of how well the variable is explaned by the component

var <- get_pca_var(Metabolite.PCA)


#create plot for Cos2
pdf("OUTPUT/Metabolite_cos2.pdf")
corrplot(
  var$cos2,
  method = "color",
  #p.mat = round(var$cos2, digits = 2),
  insig = "p-value",
  sig.level = 0.01,
  outline = T,
  #cl.lim = c(0, 1),
  is.corr = FALSE
)
dev.off()










