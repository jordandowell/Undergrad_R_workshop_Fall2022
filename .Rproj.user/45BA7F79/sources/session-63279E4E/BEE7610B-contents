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


#add a set of label ids
x.id <- rep(1:14, c(3, rep(4, 13)))



#lets visualize the data as violin plots by genotype



#lets look at 1 specific x18
i<-1
#this will make violin plots of all metabolites by genotypes 
pdf("OUTPUT/ViolinPlot_metabolite_Genotype.pdf")
for (i in 1:43) {
  print(
    ggviolin(
      y = colnames(x)[i],
      x = "Genotype",
      add = c("jitter", error.plot = "crossbar"),
      data = x,
      draw_quantiles = c(0.025, 0.5, 0.975)
    )
  )
}
dev.off()

#now lets test if there are differences among genotypes 
#bayes factor analysis

#create a for loop for ANOVAresults of Metabolites by Genotype

#empty data frame 
BF_ANOVA_results <- data.frame()


#for loop by each column e.g each metabolite  
for (i in 1:43) {
 #find the range of the metabolite for this cycle of for loop  
   Metabolite.range <-
    range(x[, i])[2] - range(x[, i])[1]
  #find the name of the metabolite for this cycle of for loop 
  Metabolite <- colnames(x)[i]
  ##use a bayesian form of anova to assess if genotype has an effect
  bf <-
    anovaBF(formula = reformulate(
      termlabels = "Genotype",
      response = colnames(x)[i]
    ),
    data = x)
  #if you want more confidence in the estimation of the bayes factor add the following line 
  #bf<- recompute(bf, iterations =100000)
  #lets get the posterior data
  bf.frame <- data.frame(bf)
  chains <- posterior(bf, iterations = 100000)
  chains.frame <- data.frame(unclass(summary(chains)))
  #save results
  
  results <- data.frame()
  results <-
    cbind(
      bf.frame$bf,
      bf.frame$error,
      #we standardize the coeficienty by the range of metabolites 
      #so that coeficients are comparable among metabolites 
      (chains.frame$statistics.Mean[2] / abs(Metabolite.range)),
      (chains.frame$statistics.SD[2] / abs(Metabolite.range)),
      (chains.frame$quantiles.2.5.[2] / abs(Metabolite.range)),
      (chains.frame$quantiles.97.5.[2] / abs(Metabolite.range))
    )
  row.names(results) <- Metabolite
  colnames(results) <-
    c("BF", "BFerror", "mean", "mean.SD", "quant.2.5", "quant.97.5")
  BF_ANOVA_results <- rbind(BF_ANOVA_results, results)
}
#take a second to view the results 
#reference Kass & raferty to interepret bayes factor
#bayes factor == the probability that the hypothesis is more likely than the null 
# do genotypes vary among all metabolites?
View(BF_ANOVA_results)

#lets begin a posthoc analysis we will do a bayesian T-test on all genotype x metabolite combinations

#explicit question asked by the model. what is the probability that Data A comes from the same distribution as DataB?
#post-hoctest

BF_TTest_results <- data.frame()


#subset data by genotypes 

subset_of_x <- split(x, x$Genotype)
#pairwise for loop for genotypes 
#set values for intial loop 
i<-1
j<-1
l<-1

for (i in 1:length(subset_of_x)) {
  
  for(j in 1:length(subset_of_x)){
    #no need to compare to two of the same genotypes 
    if(j != i){print(paste(names(subset_of_x[i]),names(subset_of_x[j]),"Starting..."))
 #   else{}
  #  }}
    #run comparisons of each metabolite 
    Comparison<-rbind(subset_of_x[[i]],subset_of_x[[j]])
    #drop extra levels so bayes factor doesnt try to estimate for the absence factors 
    Comparison[,45]<-droplevels.factor(Comparison[,45])
    for (l in 1:43) {
      
     # print(paste(names(subset_of_x[i]),names(subset_of_x[j]),colnames(Comparison)[l],"Starting..."))
      trait.range <-
        range(Comparison[, l])[2] - range(Comparison[, l])[1]
      
      Trait <- colnames(Comparison)[l]
      #begin for loop for #ttest
      #skip_to_next <- FALSE
      #try catch is a good addition when certain combinations may form an error 
     
      #tryCatch( 
      #this section is running our T-test
        bf <-
          ttestBF(formula =  reformulate(
            termlabels = "Genotype",
            response = colnames(Comparison)[l]
          ),
          data = Comparison)#, error= function(e){ skip_to_next <<- TRUE})
      #this is what happens if theres an error 
     # if(skip_to_next) { next } 
      #save our results 
        #we are sampling the posterior to get an idea of the distribution of effects 
        #and out confidence is the estimates
      bf.frame <- data.frame(bf)
      chains <- posterior(bf, iterations = 10000)
      chains.frame <- data.frame(unclass(summary(chains)))
      #save results 
      #in the current form combinations with errors have the same values as the previous line 
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
      
      #add out results for this iteration to our growing results table 
      BF_TTest_results <- rbind(BF_TTest_results, results)
     # print(paste(names(subset_of_x[i]),names(subset_of_x[j]),colnames(Comparison)[l],"Finished..."))
      
    }
    
    
    print(paste(Trait,row.names(chains.frame[2,])))
    print(paste(names(subset_of_x[i]),names(subset_of_x[j]),"Finished..."))
    }else{next}
  }
}



#lets check some results

View(BF_TTest_results)
dim(BF_TTest_results)
#at this point you could write results to a csv.



#custom function for a Volcano plot 



JordanVolcano<-function (toptable, lab, x, y, selectLab = NULL, xlim = c(min(toptable[[x]], 
                                                                             na.rm = TRUE) - 1, max(toptable[[x]], na.rm = TRUE) + 1), 
                         ylim = c(0, max((toptable[[y]]), na.rm = TRUE) + 5), 
                         xlab = bquote(~Log[2] ~ "fold change"), ylab = bquote(~-Log[10] ~ 
                                                                                 italic(P)), axisLabSize = 18, title = "Volcano plot", 
                         subtitle = "EnhancedVolcano", caption = paste0("Total = ", 
                                                                        nrow(toptable), " compounds"), titleLabSize = 18, subtitleLabSize = 14, 
                         captionLabSize = 14, pCutoff = 1e-05, FCcutoff = 1, cutoffLineType = "longdash", 
                         cutoffLineCol = "black", cutoffLineWidth = 0.4, pointSize = 2, 
                         labSize = 3, labCol = "black", labFace = "plain", labhjust = 0.5, 
                         labvjust = 1.5, boxedLabels = FALSE, shape = 19, shapeCustom = NULL, 
                         col = c("grey30", "forestgreen", "royalblue", "red2"), colCustom = NULL, 
                         colAlpha = 1/2, colGradient = NULL, colGradientBreaks = c(pCutoff, 
                                                                                   1), colGradientLabels = c("0", "1.0"), colGradientLimits = c(0, 
                                                                                                                                                1), .legend = c("NS", "Log2 FC", "P", "P & Log2 FC"), 
                         legendLabels = c("NS", expression(Log[2] ~ FC), "BF", 
                                          expression(BF ~ and ~ log[2] ~ FC)), legendPosition = "top", 
                         legendLabSize = 14, legendIconSize = 4, shade = NULL, shadeLabel = NULL, 
                         shadeAlpha = 1/2, shadeFill = "grey", shadeSize = 0.01, shadeBins = 2, 
                         drawConnectors = FALSE, widthConnectors = 0.5, typeConnectors = "closed", 
                         endsConnectors = "first", lengthConnectors = unit(0.01, "npc"), 
                         colConnectors = "grey10", hline = NULL, hlineType = "longdash", 
                         hlineCol = "black", hlineWidth = 0.4, vline = NULL, vlineType = "longdash", 
                         vlineCol = "black", vlineWidth = 0.4, gridlines.major = TRUE, 
                         gridlines.minor = TRUE, border = "partial", borderWidth = 0.8, 
                         borderColour = "black") 
{
  if (!is.numeric(toptable[[x]])) {
    stop(paste(x, " is not numeric!", sep = ""))
  }
  if (!is.numeric(toptable[[y]])) {
    stop(paste(y, " is not numeric!", sep = ""))
  }
  i <- xvals <- yvals <- Sig <- NULL
  toptable <- as.data.frame(toptable)
  toptable$Sig <- "NS"
  toptable$Sig[(abs(toptable[[x]]) > FCcutoff)] <- "FC"
  toptable$Sig[(toptable[[y]] > pCutoff)] <- "P"
  toptable$Sig[(toptable[[y]] > pCutoff) & (abs(toptable[[x]]) > 
                                              FCcutoff)] <- "FC_P"
  toptable$Sig <- factor(toptable$Sig, levels = c("NS", "FC", 
                                                  "P", "FC_P"))
  if (min(toptable[[y]], na.rm = TRUE) == 0) {
    warning(paste("One or more Bayes Factors is 0.", "Converting to 10^-1 * current", 
                  "lowest non-zero Bayes Factor..."), call. = FALSE)
    toptable[which(toptable[[y]] == 0), y] <- min(toptable[which(toptable[[y]] != 
                                                                   0), y], na.rm = TRUE) * 10^-1
  }
  toptable$lab <- lab
  toptable$xvals <- toptable[[x]]
  toptable$yvals <- toptable[[y]]
  if (!is.null(selectLab)) {
    names.new <- rep(NA, length(toptable$lab))
    indices <- which(toptable$lab %in% selectLab)
    names.new[indices] <- toptable$lab[indices]
    toptable$lab <- names.new
  }
  th <- theme_bw(base_size = 24) + theme(legend.background = element_rect(), 
                                         plot.title = element_text(angle = 0, size = titleLabSize, 
                                                                   face = "bold", vjust = 1), plot.subtitle = element_text(angle = 0, 
                                                                                                                           size = subtitleLabSize, face = "plain", vjust = 1), 
                                         plot.caption = element_text(angle = 0, size = captionLabSize, 
                                                                     face = "plain", vjust = 1), axis.text.x = element_text(angle = 0, 
                                                                                                                            size = axisLabSize, vjust = 1), axis.text.y = element_text(angle = 0, 
                                                                                                                                                                                       size = axisLabSize, vjust = 1), axis.title = element_text(size = axisLabSize), 
                                         legend.position = legendPosition, legend.key = element_blank(), 
                                         legend.key.size = unit(0.5, "cm"), legend.text = element_text(size = legendLabSize), 
                                         title = element_text(size = legendLabSize), legend.title = element_blank())
  if (!is.null(colCustom) & !is.null(shapeCustom)) {
    plot <- ggplot(toptable, aes(x = xvals, y = (yvals))) + 
      th + guides(colour = guide_legend(order = 1, override.aes = list(size = legendIconSize)), 
                  shape = guide_legend(order = 2, override.aes = list(size = legendIconSize))) + 
      geom_point(aes(color = factor(names(colCustom)), 
                     shape = factor(names(shapeCustom))), alpha = colAlpha, 
                 size = pointSize, na.rm = TRUE) + scale_color_manual(values = colCustom) + 
      scale_shape_manual(values = shapeCustom)
  }
  else if (!is.null(colCustom) & is.null(shapeCustom) & length(shape) == 
           1) {
    plot <- ggplot(toptable, aes(x = xvals, y = (yvals))) + 
      th + guides(colour = guide_legend(order = 1, override.aes = list(size = legendIconSize)), 
                  shape = guide_legend(order = 2, override.aes = list(size = legendIconSize))) + 
      geom_point(aes(color = factor(names(colCustom))), 
                 alpha = colAlpha, shape = shape, size = pointSize, 
                 na.rm = TRUE) + scale_color_manual(values = colCustom) + 
      scale_shape_manual(guide = TRUE)
  }
  else if (!is.null(colCustom) & is.null(shapeCustom) & length(shape) == 
           4) {
    plot <- ggplot(toptable, aes(x = xvals, y = (yvals))) + 
      th + guides(colour = guide_legend(order = 1, override.aes = list(size = legendIconSize)), 
                  shape = guide_legend(order = 2, override.aes = list(size = legendIconSize))) + 
      geom_point(aes(color = factor(names(colCustom)), 
                     shape = factor(Sig)), alpha = colAlpha, size = pointSize, 
                 na.rm = TRUE) + scale_color_manual(values = colCustom) + 
      scale_shape_manual(values = c(NS = shape[1], FC = shape[2], 
                                    P = shape[3], FC_P = shape[4]), labels = c(NS = legendLabels[1], 
                                                                               FC = legendLabels[2], P = legendLabels[3], FC_P = legendLabels[4]), 
                         guide = TRUE)
  }
  else if (is.null(colCustom) & !is.null(shapeCustom)) {
    if (is.null(colGradient)) {
      plot <- ggplot(toptable, aes(x = xvals, y = (yvals))) + 
        th + guides(colour = guide_legend(order = 1, 
                                          override.aes = list(size = legendIconSize)), 
                    shape = guide_legend(order = 2, override.aes = list(size = legendIconSize))) + 
        geom_point(aes(color = factor(Sig), shape = factor(names(shapeCustom))), 
                   alpha = colAlpha, size = pointSize, na.rm = TRUE) + 
        scale_color_manual(values = c(NS = col[1], FC = col[2], 
                                      P = col[3], FC_P = col[4]), labels = c(NS = legendLabels[1], 
                                                                             FC = legendLabels[2], P = legendLabels[3], 
                                                                             FC_P = legendLabels[4])) + scale_shape_manual(values = shapeCustom)
    }
    else {
      plot <- ggplot(toptable, aes(x = xvals, y = (yvals))) + 
        th + guides(shape = guide_legend(order = 2, override.aes = list(size = legendIconSize))) + 
        geom_point(aes(color = factor(Sig), shape = factor(names(shapeCustom))), 
                   alpha = colAlpha, size = pointSize, na.rm = TRUE) + 
        scale_colour_gradient(low = colGradient[1], high = colGradient[2], 
                              limits = colGradientLimits, breaks = colGradientBreaks, 
                              labels = colGradientLabels)
      scale_shape_manual(values = shapeCustom)
    }
  }
  else if (is.null(colCustom) & is.null(shapeCustom) & length(shape) == 
           1) {
    if (is.null(colGradient)) {
      plot <- ggplot(toptable, aes(x = xvals, y = (yvals))) + 
        th + guides(colour = guide_legend(order = 1, 
                                          override.aes = list(shape = shape, size = legendIconSize))) + 
        geom_point(aes(color = factor(Sig)), alpha = colAlpha, 
                   shape = shape, size = pointSize, na.rm = TRUE) + 
        scale_color_manual(values = c(NS = col[1], FC = col[2], 
                                      P = col[3], FC_P = col[4]), labels = c(NS = legendLabels[1], 
                                                                             FC = legendLabels[2], P = legendLabels[3], 
                                                                             FC_P = legendLabels[4]))
    }
    else {
      plot <- ggplot(toptable, aes(x = xvals, y = (yvals))) + 
        th + geom_point(aes(color = yvals), alpha = colAlpha, 
                        shape = shape, size = pointSize, na.rm = TRUE) + 
        scale_colour_gradient(low = colGradient[1], high = colGradient[2], 
                              limits = colGradientLimits, breaks = colGradientBreaks, 
                              labels = colGradientLabels)
    }
  }
  else if (is.null(colCustom) & is.null(shapeCustom) & length(shape) == 
           4) {
    if (is.null(colGradient)) {
      plot <- ggplot(toptable, aes(x = xvals, y = (yvals))) + 
        th + guides(colour = guide_legend(order = 1, 
                                          override.aes = list(shape = c(NS = shape[1], 
                                                                        FC = shape[2], P = shape[3], FC_P = shape[4]), 
                                                              size = legendIconSize))) + geom_point(aes(color = factor(Sig), 
                                                                                                        shape = factor(Sig)), alpha = colAlpha, size = pointSize, 
                                                                                                    na.rm = TRUE) + scale_color_manual(values = c(NS = col[1], 
                                                                                                                                                  FC = col[2], P = col[3], FC_P = col[4]), labels = c(NS = legendLabels[1], 
                                                                                                                                                                                                      FC = legendLabels[2], P = legendLabels[3], FC_P = legendLabels[4])) + 
        scale_shape_manual(values = c(NS = shape[1], 
                                      FC = shape[2], P = shape[3], FC_P = shape[4]), 
                           guide = FALSE)
    }
    else {
      plot <- ggplot(toptable, aes(x = xvals, y = (yvals))) + 
        th + geom_point(aes(color = yvals, shape = factor(Sig)), 
                        alpha = colAlpha, size = pointSize, na.rm = TRUE) + 
        scale_colour_gradient(low = colGradient[1], high = colGradient[2], 
                              limits = colGradientLimits, breaks = colGradientBreaks, 
                              labels = colGradientLabels) + scale_shape_manual(values = c(NS = shape[1], 
                                                                                          FC = shape[2], P = shape[3], FC_P = shape[4]), 
                                                                               guide = FALSE)
    }
  }
  plot <- plot + xlab(xlab) + ylab(ylab) + xlim(xlim[1], xlim[2]) + 
    ylim(ylim[1], ylim[2]) + geom_vline(xintercept = c(-FCcutoff, 
                                                       FCcutoff), linetype = cutoffLineType, colour = cutoffLineCol, 
                                        size = cutoffLineWidth) + geom_hline(yintercept = (pCutoff), 
                                                                             linetype = cutoffLineType, colour = cutoffLineCol, size = cutoffLineWidth)
  plot <- plot + labs(title = title, subtitle = subtitle, caption = caption)
  if (!is.null(vline)) {
    plot <- plot + geom_vline(xintercept = vline, linetype = vlineType, 
                              colour = vlineCol, size = vlineWidth)
  }
  if (!is.null(hline)) {
    plot <- plot + geom_hline(yintercept = (hline), 
                              linetype = hlineType, colour = hlineCol, size = hlineWidth)
  }
  if (border == "full") {
    plot <- plot + theme(panel.border = element_rect(colour = borderColour, 
                                                     fill = NA, size = borderWidth))
  }
  else if (border == "partial") {
    plot <- plot + theme(axis.line = element_line(size = borderWidth, 
                                                  colour = borderColour), panel.border = element_blank(), 
                         panel.background = element_blank())
  }
  else {
    stop("Unrecognised value passed to 'border'. Must be 'full' or 'partial'")
  }
  if (gridlines.major == TRUE) {
    plot <- plot + theme(panel.grid.major = element_line())
  }
  else {
    plot <- plot + theme(panel.grid.major = element_blank())
  }
  if (gridlines.minor == TRUE) {
    plot <- plot + theme(panel.grid.minor = element_line())
  }
  else {
    plot <- plot + theme(panel.grid.minor = element_blank())
  }
  if (boxedLabels == FALSE) {
    if (drawConnectors == TRUE && is.null(selectLab)) {
      plot <- plot + geom_text_repel(data = subset(toptable, 
                                                   toptable[[y]] > pCutoff & abs(toptable[[x]]) > 
                                                     FCcutoff), aes(label = subset(toptable, toptable[[y]] > 
                                                                                     pCutoff & abs(toptable[[x]]) > FCcutoff)[["lab"]]), 
                                     size = labSize, segment.color = colConnectors, 
                                     segment.size = widthConnectors, arrow = arrow(length = lengthConnectors, 
                                                                                   type = typeConnectors, ends = endsConnectors), 
                                     hjust = labhjust, vjust = labvjust, colour = labCol, 
                                     fontface = labFace, na.rm = TRUE)
    }
    else if (drawConnectors == TRUE && !is.null(selectLab)) {
      plot <- plot + geom_text_repel(data = subset(toptable, 
                                                   !is.na(toptable[["lab"]])), aes(label = subset(toptable, 
                                                                                                  !is.na(toptable[["lab"]]))[["lab"]]), size = labSize, 
                                     segment.color = colConnectors, segment.size = widthConnectors, 
                                     arrow = arrow(length = lengthConnectors, type = typeConnectors, 
                                                   ends = endsConnectors), hjust = labhjust, vjust = labvjust, 
                                     colour = labCol, fontface = labFace, na.rm = TRUE)
    }
    else if (drawConnectors == FALSE && !is.null(selectLab)) {
      plot <- plot + geom_text(data = subset(toptable, 
                                             !is.na(toptable[["lab"]])), aes(label = subset(toptable, 
                                                                                            !is.na(toptable[["lab"]]))[["lab"]]), size = labSize, 
                               check_overlap = TRUE, hjust = labhjust, vjust = labvjust, 
                               colour = labCol, fontface = labFace, na.rm = TRUE)
    }
    else if (drawConnectors == FALSE && is.null(selectLab)) {
      plot <- plot + geom_text(data = subset(toptable, 
                                             toptable[[y]] > pCutoff & abs(toptable[[x]]) > 
                                               FCcutoff), aes(label = subset(toptable, toptable[[y]] > 
                                                                               pCutoff & abs(toptable[[x]]) > FCcutoff)[["lab"]]), 
                               size = labSize, check_overlap = TRUE, hjust = labhjust, 
                               vjust = labvjust, colour = labCol, fontface = labFace, 
                               na.rm = TRUE)
    }
  }
  else {
    if (drawConnectors == TRUE && is.null(selectLab)) {
      plot <- plot + geom_label_repel(data = subset(toptable, 
                                                    toptable[[y]] > pCutoff & abs(toptable[[x]]) > 
                                                      FCcutoff), aes(label = subset(toptable, toptable[[y]] > 
                                                                                      pCutoff & abs(toptable[[x]]) > FCcutoff)[["lab"]]), 
                                      size = labSize, segment.color = colConnectors, 
                                      segment.size = widthConnectors, arrow = arrow(length = lengthConnectors, 
                                                                                    type = typeConnectors, ends = endsConnectors), 
                                      hjust = labhjust, vjust = labvjust, colour = labCol, 
                                      fontface = labFace, na.rm = TRUE)
    }
    else if (drawConnectors == TRUE && !is.null(selectLab)) {
      plot <- plot + geom_label_repel(data = subset(toptable, 
                                                    !is.na(toptable[["lab"]])), aes(label = subset(toptable, 
                                                                                                   !is.na(toptable[["lab"]]))[["lab"]]), size = labSize, 
                                      segment.color = colConnectors, segment.size = widthConnectors, 
                                      arrow = arrow(length = lengthConnectors, type = typeConnectors, 
                                                    ends = endsConnectors), hjust = labhjust, vjust = labvjust, 
                                      colour = labCol, fontface = labFace, na.rm = TRUE)
    }
    else if (drawConnectors == FALSE && !is.null(selectLab)) {
      plot <- plot + geom_label(data = subset(toptable, 
                                              !is.na(toptable[["lab"]])), aes(label = subset(toptable, 
                                                                                             !is.na(toptable[["lab"]]))[["lab"]]), size = labSize, 
                                hjust = labhjust, vjust = labvjust, colour = labCol, 
                                fontface = labFace, na.rm = TRUE)
    }
    else if (drawConnectors == FALSE && is.null(selectLab)) {
      plot <- plot + geom_label(data = subset(toptable, 
                                              toptable[[y]] > pCutoff & abs(toptable[[x]]) > 
                                                FCcutoff), aes(label = subset(toptable, toptable[[y]] > 
                                                                                pCutoff & abs(toptable[[x]]) > FCcutoff)[["lab"]]), 
                                size = labSize, hjust = labhjust, vjust = labvjust, 
                                colour = labCol, fontface = labFace, na.rm = TRUE)
    }
  }
  if (!is.null(shade)) {
    plot <- plot + stat_density2d(data = subset(toptable, 
                                                rownames(toptable) %in% shade), fill = shadeFill, 
                                  alpha = shadeAlpha, geom = "polygon", contour = TRUE, 
                                  size = shadeSize, bins = shadeBins, show.legend = FALSE, 
                                  na.rm = TRUE) + scale_fill_identity(name = shadeLabel, 
                                                                      labels = shadeLabel, guide = "legend")
  }
  return(plot)
}


# we will just look at one example 
#colWT and d172 have some easily viewable differences so they have been selected


#subset dataframes

ColWTT<-subset_of_x$ColWT

d1722<-subset_of_x$d172

j<-1
#empty dataframe to store species specific results
BayesfactorEvidence<-data.frame()
#for loop to go process each metabolite
for (j in 1:43){


View()
#calculate ColWT_d172 log2 fold change difference & combine with bayes factor data 
if(sum(ColWTT[,j])!=0 && sum(d1722[,j])!=0){
  #Bayesian factor analysis
  #this is a shortened version of what we already did above 
  #we could aways reference previous data this is just easier for now 
  ColWTT_d1722_bf<-ttestBF(x=ColWTT[,j],y=d1722[,j], paired = FALSE)
  #pull out evidence it needs to be saved as a vector to be paired with LogFC
  ColWTT_d1722_Evidence<-as.data.frame(ColWTT_d1722_bf)[,1:2]
  rownames(ColWTT_d1722_Evidence)<-c()
  #calculate Log2FC
  #replace 0s with 0.00001
  #in this instance because we may have negative numbers 
  #we need to increase the data so all values are above 0
  #this is unnecessay if all trait values are above 0 
  x1cmpd<-ColWTT[,j]+10
  x1cmpd[x1cmpd==0]<-0.000001
  x2cmpd<-d1722[,j]+10
  x2cmpd[x2cmpd==0]<-0.000001
  FC<- mean(x1cmpd)/mean(x2cmpd)
  ColWTT_d1722log2FC<-log(FC,2)
 #hear we are creating a standardized effect size 
  #standardized to the standard deviation of each trait X genotype 
  #this type of effectsize is called strictly standardized mean difference (SSMD) 
  # SSMD is the average fold change (on the log scale) penalized by the variability of fold change (on the log scale)
  #there are many ways to calculate effect size 
  #cohen'd for large sample size
  ColWTT_d1722_Eff_CohD<-((mean(x1cmpd)-mean(x2cmpd))/sqrt((sd(x1cmpd)^2+ sd(x2cmpd)^2)/2))
  #hedges'g for small sample size * it applies a scalar based on sample size 
  ColWTT_d1722_Eff_HedG<-((mean(x1cmpd)-mean(x2cmpd))/(sqrt( (((length(x1cmpd)-1)*sd(x1cmpd)^2)+(length(x1cmpd)-1)*sd(x2cmpd)^2)/(length(x1cmpd)+length(x2cmpd)-2))  ))
               
  
  #combine evidence and log2FC & effect size 
  ColWTT_d1722_Vectortoadd<-as.data.frame(cbind(ColWTT_d1722_Evidence,ColWTT_d1722log2FC,ColWTT_d1722_Eff_CohD,ColWTT_d1722_Eff_HedG))
  colnames(ColWTT_d1722_Vectortoadd)<-c("ColWTT_d1722_BF","ColWTT_d1722_bferror","ColWTT_d1722log2FC","ColWTT_d1722_Eff","ColWTT_d1722_Eff_HedG")
  BayesfactorEvidence<-rbind(BayesfactorEvidence,ColWTT_d1722_Vectortoadd)
}else{next}}

#View(BayesfactorEvidence)
#for any bayes factor that can't be calculated due to missing data we make it a small number 
BayesfactorEvidence[BayesfactorEvidence==0]<-0.000001

row.names(BayesfactorEvidence)<-colnames(ColWTT[,1:43])

#creates out volcano plot for bayes factor by log2 fold change 
#null is the first which in this case is d172 
#colWT is the second 
ColWTT_d1722_volcanoplot_log2FC<-JordanVolcano(BayesfactorEvidence,
                              lab = rownames(BayesfactorEvidence),
                              title =bquote(italic("A. thaliana") ~ 'd172 vs  ColWT') ,
                              subtitle = "Peak area",
                              FCcutoff = .2,
                              pCutoff = 20,
                              boxedLabels = F,
                              drawConnectors = F,
                              ylim = c(0.001,(max(BayesfactorEvidence$ColWTT_d1722_BF)+1)),
                              xlim = c(-(max(abs(BayesfactorEvidence$ColWTT_d1722log2FC))+1),(max(abs(BayesfactorEvidence$ColWTT_d1722log2FC))+1)),
                              x = 'ColWTT_d1722log2FC',
                              ylab = "Bayes Factor",
                              y = 'ColWTT_d1722_BF')
#view the plot 
ColWTT_d1722_volcanoplot_log2FC

#0.25 is about a 20% relative increase in the respective metabolite

ColWTT_d1722_volcanoplot_EFF<-JordanVolcano(BayesfactorEvidence,
                                        lab = rownames(BayesfactorEvidence),
                                        title =bquote(italic("A. thaliana") ~ 'd172 vs  ColWT') ,
                                        subtitle = "Peak area",
                                        FCcutoff = 2,#x cutoff
                                        pCutoff = 10,#y cutoff
                                        boxedLabels = F,
                                        drawConnectors = F,
                                        ylim = c(0.001,(max(BayesfactorEvidence$ColWTT_d1722_BF)+1)),
                                        xlim = c(-(max(abs(BayesfactorEvidence$ColWTT_d1722_Eff_HedG))+1),(max(abs(BayesfactorEvidence$ColWTT_d1722_Eff_HedG))+1)),
                                        x = 'ColWTT_d1722_Eff_HedG',
                                        ylab = "Bayes Factor",
                                        y = 'ColWTT_d1722_BF',
                                        xlab = paste0("Hedge's g"),
                                        legendLabels = c("NS", "Hedge's g", "BF", 
                                                         paste0("BF and Hedge's g")))
#view Plot 
ColWTT_d1722_volcanoplot_EFF

#effect size value of 1 is equal to 1 pooled standard deviation value for instance
#x18 has an Hedge's g of +8.5565481, 
#the mean of d172 is 8.5565481 pooled standard deviaiton units higher that ColWT
#the pooled original standard deviation *scaled by sample size is 0.2258408
#0.2258408*8.5565481= 1.932418: the mean difference in original units 
# is 87 times more likely that the mean of x18 in Col WT is different that d172,
# further we estimate that the mean of x18 is 1.932418 units higher than  d172,
#however in leveraging MCMC walks we estimate the mean difference to be much lower 
# the error 1.687594 units with a 97% HDPI between (0.9018002-2.220037)


#check out the violin plot to make sure 
print(
  ggviolin(
    y = colnames(x)[23],
    x = "Genotype",
    add = c("jitter", error.plot = "crossbar"),
    data = x,
    draw_quantiles = c(0.025, 0.5, 0.975)
  )
)
