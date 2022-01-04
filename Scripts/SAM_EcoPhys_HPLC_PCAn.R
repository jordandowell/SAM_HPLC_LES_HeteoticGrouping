###SAM Ecophys and HPLC data
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
getwd()
setwd("DATA/")
getwd()
ECOPHYS <- read.csv("SAMecophys_linemeans.csv", row.names = 1)

#change PPN to row names & remove last column

ECOPHYS <- ECOPHYS[, -21]



#Read in Metadata
SAMMETA <- read.csv("SAM_Metadata.csv", row.names = 1)



Combodataframe <- merge(ECOPHYS, SAMMETA, by = 0)

#make column 22-27 factors

Combodataframe[, 22:27] <- lapply(Combodataframe[, 22:27], factor)

st(Combodataframe[, c(3:21, 24)], group = 'BREED')



#create dataframe of summary statistics

#

Combodataframe %>%
  group_by(BREED) %>%
  summarise(mean_Freshmass = mean(FreshMass))




#estimate number of components necessary for imputation
#estimate based on all to narrow search
#PCestimates<- estim_ncpPCA(ECOPHYS[,-1], ncp.min = 0,ncp.max = 20, method.cv = "gcv")

#PCestimates
#min(PCestimates$criterion)


#redo with more iteractions and EM

#PCestimates_2<- estim_ncpPCA(ECOPHYS[,-1], ncp.min = 9,ncp.max = 15,method = "EM",nbsim = 1000, method.cv = "gcv")

#Number of PCs = 12



#impute missing data

#ECOPHYS_Imputation <- imputePCA(ECOPHYS[,-1],ncp=12)

#Save imputed data

#ECOPHYS_Completeset <-ECOPHYS_Imputation$completeObs

#write.csv(ECOPHYS_Completeset,"Ecophys_imputed.csv")
#read in imputed data for assessment
ECOPHYS_Completeset <- read.csv("Ecophys_imputed.csv", row.names = 1)









dim(ECOPHY_Composites)
#create dataset to store Composite Traits

ECOPHY_Composites <- data.frame()
#subset phenolics and flavonoid


FlavPhen_Traits <- ECOPHYS_Completeset[, c(18, 19)]

#   c(3,9,10,12,15,17)]

#PCA of flavenoids and phenolics


FlavPhen.PCA <- prcomp(FlavPhen_Traits, center = T, scale. = T)


#produce eigenvalue plot

pdf("OUTPUT/FlavPhen_eigen.pdf")
fviz_eig(FlavPhen.PCA, choice = "eigenvalue", addlabels = T)
dev.off()


#produce variable loading plots
pdf("OUTPUT/FlavPhen_VariableLoading_axes_1_2.pdf")
fviz_pca_var(
  FlavPhen.PCA,
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
    FlavPhen.PCA,
    geom.ind = "point",
    pointshape = 19,
    label = SAMMETA$CORE12,
    col.ind = SAMMETA$CLASS,
    mean.point = F
  )


ggplot(cbind(PCA.IND$data, SAMMETA),
       aes(x = x,
           y = y,)) +
  geom_point(#color="black",
    aes(
      fill = SAMMETA$CLASS,
      # color="black",
      shape = SAMMETA$CORE12,
      size = SAMMETA$CORE12
    )) +
  scale_size_manual(values = c(5, 2)) +
  scale_shape_manual(values = c(21, 24)) +
  ggtitle(PCA.IND$labels$title) +
  xlab(PCA.IND$labels$x) + ylab(PCA.IND$labels$y) +
  labs(fill = "Breeding group", shape = "Core 12") +
  guides(size = "none", fill = guide_legend(override.aes = list(shape =
                                                                  21))) + geom_hline(yintercept = 0,
                                                                                     linetype = "dashed",
                                                                                     color = "black") + geom_vline(xintercept = 0,
                                                                                                                   linetype = "dashed",
                                                                                                                   color = "black") + theme_bw()
ggplot(cbind(PCA.IND$data, SAMMETA),
       aes(x = x,
           y = y,)) +
  geom_point(#color="black",
    aes(
      fill = SAMMETA$BREED,
      # color="black",
      shape = SAMMETA$CORE12,
      size = SAMMETA$CORE12
    )) +
  scale_size_manual(values = c(5, 2)) +
  scale_shape_manual(values = c(21, 24)) +
  ggtitle(PCA.IND$labels$title) +
  xlab(PCA.IND$labels$x) + ylab(PCA.IND$labels$y) +
  labs(fill = "Breeding group", shape = "Core 12") +
  guides(size = "none", fill = guide_legend(override.aes = list(shape =
                                                                  21))) + geom_hline(yintercept = 0,
                                                                                     linetype = "dashed",
                                                                                     color = "black") + geom_vline(xintercept = 0,
                                                                                                                   linetype = "dashed",
                                                                                                                   color = "black") + theme_bw()
ggplot(cbind(PCA.IND$data, SAMMETA),
       aes(x = x,
           y = y,)) +
  geom_point(#color="black",
    aes(
      fill = SAMMETA$COUNTRY,
      # color="black",
      shape = SAMMETA$CORE12,
      size = SAMMETA$CORE12
    )) +
  scale_size_manual(values = c(5, 2)) +
  scale_shape_manual(values = c(21, 24)) +
  ggtitle(PCA.IND$labels$title) +
  xlab(PCA.IND$labels$x) + ylab(PCA.IND$labels$y) +
  labs(fill = "Breeding group", shape = "Core 12") +
  guides(size = "none", fill = guide_legend(override.aes = list(shape =
                                                                  21))) + geom_hline(yintercept = 0,
                                                                                     linetype = "dashed",
                                                                                     color = "black") + geom_vline(xintercept = 0,
                                                                                                                   linetype = "dashed",
                                                                                                                   color = "black") + theme_bw()
ggplot(cbind(PCA.IND$data, SAMMETA),
       aes(x = x,
           y = y,)) +
  geom_point(#color="black",
    aes(
      fill = SAMMETA$USE,
      # color="black",
      shape = SAMMETA$CORE12,
      size = SAMMETA$CORE12
    )) +
  scale_size_manual(values = c(5, 2)) +
  scale_shape_manual(values = c(21, 24)) +
  ggtitle(PCA.IND$labels$title) +
  xlab(PCA.IND$labels$x) + ylab(PCA.IND$labels$y) +
  labs(fill = "Breeding group", shape = "Core 12") +
  guides(size = "none", fill = guide_legend(override.aes = list(shape =
                                                                  21))) + geom_hline(yintercept = 0,
                                                                                     linetype = "dashed",
                                                                                     color = "black") + geom_vline(xintercept = 0,
                                                                                                                   linetype = "dashed",
                                                                                                                   color = "black") + theme_bw()



dev.off()

#assess cos2

var <- get_pca_var(FlavPhen.PCA)


#create plot for Cos2
pdf("OUTPUT/FlavPhen_cos2.pdf")
corrplot(
  var$cos2,
  method = "color",
  p.mat = round(var$cos2, digits = 2),
  insig = "p-value",
  sig.level = 0.01,
  outline = T,
  cl.lim = c(0, 1),
  is.corr = FALSE
)
dev.off()




#add data to composite data set

ECOPHY_Composites <- FlavPhen.PCA$x



#subset LES traits


LES_Traits <- ECOPHYS_Completeset[, c(1:17)]
LES_Traits <- LES_Traits[, -c(5, 8, 13)]

View(LES_Traits)


#Bayes factor of each group





#count()




#continue
Combodataframe$BREED <- as.factor(Combodataframe$BREED)

breeds <- c("HA", "RHA")

Combodataframe <- Combodataframe[which(Combodataframe$BREED == breeds), ]

#start for loop for BF and violin plots


pdf("OUTPUT/Ecophys_violinplots_breed.pdf")
for (i in 2:20) {
  print(
    ggviolin(
      y = colnames(Combodataframe)[i],
      x = "BREED",
      add = c("jitter", error.plot = "crossbar"),
      data = Combodataframe,
      draw_quantiles = c(0.025, 0.5, 0.975)
    )
  )
}
dev.off()


#bayes factor analysis

#create a for loop for results of LES by breed

View(Combodataframe)
ECOPHYS_results <- data.frame()
#colnames(ECOPHYS_results)<-c("BF","BFerror","mean","mean.SD","quant.2.5","quant.97.5")

for (i in 2:20) {
  trait.range <-
    range(Combodataframe[, i])[2] - range(Combodataframe[, i])[1]
  
  Trait <- colnames(Combodataframe)[i]
  #bftest has a error but we need the loop to continue
  bf <-
    ttestBF(formula = reformulate(
      termlabels = "BREED",
      response = colnames(Combodataframe)[i]
    ),
    data = Combodataframe)
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
  row.names(results) <- Trait
  colnames(results) <-
    c("BF", "BFerror", "mean", "mean.SD", "quant.2.5", "quant.97.5")
  ECOPHYS_results <- rbind(ECOPHYS_results, results)
}
#add them to large results file

View(ECOPHYS_results)

write.csv(ECOPHYS_results, "OUTPUT/Ecophys_violinplots_breed.csv")



#Bayes factor of each group

Combodataframe <- merge(ECOPHYS_Completeset, SAMMETA, by = 0)

Combodataframe$USE <- as.factor(Combodataframe$USE)


breeds <- c("Oil", "NonOil")

Combodataframe <- Combodataframe[which(Combodataframe$USE == breeds), ]

#start for loop for BF and violin plots


pdf("OUTPUT/Ecophys_violinplots_USE.pdf")
for (i in 2:20) {
  print(ggviolin(
    y = colnames(Combodataframe)[i],
    x = "USE",
    add = c("jitter", error.plot = "crossbar"),
    data = Combodataframe,
    draw_quantiles = c(0.025, 0.5, 0.975)
  ))
}
dev.off()


#bayes factor analysis

#create a for loop for results of LES by breed


ECOPHYS_results <- data.frame()
#colnames(ECOPHYS_results)<-c("BF","BFerror","mean","mean.SD","quant.2.5","quant.97.5")

for (i in 2:20) {
  trait.range <-
    range(Combodataframe[, i])[2] - range(Combodataframe[, i])[1]
  
  Trait <- colnames(Combodataframe)[i]
  #bftest has a error but we need the loop to continue
  bf <-
    ttestBF(formula = reformulate(
      termlabels = "USE",
      response = colnames(Combodataframe)[i]
    ),
    data = Combodataframe)
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
  row.names(results) <- Trait
  colnames(results) <-
    c("BF", "BFerror", "mean", "mean.SD", "quant.2.5", "quant.97.5")
  ECOPHYS_results <- rbind(ECOPHYS_results, results)
}
#add them to large results file

View(ECOPHYS_results)

write.csv(ECOPHYS_results, "OUTPUT/Ecophys_violinplots_USE.csv")







#PCA of LES TRaits


LES_Traits.PCA <- prcomp(LES_Traits, center = T, scale. = T)
PCA.IND <-
  fviz_pca_ind(
    LES_Traits.PCA,
    geom = "point",
    axes = c(1, 2),
    label = "none",
    habillage = SAMMETA$CLASS,
    mean.point = F
  )

#Bayes factor of each group

Combodataframe <- cbind(PCA.IND$data, SAMMETA)

Combodataframe$BREED <- as.factor(Combodataframe$BREED)

breeds <- c("HA", "RHA")

Combodataframe <- Combodataframe[which(Combodataframe$BREED == breeds), ]
#change between PC1 & PC2
View(Combodataframe)
#trait.range<-range(Combodataframe$PC1)[2]-range(Combodataframe$PC1)[1]

pdf("OUTPUT/LES_violinplots_breed.pdf")
ggviolin(
  y = "x",
  ylab = "PC1",
  x = "BREED",
  add = c("jitter", error.plot = "crossbar"),
  data = Combodataframe,
  draw_quantiles = c(0.025, 0.5, 0.975)
)
ggviolin(
  y = "y",
  ylab = "PC2",
  x = "BREED",
  add = c("jitter", error.plot = "crossbar"),
  data = Combodataframe,
  draw_quantiles = c(0.025, 0.5, 0.975)
)
dev.off()


#bayes factor analysis

#create a for loop for results of LES by breed


ECOPHYS_results <- data.frame()
#colnames(ECOPHYS_results)<-c("BF","BFerror","mean","mean.SD","quant.2.5","quant.97.5")

for (i in 3:4) {
  trait.range <-
    range(Combodataframe[, i])[2] - range(Combodataframe[, i])[1]
  
  Trait <- colnames(Combodataframe)[i]
  #bftest has a error but we need the loop to continue
  bf <-
    ttestBF(formula = reformulate(
      termlabels = "BREED",
      response = colnames(Combodataframe)[i]
    ),
    data = Combodataframe)
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
  row.names(results) <- Trait
  colnames(results) <-
    c("BF", "BFerror", "mean", "mean.SD", "quant.2.5", "quant.97.5")
  ECOPHYS_results <- rbind(ECOPHYS_results, results)
}
#add them to large results file

View(ECOPHYS_results)

write.csv(ECOPHYS_results, "OUTPUT/LES_bayesFactor_BREED.csv")



#violin plots and bf of oil vs nonoil

Combodataframe <- cbind(PCA.IND$data, SAMMETA)

Combodataframe$USE <- as.factor(Combodataframe$USE)

USE <- c("Oil", "NonOil")

Combodataframe <- Combodataframe[which(Combodataframe$USE == USE), ]


unique(Combodataframe$USE)

#violin plots and bf of oil vs nonoil

pdf("OUTPUT/LES_violinplots_use.pdf")
ggviolin(
  y = "x",
  ylab = "PC1",
  x = "USE",
  add = c("jitter", error.plot = "crossbar"),
  data = Combodataframe,
  draw_quantiles = c(0.025, 0.5, 0.975)
)
ggviolin(
  y = "y",
  ylab = "PC2",
  x = "USE",
  add = c("jitter", error.plot = "crossbar"),
  data = Combodataframe,
  draw_quantiles = c(0.025, 0.5, 0.975)
)
dev.off()
#create a for loop for results of LES by breed


ECOPHYS_results <- data.frame()
#colnames(ECOPHYS_results)<-c("BF","BFerror","mean","mean.SD","quant.2.5","quant.97.5")

for (i in 3:4) {
  trait.range <-
    range(Combodataframe[, i])[2] - range(Combodataframe[, i])[1]
  
  Trait <- colnames(Combodataframe)[i]
  #bftest has a error but we need the loop to continue
  bf <-
    ttestBF(formula = reformulate(
      termlabels = "USE",
      response = colnames(Combodataframe)[i]
    ),
    data = Combodataframe)
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
  row.names(results) <- Trait
  colnames(results) <-
    c("BF", "BFerror", "mean", "mean.SD", "quant.2.5", "quant.97.5")
  ECOPHYS_results <- rbind(ECOPHYS_results, results)
}
#add them to large results file

View(ECOPHYS_results)

write.csv(ECOPHYS_results, "OUTPUT/LES_bayesFactor_USE.csv")



#produce eigenvalue plot

pdf("OUTPUT/LES_eigen.pdf")
fviz_eig(LES_Traits.PCA, choice = "eigenvalue", addlabels = T)
dev.off()


#produce variable loading plots
pdf("OUTPUT/LES_VariableLoading_axes_1_2.pdf")
fviz_pca_var(
  LES_Traits.PCA,
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
pdf("OUTPUT/LES_IndividualLoading_axes_1_2.pdf")
PCA.IND <-
  fviz_pca_ind(
    LES_Traits.PCA,
    geom = "point",
    axes = c(1, 2),
    label = "none",
    habillage = SAMMETA$CLASS,
    mean.point = F
  )

ggplot(cbind(PCA.IND$data, SAMMETA),
       aes(x = x,
           y = y,)) +
  geom_point(#color="black",
    aes(
      fill = SAMMETA$CLASS,
      # color="black",
      shape = SAMMETA$CORE12,
      size = SAMMETA$CORE12
    )) +
  scale_size_manual(values = c(5, 2)) +
  scale_shape_manual(values = c(21, 24)) +
  ggtitle(PCA.IND$labels$title) +
  xlab(PCA.IND$labels$x) + ylab(PCA.IND$labels$y) +
  labs(fill = "Breeding group", shape = "Core 12") +
  guides(size = "none", fill = guide_legend(override.aes = list(shape =
                                                                  21))) + geom_hline(yintercept = 0,
                                                                                     linetype = "dashed",
                                                                                     color = "black") + geom_vline(xintercept = 0,
                                                                                                                   linetype = "dashed",
                                                                                                                   color = "black") + theme_bw()
ggplot(cbind(PCA.IND$data, SAMMETA),
       aes(x = x,
           y = y,)) +
  geom_point(#color="black",
    aes(
      fill = SAMMETA$BREED,
      # color="black",
      shape = SAMMETA$CORE12,
      size = SAMMETA$CORE12
    )) +
  scale_size_manual(values = c(5, 2)) +
  scale_shape_manual(values = c(21, 24)) +
  ggtitle(PCA.IND$labels$title) +
  xlab(PCA.IND$labels$x) + ylab(PCA.IND$labels$y) +
  labs(fill = "Breeding group", shape = "Core 12") +
  guides(size = "none", fill = guide_legend(override.aes = list(shape =
                                                                  21))) + geom_hline(yintercept = 0,
                                                                                     linetype = "dashed",
                                                                                     color = "black") + geom_vline(xintercept = 0,
                                                                                                                   linetype = "dashed",
                                                                                                                   color = "black") + theme_bw()
ggplot(cbind(PCA.IND$data, SAMMETA),
       aes(x = x,
           y = y,)) +
  geom_point(#color="black",
    aes(
      fill = SAMMETA$COUNTRY,
      # color="black",
      shape = SAMMETA$CORE12,
      size = SAMMETA$CORE12
    )) +
  scale_size_manual(values = c(5, 2)) +
  scale_shape_manual(values = c(21, 24)) +
  ggtitle(PCA.IND$labels$title) +
  xlab(PCA.IND$labels$x) + ylab(PCA.IND$labels$y) +
  labs(fill = "Breeding group", shape = "Core 12") +
  guides(size = "none", fill = guide_legend(override.aes = list(shape =
                                                                  21))) + geom_hline(yintercept = 0,
                                                                                     linetype = "dashed",
                                                                                     color = "black") + geom_vline(xintercept = 0,
                                                                                                                   linetype = "dashed",
                                                                                                                   color = "black") + theme_bw()
ggplot(cbind(PCA.IND$data, SAMMETA),
       aes(x = x,
           y = y,)) +
  geom_point(#color="black",
    aes(
      fill = SAMMETA$USE,
      # color="black",
      shape = SAMMETA$CORE12,
      size = SAMMETA$CORE12
    )) +
  scale_size_manual(values = c(5, 2)) +
  scale_shape_manual(values = c(21, 24)) +
  ggtitle(PCA.IND$labels$title) +
  xlab(PCA.IND$labels$x) + ylab(PCA.IND$labels$y) +
  labs(fill = "Breeding group", shape = "Core 12") +
  guides(size = "none", fill = guide_legend(override.aes = list(shape =
                                                                  21))) + geom_hline(yintercept = 0,
                                                                                     linetype = "dashed",
                                                                                     color = "black") + geom_vline(xintercept = 0,
                                                                                                                   linetype = "dashed",
                                                                                                                   color = "black") + theme_bw()




dev.off()

#assess cos2

var <- get_pca_var(LES_Traits.PCA)


#create plot for Cos2
pdf("OUTPUT/LES_cos2.pdf")
corrplot(
  var$cos2[, 1:10],
  method = "color",
  p.mat = round(var$cos2[, 1:10], digits = 2),
  insig = "p-value",
  sig.level = (1 / nrow(var$cos2)),
  outline = T,
  cl.lim = c(0, 1),
  is.corr = FALSE
)
dev.off()



#add data to composite data set

ECOPHY_Composites <- cbind(FlavPhen.PCA$x, LES_Traits.PCA$x[, 1:10])




colnames(ECOPHY_Composites) <-
  c(
    "FLAVPHEN_PC1",
    "FLAVPHEN_PC2",
    "LES_PC1",
    "LES_PC2",
    "LES_PC3",
    "LES_PC4",
    "LES_PC5",
    "LES_PC6",
    "LES_PC7",
    "LES_PC8",
    "LES_PC9",
    "LES_PC10"
  )

#save ecophys composite traits
write.csv(ECOPHY_Composites, "Ecophys_CompositeTraits.csv")

#create a summary stats table of data
Combodataframe <- merge(SAMMETA, ECOPHYS, by = 0)
rownames(Combodataframe) <- Combodataframe[, 1]
Combodataframe <- merge(Combodataframe[, -1], ECOPHY_Composites, by = 0)



#make column 1:8 factors

Combodataframe[, 1:8] <- lapply(Combodataframe[, 1:8], factor)
#write summary statistics to CSV


Ecophysummary <-
  describeBy(Combodataframe[, c(4, 9:ncol(Combodataframe))], group = 'BREED', mat =
               TRUE)

write.csv(Ecophysummary, "OUTPUT/EcophysSummarySTATS_BREED.csv")

Ecophysummary <-
  describeBy(Combodataframe[, c(7, 9:ncol(Combodataframe))], group = 'USE', mat =
               TRUE)

write.csv(Ecophysummary, "OUTPUT/EcophysSummarySTATS_USE.csv")





###Sam HPLC PCA


#IMPORT ALL DATA SETS
PeakAreas <-
  read.csv("SAMHPLC_PLOTPEAKAREAS_ROUNDED.csv",
           header = T,
           row.names = 1)
PeakAreas200 <-
  read.csv(
    "SAMHPLC_Retentiontimes200_PlotPeakAreas.csv",
    header = T,
    row.names = 1
  )
RelativeRatios <-
  read.csv("SAMHPLC_RelativeRatios.csv",
           header = T,
           row.names = 1)
RelativeRatios200 <-
  read.csv(
    "SAMHPLC_Retentiontimes200_RelativeRatios.csv",
    header = T,
    row.names = 1
  )

#create Composite dataset to store 10axes each

HPLC_Composites <- data.frame()

#change all NA to 0

PeakAreas[is.na(PeakAreas)] <- 0
PeakAreas200[is.na(PeakAreas200)] <- 0
RelativeRatios[is.na(RelativeRatios)] <- 0
RelativeRatios200[is.na(RelativeRatios200)] <- 0

#change column names to incluse dataframe
colnames(PeakAreas) <- paste0("PeakAreas_RT", colnames(PeakAreas))
colnames(PeakAreas200) <-
  paste0("PeakAreas_RT", colnames(PeakAreas200))
colnames(RelativeRatios) <-
  paste0("PeakAreas_RT", colnames(RelativeRatios))
colnames(RelativeRatios200) <-
  paste0("PeakAreas_RT", colnames(RelativeRatios200))
#cbind because rows are same order
Combodataframe <-
  cbind.data.frame(PeakAreas, PeakAreas200, RelativeRatios, RelativeRatios200)
View(Combodataframe)
#create a summary stats table of data
Combodataframe <- merge(SAMMETA, Combodataframe, by = 0)



#make column 1:7 factors

Combodataframe[, 1:7] <- lapply(Combodataframe[, 1:7], factor)
#write summary statistics to CSV


Ecophysummary <-
  describeBy(Combodataframe[, c(4, 9:ncol(Combodataframe))], group = 'BREED', mat =
               TRUE)

write.csv(Ecophysummary, "OUTPUT/HPLC_PEAKS_SummarySTATS_BREED.csv")

Ecophysummary <-
  describeBy(Combodataframe[, c(7, 9:ncol(Combodataframe))], group = 'USE', mat =
               TRUE)

write.csv(Ecophysummary, "OUTPUT/HPLC_PEAKS_STATS_USE.csv")






#compute PCA

PeakAreas.PCA <- prcomp(PeakAreas, center = T, scale. = T)
PCA.IND <-
  fviz_pca_ind(
    PeakAreas.PCA,
    geom = "point",
    axes = c(1, 2),
    label = "none",
    habillage = SAMMETA$CLASS,
    mean.point = F
  )

#produce eigenvalue plot

pdf("OUTPUT/PeakAreas_eigen.pdf")
fviz_eig(PeakAreas.PCA, choice = "eigenvalue", addlabels = T)
dev.off()

#Bayes factor of each group

Combodataframe <- cbind(PCA.IND$data, SAMMETA)

Combodataframe$BREED <- as.factor(Combodataframe$BREED)

breeds <- c("HA", "RHA")

Combodataframe <- Combodataframe[which(Combodataframe$BREED == breeds), ]
#change between PC1 & PC2
View(Combodataframe)
#trait.range<-range(Combodataframe$PC1)[2]-range(Combodataframe$PC1)[1]

pdf("OUTPUT/PeakAreas_violinplots_breed.pdf")
ggviolin(
  y = "x",
  ylab = "PC1",
  x = "BREED",
  add = c("jitter", error.plot = "crossbar"),
  data = Combodataframe,
  draw_quantiles = c(0.025, 0.5, 0.975)
)
ggviolin(
  y = "y",
  ylab = "PC2",
  x = "BREED",
  add = c("jitter", error.plot = "crossbar"),
  data = Combodataframe,
  draw_quantiles = c(0.025, 0.5, 0.975)
)
dev.off()


#bayes factor analysis

#create a for loop for results of LES by breed


ECOPHYS_results <- data.frame()
#colnames(ECOPHYS_results)<-c("BF","BFerror","mean","mean.SD","quant.2.5","quant.97.5")

for (i in 3:4) {
  trait.range <-
    range(Combodataframe[, i])[2] - range(Combodataframe[, i])[1]
  
  Trait <- colnames(Combodataframe)[i]
  #bftest has a error but we need the loop to continue
  bf <-
    ttestBF(formula = reformulate(
      termlabels = "BREED",
      response = colnames(Combodataframe)[i]
    ),
    data = Combodataframe)
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
  row.names(results) <- Trait
  colnames(results) <-
    c("BF", "BFerror", "mean", "mean.SD", "quant.2.5", "quant.97.5")
  ECOPHYS_results <- rbind(ECOPHYS_results, results)
}
#add them to large results file

View(ECOPHYS_results)

write.csv(ECOPHYS_results, "OUTPUT/PeakAreas_bayesFactor_BREED.csv")



#violin plots and bf of oil vs nonoil

Combodataframe <- cbind(PCA.IND$data, SAMMETA)

Combodataframe$USE <- as.factor(Combodataframe$USE)

USE <- c("Oil", "NonOil")

Combodataframe <- Combodataframe[which(Combodataframe$USE == USE), ]


unique(Combodataframe$USE)

#violin plots and bf of oil vs nonoil

pdf("OUTPUT/PeakAreas_violinplots_use.pdf")
ggviolin(
  y = "x",
  ylab = "PC1",
  x = "USE",
  add = c("jitter", error.plot = "crossbar"),
  data = Combodataframe,
  draw_quantiles = c(0.025, 0.5, 0.975)
)
ggviolin(
  y = "y",
  ylab = "PC2",
  x = "USE",
  add = c("jitter", error.plot = "crossbar"),
  data = Combodataframe,
  draw_quantiles = c(0.025, 0.5, 0.975)
)
dev.off()
#create a for loop for results of LES by breed


ECOPHYS_results <- data.frame()
#colnames(ECOPHYS_results)<-c("BF","BFerror","mean","mean.SD","quant.2.5","quant.97.5")

for (i in 3:4) {
  trait.range <-
    range(Combodataframe[, i])[2] - range(Combodataframe[, i])[1]
  
  Trait <- colnames(Combodataframe)[i]
  #bftest has a error but we need the loop to continue
  bf <-
    ttestBF(formula = reformulate(
      termlabels = "USE",
      response = colnames(Combodataframe)[i]
    ),
    data = Combodataframe)
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
  row.names(results) <- Trait
  colnames(results) <-
    c("BF", "BFerror", "mean", "mean.SD", "quant.2.5", "quant.97.5")
  ECOPHYS_results <- rbind(ECOPHYS_results, results)
}
#add them to large results file

View(ECOPHYS_results)

write.csv(ECOPHYS_results, "OUTPUT/PeakAreas_bayesFactor_USE.csv")



#produce plot of individuals labels by core 12
pdf("OUTPUT/PeakAreas_IndividualLoading_axes_1_2.pdf")
PCA.IND <-
  fviz_pca_ind(
    PeakAreas.PCA,
    geom = "point",
    axes = c(1, 2),
    label = "none",
    habillage = SAMMETA$CLASS,
    mean.point = F
  )


ggplot(cbind(PCA.IND$data, SAMMETA),
       aes(x = x,
           y = y,)) +
  geom_point(#color="black",
    aes(
      fill = SAMMETA$CLASS,
      # color="black",
      shape = SAMMETA$CORE12,
      size = SAMMETA$CORE12
    )) +
  scale_size_manual(values = c(5, 2)) +
  scale_shape_manual(values = c(21, 24)) +
  ggtitle(PCA.IND$labels$title) +
  xlab(PCA.IND$labels$x) + ylab(PCA.IND$labels$y) +
  labs(fill = "Breeding group", shape = "Core 12") +
  guides(size = "none", fill = guide_legend(override.aes = list(shape =
                                                                  21))) + geom_hline(yintercept = 0,
                                                                                     linetype = "dashed",
                                                                                     color = "black") + geom_vline(xintercept = 0,
                                                                                                                   linetype = "dashed",
                                                                                                                   color = "black") + theme_bw()
ggplot(cbind(PCA.IND$data, SAMMETA),
       aes(x = x,
           y = y,)) +
  geom_point(#color="black",
    aes(
      fill = SAMMETA$BREED,
      # color="black",
      shape = SAMMETA$CORE12,
      size = SAMMETA$CORE12
    )) +
  scale_size_manual(values = c(5, 2)) +
  scale_shape_manual(values = c(21, 24)) +
  ggtitle(PCA.IND$labels$title) +
  xlab(PCA.IND$labels$x) + ylab(PCA.IND$labels$y) +
  labs(fill = "Breeding group", shape = "Core 12") +
  guides(size = "none", fill = guide_legend(override.aes = list(shape =
                                                                  21))) + geom_hline(yintercept = 0,
                                                                                     linetype = "dashed",
                                                                                     color = "black") + geom_vline(xintercept = 0,
                                                                                                                   linetype = "dashed",
                                                                                                                   color = "black") + theme_bw()
ggplot(cbind(PCA.IND$data, SAMMETA),
       aes(x = x,
           y = y,)) +
  geom_point(#color="black",
    aes(
      fill = SAMMETA$COUNTRY,
      # color="black",
      shape = SAMMETA$CORE12,
      size = SAMMETA$CORE12
    )) +
  scale_size_manual(values = c(5, 2)) +
  scale_shape_manual(values = c(21, 24)) +
  ggtitle(PCA.IND$labels$title) +
  xlab(PCA.IND$labels$x) + ylab(PCA.IND$labels$y) +
  labs(fill = "Breeding group", shape = "Core 12") +
  guides(size = "none", fill = guide_legend(override.aes = list(shape =
                                                                  21))) + geom_hline(yintercept = 0,
                                                                                     linetype = "dashed",
                                                                                     color = "black") + geom_vline(xintercept = 0,
                                                                                                                   linetype = "dashed",
                                                                                                                   color = "black") + theme_bw()
ggplot(cbind(PCA.IND$data, SAMMETA),
       aes(x = x,
           y = y,)) +
  geom_point(#color="black",
    aes(
      fill = SAMMETA$USE,
      # color="black",
      shape = SAMMETA$CORE12,
      size = SAMMETA$CORE12
    )) +
  scale_size_manual(values = c(5, 2)) +
  scale_shape_manual(values = c(21, 24)) +
  ggtitle(PCA.IND$labels$title) +
  xlab(PCA.IND$labels$x) + ylab(PCA.IND$labels$y) +
  labs(fill = "Breeding group", shape = "Core 12") +
  guides(size = "none", fill = guide_legend(override.aes = list(shape =
                                                                  21))) + geom_hline(yintercept = 0,
                                                                                     linetype = "dashed",
                                                                                     color = "black") + geom_vline(xintercept = 0,
                                                                                                                   linetype = "dashed",
                                                                                                                   color = "black") + theme_bw()

dev.off()
#addto HPLC composite Traits

HPLC_Composites <- PeakAreas.PCA$x[, 1:10]

###PEAK AREAS >200

#compute PCA

PeakAreas200.PCA <- prcomp(PeakAreas200, center = T, scale. = T)
PCA.IND <-
  fviz_pca_ind(
    PeakAreas200.PCA,
    geom = "point",
    axes = c(1, 2),
    label = "none",
    habillage = SAMMETA$CLASS,
    mean.point = F
  )

#Bayes factor of each group

Combodataframe <- cbind(PCA.IND$data, SAMMETA)

Combodataframe$BREED <- as.factor(Combodataframe$BREED)

breeds <- c("HA", "RHA")

Combodataframe <- Combodataframe[which(Combodataframe$BREED == breeds), ]
#change between PC1 & PC2
View(Combodataframe)
#trait.range<-range(Combodataframe$PC1)[2]-range(Combodataframe$PC1)[1]

pdf("OUTPUT/PeakAreas200_violinplots_breed.pdf")
ggviolin(
  y = "x",
  ylab = "PC1",
  x = "BREED",
  add = c("jitter", error.plot = "crossbar"),
  data = Combodataframe,
  draw_quantiles = c(0.025, 0.5, 0.975)
)
ggviolin(
  y = "y",
  ylab = "PC2",
  x = "BREED",
  add = c("jitter", error.plot = "crossbar"),
  data = Combodataframe,
  draw_quantiles = c(0.025, 0.5, 0.975)
)
dev.off()


#bayes factor analysis

#create a for loop for results of LES by breed


ECOPHYS_results <- data.frame()
#colnames(ECOPHYS_results)<-c("BF","BFerror","mean","mean.SD","quant.2.5","quant.97.5")

for (i in 3:4) {
  trait.range <-
    range(Combodataframe[, i])[2] - range(Combodataframe[, i])[1]
  
  Trait <- colnames(Combodataframe)[i]
  #bftest has a error but we need the loop to continue
  bf <-
    ttestBF(formula = reformulate(
      termlabels = "BREED",
      response = colnames(Combodataframe)[i]
    ),
    data = Combodataframe)
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
  row.names(results) <- Trait
  colnames(results) <-
    c("BF", "BFerror", "mean", "mean.SD", "quant.2.5", "quant.97.5")
  ECOPHYS_results <- rbind(ECOPHYS_results, results)
}
#add them to large results file

View(ECOPHYS_results)

write.csv(ECOPHYS_results,
          "OUTPUT/PeakAreas200_bayesFactor_BREED.csv")



#violin plots and bf of oil vs nonoil

Combodataframe <- cbind(PCA.IND$data, SAMMETA)

Combodataframe$USE <- as.factor(Combodataframe$USE)

USE <- c("Oil", "NonOil")

Combodataframe <- Combodataframe[which(Combodataframe$USE == USE), ]


unique(Combodataframe$USE)

#violin plots and bf of oil vs nonoil

pdf("OUTPUT/PeakAreas200_violinplots_use.pdf")
ggviolin(
  y = "x",
  ylab = "PC1",
  x = "USE",
  add = c("jitter", error.plot = "crossbar"),
  data = Combodataframe,
  draw_quantiles = c(0.025, 0.5, 0.975)
)
ggviolin(
  y = "y",
  ylab = "PC2",
  x = "USE",
  add = c("jitter", error.plot = "crossbar"),
  data = Combodataframe,
  draw_quantiles = c(0.025, 0.5, 0.975)
)
dev.off()
#create a for loop for results of LES by breed


ECOPHYS_results <- data.frame()
#colnames(ECOPHYS_results)<-c("BF","BFerror","mean","mean.SD","quant.2.5","quant.97.5")

for (i in 3:4) {
  trait.range <-
    range(Combodataframe[, i])[2] - range(Combodataframe[, i])[1]
  
  Trait <- colnames(Combodataframe)[i]
  #bftest has a error but we need the loop to continue
  bf <-
    ttestBF(formula = reformulate(
      termlabels = "USE",
      response = colnames(Combodataframe)[i]
    ),
    data = Combodataframe)
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
  row.names(results) <- Trait
  colnames(results) <-
    c("BF", "BFerror", "mean", "mean.SD", "quant.2.5", "quant.97.5")
  ECOPHYS_results <- rbind(ECOPHYS_results, results)
}
#add them to large results file

View(ECOPHYS_results)

write.csv(ECOPHYS_results, "OUTPUT/PeakAreas200_bayesFactor_USE.csv")




#produce eigenvalue plot

pdf("OUTPUT/PeakAreas200_eigen.pdf")
fviz_eig(PeakAreas200.PCA, choice = "eigenvalue", addlabels = T)
dev.off()


#produce plot of individuals labels by core 12
pdf("OUTPUT/PeakAreas200_IndividualLoading_axes_1_2.pdf")
PCA.IND <-
  fviz_pca_ind(
    PeakAreas200.PCA,
    geom = "point",
    axes = c(1, 2),
    label = "none",
    habillage = SAMMETA$CLASS,
    mean.point = F
  )
ggplot(cbind(PCA.IND$data, SAMMETA),
       aes(x = x,
           y = y,)) +
  geom_point(#color="black",
    aes(
      fill = SAMMETA$CLASS,
      # color="black",
      shape = SAMMETA$CORE12,
      size = SAMMETA$CORE12
    )) +
  scale_size_manual(values = c(5, 2)) +
  scale_shape_manual(values = c(21, 24)) +
  ggtitle(PCA.IND$labels$title) +
  xlab(PCA.IND$labels$x) + ylab(PCA.IND$labels$y) +
  labs(fill = "Breeding group", shape = "Core 12") +
  guides(size = "none", fill = guide_legend(override.aes = list(shape =
                                                                  21))) + geom_hline(yintercept = 0,
                                                                                     linetype = "dashed",
                                                                                     color = "black") + geom_vline(xintercept = 0,
                                                                                                                   linetype = "dashed",
                                                                                                                   color = "black") + theme_bw()
ggplot(cbind(PCA.IND$data, SAMMETA),
       aes(x = x,
           y = y,)) +
  geom_point(#color="black",
    aes(
      fill = SAMMETA$BREED,
      # color="black",
      shape = SAMMETA$CORE12,
      size = SAMMETA$CORE12
    )) +
  scale_size_manual(values = c(5, 2)) +
  scale_shape_manual(values = c(21, 24)) +
  ggtitle(PCA.IND$labels$title) +
  xlab(PCA.IND$labels$x) + ylab(PCA.IND$labels$y) +
  labs(fill = "Breeding group", shape = "Core 12") +
  guides(size = "none", fill = guide_legend(override.aes = list(shape =
                                                                  21))) + geom_hline(yintercept = 0,
                                                                                     linetype = "dashed",
                                                                                     color = "black") + geom_vline(xintercept = 0,
                                                                                                                   linetype = "dashed",
                                                                                                                   color = "black") + theme_bw()
ggplot(cbind(PCA.IND$data, SAMMETA),
       aes(x = x,
           y = y,)) +
  geom_point(#color="black",
    aes(
      fill = SAMMETA$COUNTRY,
      # color="black",
      shape = SAMMETA$CORE12,
      size = SAMMETA$CORE12
    )) +
  scale_size_manual(values = c(5, 2)) +
  scale_shape_manual(values = c(21, 24)) +
  ggtitle(PCA.IND$labels$title) +
  xlab(PCA.IND$labels$x) + ylab(PCA.IND$labels$y) +
  labs(fill = "Breeding group", shape = "Core 12") +
  guides(size = "none", fill = guide_legend(override.aes = list(shape =
                                                                  21))) + geom_hline(yintercept = 0,
                                                                                     linetype = "dashed",
                                                                                     color = "black") + geom_vline(xintercept = 0,
                                                                                                                   linetype = "dashed",
                                                                                                                   color = "black") + theme_bw()
ggplot(cbind(PCA.IND$data, SAMMETA),
       aes(x = x,
           y = y,)) +
  geom_point(#color="black",
    aes(
      fill = SAMMETA$USE,
      # color="black",
      shape = SAMMETA$CORE12,
      size = SAMMETA$CORE12
    )) +
  scale_size_manual(values = c(5, 2)) +
  scale_shape_manual(values = c(21, 24)) +
  ggtitle(PCA.IND$labels$title) +
  xlab(PCA.IND$labels$x) + ylab(PCA.IND$labels$y) +
  labs(fill = "Breeding group", shape = "Core 12") +
  guides(size = "none", fill = guide_legend(override.aes = list(shape =
                                                                  21))) + geom_hline(yintercept = 0,
                                                                                     linetype = "dashed",
                                                                                     color = "black") + geom_vline(xintercept = 0,
                                                                                                                   linetype = "dashed",
                                                                                                                   color = "black") + theme_bw()

dev.off()


#addto HPLC composite Traits

HPLC_Composites <- cbind(HPLC_Composites, PeakAreas200.PCA$x[, 1:10])


###RelativeRatios > 200

#compute PCA

RelativeRatios.PCA <- prcomp(RelativeRatios, center = T, scale. = T)
PCA.IND <-
  fviz_pca_ind(
    RelativeRatios.PCA,
    geom = "point",
    axes = c(1, 2),
    label = "none",
    habillage = SAMMETA$CLASS,
    mean.point = F
  )

#Bayes factor of each group

Combodataframe <- cbind(PCA.IND$data, SAMMETA)

Combodataframe$BREED <- as.factor(Combodataframe$BREED)

breeds <- c("HA", "RHA")

Combodataframe <- Combodataframe[which(Combodataframe$BREED == breeds), ]
#change between PC1 & PC2
View(Combodataframe)
#trait.range<-range(Combodataframe$x)[2]-range(Combodataframe$y)[1]

pdf("OUTPUT/RelativeRatios_violinplots_breed.pdf")
ggviolin(
  y = "x",
  ylab = "PC1",
  x = "BREED",
  add = c("jitter", error.plot = "crossbar"),
  data = Combodataframe,
  draw_quantiles = c(0.025, 0.5, 0.975)
)
ggviolin(
  y = "y",
  ylab = "PC1",
  x = "BREED",
  add = c("jitter", error.plot = "crossbar"),
  data = Combodataframe,
  draw_quantiles = c(0.025, 0.5, 0.975)
)
dev.off()


#bayes factor analysis

#create a for loop for results of LES by breed


ECOPHYS_results <- data.frame()
#colnames(ECOPHYS_results)<-c("BF","BFerror","mean","mean.SD","quant.2.5","quant.97.5")

for (i in 3:4) {
  trait.range <-
    range(Combodataframe[, i])[2] - range(Combodataframe[, i])[1]
  
  Trait <- colnames(Combodataframe)[i]
  #bftest has a error but we need the loop to continue
  bf <-
    ttestBF(formula = reformulate(
      termlabels = "BREED",
      response = colnames(Combodataframe)[i]
    ),
    data = Combodataframe)
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
  row.names(results) <- Trait
  colnames(results) <-
    c("BF", "BFerror", "mean", "mean.SD", "quant.2.5", "quant.97.5")
  ECOPHYS_results <- rbind(ECOPHYS_results, results)
}
#add them to large results file

View(ECOPHYS_results)

write.csv(ECOPHYS_results,
          "OUTPUT/RelativeRatios_bayesFactor_BREED.csv")



#violin plots and bf of oil vs nonoil

Combodataframe <- cbind(PCA.IND$data, SAMMETA)

Combodataframe$USE <- as.factor(Combodataframe$USE)

USE <- c("Oil", "NonOil")

Combodataframe <- Combodataframe[which(Combodataframe$USE == USE), ]



#violin plots and bf of oil vs nonoil

pdf("OUTPUT/RelativeRatios_violinplots_use.pdf")
ggviolin(
  y = "x",
  ylab = "PC1",
  x = "USE",
  add = c("jitter", error.plot = "crossbar"),
  data = Combodataframe,
  draw_quantiles = c(0.025, 0.5, 0.975)
)
ggviolin(
  y = "y",
  ylab = "PC2",
  x = "USE",
  add = c("jitter", error.plot = "crossbar"),
  data = Combodataframe,
  draw_quantiles = c(0.025, 0.5, 0.975)
)
dev.off()
#create a for loop for results of LES by breed


ECOPHYS_results <- data.frame()
#colnames(ECOPHYS_results)<-c("BF","BFerror","mean","mean.SD","quant.2.5","quant.97.5")

for (i in 3:4) {
  trait.range <-
    range(Combodataframe[, i])[2] - range(Combodataframe[, i])[1]
  
  Trait <- colnames(Combodataframe)[i]
  #bftest has a error but we need the loop to continue
  bf <-
    ttestBF(formula = reformulate(
      termlabels = "USE",
      response = colnames(Combodataframe)[i]
    ),
    data = Combodataframe)
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
  row.names(results) <- Trait
  colnames(results) <-
    c("BF", "BFerror", "mean", "mean.SD", "quant.2.5", "quant.97.5")
  ECOPHYS_results <- rbind(ECOPHYS_results, results)
}
#add them to large results file

View(ECOPHYS_results)

write.csv(ECOPHYS_results,
          "OUTPUT/RelativeRatios_bayesFactor_USE.csv")




#produce eigenvalue plot

pdf("OUTPUT/RelativeRatios_eigen.pdf")
fviz_eig(RelativeRatios.PCA, choice = "eigenvalue", addlabels = T)
dev.off()


#produce plot of individuals labels by core 12
pdf("OUTPUT/RelativeRatios_IndividualLoading_axes_1_2.pdf")
PCA.IND <-
  fviz_pca_ind(
    RelativeRatios.PCA,
    geom = "point",
    axes = c(1, 2),
    label = "none",
    habillage = SAMMETA$CLASS,
    mean.point = F
  )
ggplot(cbind(PCA.IND$data, SAMMETA),
       aes(x = x,
           y = y,)) +
  geom_point(#color="black",
    aes(
      fill = SAMMETA$CLASS,
      # color="black",
      shape = SAMMETA$CORE12,
      size = SAMMETA$CORE12
    )) +
  scale_size_manual(values = c(5, 2)) +
  scale_shape_manual(values = c(21, 24)) +
  ggtitle(PCA.IND$labels$title) +
  xlab(PCA.IND$labels$x) + ylab(PCA.IND$labels$y) +
  labs(fill = "Breeding group", shape = "Core 12") +
  guides(size = "none", fill = guide_legend(override.aes = list(shape =
                                                                  21))) + geom_hline(yintercept = 0,
                                                                                     linetype = "dashed",
                                                                                     color = "black") + geom_vline(xintercept = 0,
                                                                                                                   linetype = "dashed",
                                                                                                                   color = "black") + theme_bw()
ggplot(cbind(PCA.IND$data, SAMMETA),
       aes(x = x,
           y = y,)) +
  geom_point(#color="black",
    aes(
      fill = SAMMETA$BREED,
      # color="black",
      shape = SAMMETA$CORE12,
      size = SAMMETA$CORE12
    )) +
  scale_size_manual(values = c(5, 2)) +
  scale_shape_manual(values = c(21, 24)) +
  ggtitle(PCA.IND$labels$title) +
  xlab(PCA.IND$labels$x) + ylab(PCA.IND$labels$y) +
  labs(fill = "Breeding group", shape = "Core 12") +
  guides(size = "none", fill = guide_legend(override.aes = list(shape =
                                                                  21))) + geom_hline(yintercept = 0,
                                                                                     linetype = "dashed",
                                                                                     color = "black") + geom_vline(xintercept = 0,
                                                                                                                   linetype = "dashed",
                                                                                                                   color = "black") + theme_bw()
ggplot(cbind(PCA.IND$data, SAMMETA),
       aes(x = x,
           y = y,)) +
  geom_point(#color="black",
    aes(
      fill = SAMMETA$COUNTRY,
      # color="black",
      shape = SAMMETA$CORE12,
      size = SAMMETA$CORE12
    )) +
  scale_size_manual(values = c(5, 2)) +
  scale_shape_manual(values = c(21, 24)) +
  ggtitle(PCA.IND$labels$title) +
  xlab(PCA.IND$labels$x) + ylab(PCA.IND$labels$y) +
  labs(fill = "Breeding group", shape = "Core 12") +
  guides(size = "none", fill = guide_legend(override.aes = list(shape =
                                                                  21))) + geom_hline(yintercept = 0,
                                                                                     linetype = "dashed",
                                                                                     color = "black") + geom_vline(xintercept = 0,
                                                                                                                   linetype = "dashed",
                                                                                                                   color = "black") + theme_bw()
ggplot(cbind(PCA.IND$data, SAMMETA),
       aes(x = x,
           y = y,)) +
  geom_point(#color="black",
    aes(
      fill = SAMMETA$USE,
      # color="black",
      shape = SAMMETA$CORE12,
      size = SAMMETA$CORE12
    )) +
  scale_size_manual(values = c(5, 2)) +
  scale_shape_manual(values = c(21, 24)) +
  ggtitle(PCA.IND$labels$title) +
  xlab(PCA.IND$labels$x) + ylab(PCA.IND$labels$y) +
  labs(fill = "Breeding group", shape = "Core 12") +
  guides(size = "none", fill = guide_legend(override.aes = list(shape =
                                                                  21))) + geom_hline(yintercept = 0,
                                                                                     linetype = "dashed",
                                                                                     color = "black") + geom_vline(xintercept = 0,
                                                                                                                   linetype = "dashed",
                                                                                                                   color = "black") + theme_bw()

dev.off()


#addto HPLC composite Traits

HPLC_Composites <- cbind(HPLC_Composites, RelativeRatios.PCA$x[, 1:10])




###RelativeRatios

#compute PCA

RelativeRatios200.PCA <-
  prcomp(RelativeRatios200, center = T, scale. = T)
PCA.IND <-
  fviz_pca_ind(
    RelativeRatios200.PCA,
    geom = "point",
    axes = c(1, 2),
    label = "none",
    habillage = SAMMETA$CLASS,
    mean.point = F
  )
Combodataframe <- cbind(PCA.IND$data, SAMMETA)
#Bayes factor of each group

Combodataframe <- cbind(PCA.IND$data, SAMMETA)


Combodataframe$BREED <- as.factor(Combodataframe$BREED)

breeds <- c("HA", "RHA")

Combodataframe <- Combodataframe[which(Combodataframe$BREED == breeds), ]
#change between PC1 & PC2
#View(Combodataframe)
#trait.range<-range(Combodataframe$PC1)[2]-range(Combodataframe$PC1)[1]

pdf("OUTPUT/RelativeRatios200_violinplots_breed.pdf")
ggviolin(
  y = "x",
  ylab = "PC1",
  x = "BREED",
  add = c("jitter", error.plot = "crossbar"),
  data = Combodataframe,
  draw_quantiles = c(0.025, 0.5, 0.975)
)
ggviolin(
  y = "y",
  ylab = "PC2",
  x = "BREED",
  add = c("jitter", error.plot = "crossbar"),
  data = Combodataframe,
  draw_quantiles = c(0.025, 0.5, 0.975)
)
dev.off()


#bayes factor analysis

#create a for loop for results of LES by breed


ECOPHYS_results <- data.frame()
#colnames(ECOPHYS_results)<-c("BF","BFerror","mean","mean.SD","quant.2.5","quant.97.5")

for (i in 3:4) {
  trait.range <-
    range(Combodataframe[, i])[2] - range(Combodataframe[, i])[1]
  
  Trait <- colnames(Combodataframe)[i]
  #bftest has a error but we need the loop to continue
  bf <-
    ttestBF(formula = reformulate(
      termlabels = "BREED",
      response = colnames(Combodataframe)[i]
    ),
    data = Combodataframe)
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
  row.names(results) <- Trait
  colnames(results) <-
    c("BF", "BFerror", "mean", "mean.SD", "quant.2.5", "quant.97.5")
  ECOPHYS_results <- rbind(ECOPHYS_results, results)
}
#add them to large results file

View(ECOPHYS_results)

write.csv(ECOPHYS_results,
          "OUTPUT/RelativeRatios200_bayesFactor_BREED.csv")



#violin plots and bf of oil vs nonoil
PCA.IND <-
  fviz_pca_ind(
    RelativeRatios200.PCA,
    geom = "point",
    axes = c(1, 2),
    label = "none",
    habillage = SAMMETA$CLASS,
    mean.point = F
  )

Combodataframe <- cbind(PCA.IND$data, SAMMETA)

Combodataframe$USE <- as.factor(Combodataframe$USE)

USE <- c("Oil", "NonOil")

Combodataframe <- Combodataframe[which(Combodataframe$USE == USE), ]



#violin plots and bf of oil vs nonoil

pdf("OUTPUT/RelativeRatios200_violinplots_use.pdf")
ggviolin(
  y = "x",
  ylab = "PC1",
  x = "USE",
  add = c("jitter", error.plot = "crossbar"),
  data = Combodataframe,
  draw_quantiles = c(0.025, 0.5, 0.975)
)
ggviolin(
  y = "y",
  ylab = "PC2",
  x = "USE",
  add = c("jitter", error.plot = "crossbar"),
  data = Combodataframe,
  draw_quantiles = c(0.025, 0.5, 0.975)
)
dev.off()
#create a for loop for results of LES by breed


ECOPHYS_results <- data.frame()
#colnames(ECOPHYS_results)<-c("BF","BFerror","mean","mean.SD","quant.2.5","quant.97.5")

for (i in 3:4) {
  trait.range <-
    range(Combodataframe[, i])[2] - range(Combodataframe[, i])[1]
  
  Trait <- colnames(Combodataframe)[i]
  #bftest has a error but we need the loop to continue
  bf <-
    ttestBF(formula = reformulate(
      termlabels = "USE",
      response = colnames(Combodataframe)[i]
    ),
    data = Combodataframe)
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
  row.names(results) <- Trait
  colnames(results) <-
    c("BF", "BFerror", "mean", "mean.SD", "quant.2.5", "quant.97.5")
  ECOPHYS_results <- rbind(ECOPHYS_results, results)
}
#add them to large results file

View(ECOPHYS_results)

write.csv(ECOPHYS_results,
          "OUTPUT/RelativeRatios200_bayesFactor_USE.csv")





#produce eigenvalue plot

pdf("OUTPUT/RelativeRatios200_eigen.pdf")
fviz_eig(RelativeRatios200.PCA,
         choice = "eigenvalue",
         addlabels = T)
dev.off()




#produce plot of individuals labels by core 12
pdf("OUTPUT/RelativeRatios200_IndividualLoading_axes_1_2.pdf")


PCA.IND <-
  fviz_pca_ind(
    RelativeRatios200.PCA,
    geom = "point",
    axes = c(1, 2),
    label = "none",
    habillage = SAMMETA$CLASS,
    mean.point = F
  )
ggplot(cbind(PCA.IND$data, SAMMETA),
       aes(x = x,
           y = y,)) +
  geom_point(#color="black",
    aes(
      fill = SAMMETA$CLASS,
      # color="black",
      shape = SAMMETA$CORE12,
      size = SAMMETA$CORE12
    )) +
  scale_size_manual(values = c(5, 2)) +
  scale_shape_manual(values = c(21, 24)) +
  ggtitle(PCA.IND$labels$title) +
  xlab(PCA.IND$labels$x) + ylab(PCA.IND$labels$y) +
  labs(fill = "Breeding group", shape = "Core 12") +
  guides(size = "none", fill = guide_legend(override.aes = list(shape =
                                                                  21))) + geom_hline(yintercept = 0,
                                                                                     linetype = "dashed",
                                                                                     color = "black") + geom_vline(xintercept = 0,
                                                                                                                   linetype = "dashed",
                                                                                                                   color = "black") + theme_bw()
ggplot(cbind(PCA.IND$data, SAMMETA),
       aes(x = x,
           y = y,)) +
  geom_point(#color="black",
    aes(
      fill = SAMMETA$BREED,
      # color="black",
      shape = SAMMETA$CORE12,
      size = SAMMETA$CORE12
    )) +
  scale_size_manual(values = c(5, 2)) +
  scale_shape_manual(values = c(21, 24)) +
  ggtitle(PCA.IND$labels$title) +
  xlab(PCA.IND$labels$x) + ylab(PCA.IND$labels$y) +
  labs(fill = "Breeding group", shape = "Core 12") +
  guides(size = "none", fill = guide_legend(override.aes = list(shape =
                                                                  21))) + geom_hline(yintercept = 0,
                                                                                     linetype = "dashed",
                                                                                     color = "black") + geom_vline(xintercept = 0,
                                                                                                                   linetype = "dashed",
                                                                                                                   color = "black") + theme_bw()
ggplot(cbind(PCA.IND$data, SAMMETA),
       aes(x = x,
           y = y,)) +
  geom_point(#color="black",
    aes(
      fill = SAMMETA$COUNTRY,
      # color="black",
      shape = SAMMETA$CORE12,
      size = SAMMETA$CORE12
    )) +
  scale_size_manual(values = c(5, 2)) +
  scale_shape_manual(values = c(21, 24)) +
  ggtitle(PCA.IND$labels$title) +
  xlab(PCA.IND$labels$x) + ylab(PCA.IND$labels$y) +
  labs(fill = "Breeding group", shape = "Core 12") +
  guides(size = "none", fill = guide_legend(override.aes = list(shape =
                                                                  21))) + geom_hline(yintercept = 0,
                                                                                     linetype = "dashed",
                                                                                     color = "black") + geom_vline(xintercept = 0,
                                                                                                                   linetype = "dashed",
                                                                                                                   color = "black") + theme_bw()
ggplot(cbind(PCA.IND$data, SAMMETA),
       aes(x = x,
           y = y,)) +
  geom_point(#color="black",
    aes(
      fill = SAMMETA$USE,
      # color="black",
      shape = SAMMETA$CORE12,
      size = SAMMETA$CORE12
    )) +
  scale_size_manual(values = c(5, 2)) +
  scale_shape_manual(values = c(21, 24)) +
  ggtitle(PCA.IND$labels$title) +
  xlab(PCA.IND$labels$x) + ylab(PCA.IND$labels$y) +
  labs(fill = "Breeding group", shape = "Core 12") +
  guides(size = "none", fill = guide_legend(override.aes = list(shape =
                                                                  21))) + geom_hline(yintercept = 0,
                                                                                     linetype = "dashed",
                                                                                     color = "black") + geom_vline(xintercept = 0,
                                                                                                                   linetype = "dashed",
                                                                                                                   color = "black") + theme_bw()



dev.off()


#add to HPLC composite Traits

HPLC_Composites <-
  cbind(HPLC_Composites, RelativeRatios200.PCA$x[, 1:10])


#GiveRow names to HPLC composite traits
HPLCCOMPOSITENAMES <-
  c(
    paste0("PEAKAREAS_PC", (seq(
      from = 1, to = 10, by = 1
    ))),
    paste0("PeakAreas200_PC", (seq(
      from = 1, to = 10, by = 1
    ))),
    paste0("RelativeRatios_PC", (seq(
      from = 1, to = 10, by = 1
    ))),
    paste0("RelativeRatios200_PC", (seq(
      from = 1, to = 10, by = 1
    )))
  )

colnames(HPLC_Composites) <- HPLCCOMPOSITENAMES

write.csv(HPLC_Composites, "HPLC_Composites.csv")


#create a summary stats table of data
Combodataframe <- merge(SAMMETA, HPLC_Composites, by = 0)



#make column 1:7 factors

Combodataframe[, 1:7] <- lapply(Combodataframe[, 1:7], factor)
#write summary statistics to CSV


Ecophysummary <-
  describeBy(Combodataframe[, c(4, 9:ncol(Combodataframe))], group = 'BREED', mat =
               TRUE)

write.csv(Ecophysummary,
          "OUTPUT/HPLC_composites_SummarySTATS_BREED.csv")

Ecophysummary <-
  describeBy(Combodataframe[, c(7, 9:ncol(Combodataframe))], group = 'USE', mat =
               TRUE)

write.csv(Ecophysummary, "OUTPUT/HPLC_composites_STATS_USE.csv")
