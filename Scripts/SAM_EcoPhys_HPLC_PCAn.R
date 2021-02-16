###SAM Ecophys and HPLC data 
library(factoextra)
library(missMDA)
library(corrplot)
library(ggfortify)
library(BayesFactor)
library(ggpubr)
#import data
getwd()
setwd("DATA/")
getwd()
ECOPHYS<- read.csv("SAMecophys_linemeans.csv", row.names = 1)

#change PPN to row names & remove last column

ECOPHYS<-ECOPHYS[,-21]



#Read in Metadata
SAMMETA<-read.csv("SAMPANEL_Metadatalabels.csv", row.names = 1) 



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
ECOPHYS_Completeset<-read.csv("Ecophys_imputed.csv",row.names = 1)


#create dataset to store Composite Traits

ECOPHY_Composites<- data.frame()
#subset phenolics and flavonoid


FlavPhen_Traits<- ECOPHYS_Completeset[,c(18,19)]

                              #   c(3,9,10,12,15,17)]

#PCA of flavenoids and phenolics


FlavPhen.PCA<-prcomp(FlavPhen_Traits,center = T,scale. =T)


#produce eigenvalue plot

pdf("OUTPUT/FlavPhen_eigen.pdf")
fviz_eig(FlavPhen.PCA, choice = "eigenvalue", addlabels = T)
dev.off()


#produce variable loading plots
pdf("OUTPUT/FlavPhen_VariableLoading_axes_1_2.pdf")
fviz_pca_var(FlavPhen.PCA,axes = c(1,2), label="var",
             alpha.var = "cos2",
             col.var = "cos2", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)
dev.off()

#produce plot of individuals labels by core 12
pdf("OUTPUT/FlavPhen_IndividualLoading_axes_1_2.pdf")
PCA.IND<-fviz_pca_ind(FlavPhen.PCA,geom.ind = "point",pointshape=19, label=SAMMETA$CORE12,col.ind = SAMMETA$GROUP, mean.point=F)
ggplot(cbind(PCA.IND$data,SAMMETA),aes(x=x,y=y,col=SAMMETA$GROUP,shape=SAMMETA$CORE12,size=SAMMETA$CORE12)) + geom_point() + scale_size_manual(values=c(4,2))+
  ggtitle(PCA.IND$labels$title)+xlab(PCA.IND$labels$x)+ylab(PCA.IND$labels$y)+labs(col="Breeding group",shape="Core 12")+guides(size=FALSE)+ geom_hline(yintercept=0, linetype="dashed", color = "black")+ geom_vline(xintercept=0, linetype="dashed", color = "black")+theme_bw()
dev.off()

#assess cos2 

var <- get_pca_var(FlavPhen.PCA)


#create plot for Cos2
pdf("OUTPUT/FlavPhen_cos2.pdf")
corrplot(var$cos2,method="color",p.mat =round(var$cos2,digits = 2), insig = "p-value", sig.level = 0.01,  outline = T,cl.lim = c(0,1), is.corr = FALSE)
dev.off()




#add data to composite data set

ECOPHY_Composites <- FlavPhen.PCA$x



#subset LES traits


LES_Traits<- ECOPHYS_Completeset[,c(1:17)]
LES_Traits <- LES_Traits[,-c(5,8,13)]


#PCA of LES TRaits


LES_Traits.PCA<-prcomp(LES_Traits,center = T,scale. =T)
summary(LES_Traits.PCA)

#Bayes factor of each group

Combodataframe<-merge(LES_Traits.PCA$x[,1:3],SAMMETA, by=0)

Combodataframe$BREED<-as.factor(Combodataframe$BREED)

breeds<-c("HA", "RHA")

Combodataframe<-Combodataframe[which(Combodataframe$BREED==breeds),]
#change between PC1 & PC2
View(Combodataframe)
range(Combodataframe$PC1)[1]-range(Combodataframe$PC1)[2]

ggviolin(y="PC1",x= "BREED",add = c("jitter", error.plot = "crossbar"), data = Combodataframe, draw_quantiles = c(0.025,0.5,0.975))
bf <- ttestBF(formula = PC1 ~ BREED, data = Combodataframe)
bf
chains <- posterior(bf, iterations = 10000)
summary(chains)
plot(chains)




#produce eigenvalue plot

pdf("OUTPUT/LES_eigen.pdf")
fviz_eig(LES_Traits.PCA, choice = "eigenvalue", addlabels = T)
dev.off()


#produce variable loading plots
pdf("OUTPUT/LES_VariableLoading_axes_1_4.pdf")
fviz_pca_var(LES_Traits.PCA,axes = c(1,2), label="var",
             alpha.var = "cos2",
             col.var = "cos2", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)
dev.off()

#produce plot of individuals labels by core 12
pdf("OUTPUT/LES_IndividualLoading_axes_1_2.pdf")
PCA.IND <-
  fviz_pca_ind(
    LES_Traits.PCA,
    axes = c(1, 2),
    geom.ind = "point",
    pointshape = 19,
    label = SAMMETA$CORE12,
    col.ind = SAMMETA$GROUP,
    mean.point = F
  )
ggplot(
  cbind(PCA.IND$data, SAMMETA),
  aes(
    x = x,
    y = y,
    fill = SAMMETA$GROUP,
    colour="black",
    shape = SAMMETA$CORE12,
    size = SAMMETA$CORE12
  )
) + 
  geom_point(colour="black") + 
  scale_size_manual(values = c(5, 2)) +
  scale_shape_manual(values = c(21, 24))+
  ggtitle(PCA.IND$labels$title) + 
  xlab(PCA.IND$labels$x) + ylab(PCA.IND$labels$y) +
  labs(fill = "Breeding group", shape = "Core 12") + 
  guides(size = FALSE) + geom_hline(yintercept =0,linetype = "dashed",color = "black") + geom_vline(xintercept = 0,linetype = "dashed",color = "black") + theme_bw()
dev.off()

#assess cos2 

var <- get_pca_var(LES_Traits.PCA)


#create plot for Cos2
pdf("OUTPUT/LES_cos2.pdf")
corrplot(var$cos2[,1:10],method="color",p.mat =round(var$cos2[,1:10],digits = 2), insig = "p-value", sig.level = (1/nrow(var$cos2)),  outline = T,cl.lim = c(0,1), is.corr = FALSE)
dev.off()



#add data to composite data set

ECOPHY_Composites <- cbind(FlavPhen.PCA$x,LES_Traits.PCA$x[,1:10])




colnames(ECOPHY_Composites) <- c("FLAVPHEN_PC1","FLAVPHEN_PC2","LES_PC1","LES_PC2","LES_PC3","LES_PC4","LES_PC5","LES_PC6","LES_PC7","LES_PC8","LES_PC9","LES_PC10")

#save ecophys composite traits
write.csv(ECOPHY_Composites, "Ecophys_CompositeTraits.csv")









###Sam HPLC PCA


#IMPORT ALL DATA SETS
PeakAreas<- read.csv("SAMHPLC_PLOTPEAKAREAS_ROUNDED.csv", header = T, row.names = 1)
PeakAreas200<-read.csv("SAMHPLC_Retentiontimes200_PlotPeakAreas.csv", header = T, row.names = 1)
RelativeRatios<-read.csv("SAMHPLC_RelativeRatios.csv", header = T, row.names = 1)
RelativeRatios200<- read.csv("SAMHPLC_Retentiontimes200_RelativeRatios.csv", header = T, row.names = 1)

#create Composite dataset to store 10axes each 

HPLC_Composites<- data.frame()

#change all NA to 0

PeakAreas[is.na(PeakAreas)] <- 0
PeakAreas200[is.na(PeakAreas200)] <- 0
RelativeRatios[is.na(RelativeRatios)] <- 0
RelativeRatios200[is.na(RelativeRatios200)] <- 0



#compute PCA 

PeakAreas.PCA<-prcomp(PeakAreas,center = T,scale. =T)

#produce eigenvalue plot

pdf("OUTPUT/PeakAreas_eigen.pdf")
fviz_eig(PeakAreas.PCA, choice = "eigenvalue", addlabels = T)
dev.off()

#Bayes factor of each group

Combodataframe<-merge(PeakAreas.PCA$x[,1:3],SAMMETA, by=0)

Combodataframe$BREED<-as.factor(Combodataframe$BREED)

breeds<-c("HA", "RHA")

Combodataframe<-Combodataframe[which(Combodataframe$BREED==breeds),]
#change between PC1 & PC2

range(Combodataframe$PC1)[1]-range(Combodataframe$PC1)[2]

ggviolin(y="PC1",x= "BREED",add = c("jitter", error.plot = "crossbar"), data = Combodataframe, draw_quantiles = c(0.025,0.5,0.975))
bf <- ttestBF(formula = PC1 ~ BREED, data = Combodataframe)
bf
chains <- posterior(bf, iterations = 10000)
summary(chains)
plot(chains)






#produce plot of individuals labels by core 12
pdf("OUTPUT/PeakAreas_IndividualLoading_axes_1_2.pdf")
PCA.IND<-fviz_pca_ind(PeakAreas.PCA,geom = "point",axes = c(1,2), label="none", habillage= ECOPHYS$CORE12, mean.point=F)
ggplot(
  cbind(PCA.IND$data, SAMMETA),
  aes(
    x = x,
    y = y,
    fill = SAMMETA$GROUP,
    colour="black",
    shape = SAMMETA$CORE12,
    size = SAMMETA$CORE12
  )
) + 
  geom_point(colour="black") + 
  scale_size_manual(values = c(5, 2)) +
  scale_shape_manual(values = c(21, 24))+
  ggtitle(PCA.IND$labels$title) + 
  xlab(PCA.IND$labels$x) + ylab(PCA.IND$labels$y) +
  labs(fill = "Breeding group", shape = "Core 12") + 
  guides(size = FALSE) + geom_hline(yintercept =0,linetype = "dashed",color = "black") + geom_vline(xintercept = 0,linetype = "dashed",color = "black") + theme_bw()
dev.off()
#addto HPLC composite Traits 

HPLC_Composites<- PeakAreas.PCA$x[,1:10]

###PEAK AREAS >200

#compute PCA 

PeakAreas200.PCA<-prcomp(PeakAreas200,center = T,scale. =T)

#Bayes factor of each group

Combodataframe<-merge(PeakAreas200.PCA$x[,1:3],SAMMETA, by=0)

Combodataframe$BREED<-as.factor(Combodataframe$BREED)

breeds<-c("HA", "RHA")

Combodataframe<-Combodataframe[which(Combodataframe$BREED==breeds),]
#change between PC1 & PC2
View(Combodataframe)
range(Combodataframe$PC1)[1]-range(Combodataframe$PC1)[2]

ggviolin(y="PC1",x= "BREED",add = c("jitter", error.plot = "crossbar"), data = Combodataframe, draw_quantiles = c(0.025,0.5,0.975))
bf <- ttestBF(formula = PC1 ~ BREED, data = Combodataframe)
bf
chains <- posterior(bf, iterations = 10000)
summary(chains)
plot(chains)



#produce eigenvalue plot

pdf("PeakAreas200_eigen.pdf")
fviz_eig(PeakAreas200.PCA, choice = "eigenvalue", addlabels = T)
dev.off()


#produce plot of individuals labels by core 12
pdf("PeakAreas200_IndividualLoading_axes_1_2.pdf")
PCA.IND<-fviz_pca_ind(PeakAreas200.PCA,geom = "point",axes = c(1,2), label="none", habillage= ECOPHYS$CORE12, mean.point=F)
ggplot(
  cbind(PCA.IND$data, SAMMETA),
  aes(
    x = x,
    y = y,
    fill = SAMMETA$GROUP,
    colour="black",
    shape = SAMMETA$CORE12,
    size = SAMMETA$CORE12
  )
) + 
  geom_point(colour="black") + 
  scale_size_manual(values = c(5, 2)) +
  scale_shape_manual(values = c(21, 24))+
  ggtitle(PCA.IND$labels$title) + 
  xlab(PCA.IND$labels$x) + ylab(PCA.IND$labels$y) +
  labs(fill = "Breeding group", shape = "Core 12") + 
  guides(size = FALSE) + geom_hline(yintercept =0,linetype = "dashed",color = "black") + geom_vline(xintercept = 0,linetype = "dashed",color = "black") + theme_bw()

dev.off()


#addto HPLC composite Traits 

HPLC_Composites<- cbind(HPLC_Composites,PeakAreas200.PCA$x[,1:10])


###RelativeRatios > 200

#compute PCA 

RelativeRatios.PCA<-prcomp(RelativeRatios,center = T,scale. =T)

#Bayes factor of each group

Combodataframe<-merge(RelativeRatios.PCA$x[,1:3],SAMMETA, by=0)

Combodataframe$BREED<-as.factor(Combodataframe$BREED)

breeds<-c("HA", "RHA")

Combodataframe<-Combodataframe[which(Combodataframe$BREED==breeds),]
#change between PC1 & PC2
View(Combodataframe)
range(Combodataframe$PC1)[1]-range(Combodataframe$PC1)[2]

ggviolin(y="PC1",x= "BREED",add = c("jitter", error.plot = "crossbar"), data = Combodataframe, draw_quantiles = c(0.025,0.5,0.975))
bf <- ttestBF(formula = PC1 ~ BREED, data = Combodataframe)
bf
chains <- posterior(bf, iterations = 10000)
summary(chains)
plot(chains)


#produce eigenvalue plot

pdf("OUTPUT/RelativeRatios_eigen.pdf")
fviz_eig(RelativeRatios.PCA, choice = "eigenvalue", addlabels = T)
dev.off()


#produce plot of individuals labels by core 12
pdf("OUTPUT/RelativeRatios_IndividualLoading_axes_1_2.pdf")
PCA.IND<-fviz_pca_ind(RelativeRatios.PCA,geom = "point",axes = c(1,2), label="none", habillage= ECOPHYS$CORE12, mean.point=F)
ggplot(
  cbind(PCA.IND$data, SAMMETA),
  aes(
    x = x,
    y = y,
    fill = SAMMETA$GROUP,
    colour="black",
    shape = SAMMETA$CORE12,
    size = SAMMETA$CORE12
  )
) + 
  geom_point(colour="black") + 
  scale_size_manual(values = c(5, 2)) +
  scale_shape_manual(values = c(21, 24))+
  ggtitle(PCA.IND$labels$title) + 
  xlab(PCA.IND$labels$x) + ylab(PCA.IND$labels$y) +
  labs(fill = "Breeding group", shape = "Core 12") + 
  guides(size = FALSE) + geom_hline(yintercept =0,linetype = "dashed",color = "black") + geom_vline(xintercept = 0,linetype = "dashed",color = "black") + theme_bw()

dev.off()


#addto HPLC composite Traits 

HPLC_Composites<- cbind(HPLC_Composites,RelativeRatios.PCA$x[,1:10])




###RelativeRatios

#compute PCA 

RelativeRatios200.PCA<-prcomp(RelativeRatios200,center = T,scale. =T)

#Bayes factor of each group

Combodataframe<-merge(RelativeRatios200.PCA$x[,1:3],SAMMETA, by=0)

Combodataframe$BREED<-as.factor(Combodataframe$BREED)

breeds<-c("HA", "RHA")

Combodataframe<-Combodataframe[which(Combodataframe$BREED==breeds),]
#change between PC1 & PC2
View(Combodataframe)
range(Combodataframe$PC1)[1]-range(Combodataframe$PC1)[2]

ggviolin(y="PC1",x= "BREED",add = c("jitter", error.plot = "crossbar"), data = Combodataframe, draw_quantiles = c(0.025,0.5,0.975))
bf <- ttestBF(formula = PC1 ~ BREED, data = Combodataframe)
bf
chains <- posterior(bf, iterations = 10000)
summary(chains)
plot(chains)




#produce eigenvalue plot

pdf("OUTPUT/RelativeRatios200_eigen.pdf")
fviz_eig(RelativeRatios200.PCA, choice = "eigenvalue", addlabels = T)
dev.off()


#produce plot of individuals labels by core 12
pdf("OUTPUT/RelativeRatios200_IndividualLoading_axes_1_2.pdf")


PCA.IND<-fviz_pca_ind(RelativeRatios200.PCA,geom = "point",axes = c(1,2), label="none", habillage = SAMMETA$GROUP, mean.point=F)
ggplot(
  cbind(PCA.IND$data, SAMMETA),
  aes(
    x = x,
    y = y,
    fill = SAMMETA$GROUP,
    colour="black",
    shape = SAMMETA$CORE12,
    size = SAMMETA$CORE12
  )
) + 
  geom_point(colour="black") + 
  scale_size_manual(values = c(5, 2)) +
  scale_shape_manual(values = c(21, 24))+
  ggtitle(PCA.IND$labels$title) + 
  xlab(PCA.IND$labels$x) + ylab(PCA.IND$labels$y) +
  labs(fill = "Breeding group", shape = "Core 12") + 
  guides(size = FALSE) + geom_hline(yintercept =0,linetype = "dashed",color = "black") + geom_vline(xintercept = 0,linetype = "dashed",color = "black") + theme_bw()



dev.off()


#add to HPLC composite Traits 

HPLC_Composites<- cbind(HPLC_Composites,RelativeRatios200.PCA$x[,1:10])


#GiveRow names to HPLC composite traits
HPLCCOMPOSITENAMES <- c(paste0("PEAKAREAS_PC",(seq(from=1,to=10,by=1))),paste0("PeakAreas200_PC",(seq(from=1,to=10,by=1))),
                        paste0("RelativeRatios_PC",(seq(from=1,to=10,by=1))),paste0("RelativeRatios200_PC",(seq(from=1,to=10,by=1))))

colnames(HPLC_Composites)<-HPLCCOMPOSITENAMES

write.csv(HPLC_Composites,"HPLC_Composites.csv")



