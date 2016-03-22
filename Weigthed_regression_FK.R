
library(ggplot2)
library(zoo)
library(corrplot)
library(lattice)
require(graphics)
require(grDevices)
library(pheatmap)   
library(gplots)

## Set working directory
setwd("C:/RWEM_processed_calvin/R")

mydata <- read.csv("Quarterly_Data_V2.csv", na.string = 0, stringsAsFactors = FALSE)
# na.string = 0
names(mydata)[2] <- paste("STOPPED")  ### sites stopped
names(mydata)[3] <- paste("IDENTIFIED") ### site identified 
names(mydata)[4] <- paste("ILLEGAL") ### ILLEGAL_EXP_EVENT
names(mydata)[5] <- paste("MIXE")  ### recycled mixed paper
names(mydata)[6] <- paste("NEWS")  ### recycled newspapers
names(mydata)[7] <- paste("GLAS")  ### recycled glas
names(mydata)[8] <- paste("STCN")  ### recycled steel cans
names(mydata)[9] <- paste("ALCN")  ### recycled Alu cans
names(mydata)[10] <- paste("PLBT")  ### recycled plastic bottles
names(mydata)[11] <- paste("HOUSE")  ### house waste
names(mydata)[12] <- paste("RECY")  ### RECYCLE_P
names(mydata)[13] <- paste("DRY")  ### Dry recycle
names(mydata)[14] <- paste("DISP")  ### disposal
names(mydata)[15] <- paste("LA_W")  ### LA collected waste
names(mydata)[16] <- paste("CPI")  ### CPI_NOTNATIONAL_INDICIES_BASE2005
names(mydata)[17] <- paste("RPI")  ### RPI_NOTNATIONAL_INDICIES_BASE1987
names(mydata)[18] <- paste("GDP")  ### GDP_CURRENT_MARKETPRICE
names(mydata)[19] <- paste("POP")  ### POPULATION_ENGLAND
names(mydata)[20] <- paste("DWEL")  ### DWELLINGS_ENGLAND
names(mydata)[21] <- paste("CTDR")  ### TURNOVER(£M)_COLLECTION&TREATMENT&DISPOSAL&RECOVERY
names(mydata)[22] <- paste("WATE")  ### Water Supply; Sewerage, Waste Management and Remediation Activities (Index) (Seasonally adjusted) (2012 index year)
names(mydata)[23] <- paste("WATP")  ### Water Supply; Sewerage,Waste Management & Remediation Act (period-period growth (Seasonally adjusted) index 2012

mydata$YYYY.QQ <- as.yearqtr(mydata$YYYY.QQ)


###### Function to be used in coorrelograms plots ##########################
############################################################################

# mat : is a matrix of data
# ... : further arguments to pass to the native R cor.test function
cor.mtest <- function(mat, ...) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat<- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], ...)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}

##########################################################################

library(corrplot)
pdf(paste("Correlation_Matrix",".pdf",sep=""), width = 17, height = 16.53)
mydata_standard <- data.frame(scale(mydata[2:23]))
M <- cor(mydata[2:23], use = "pairwise.complete.obs")
# M <- cor(mydata_standard, use = "pairwise.complete.obs")
M[is.na(M)]=0
p.mat <- cor.mtest(M)
# head(p.mat[, 1:5])
col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))

corrplot(M, method="color", col=col(200),  
         type="upper", order="hclust", 
         addCoef.col = "black", # Add coefficient of correlation
         tl.col="black", tl.srt=45, tl.cex = 1.2, #Text label color and rotation
         # Combine with significance
         p.mat = p.mat, sig.level = 0.01, insig = "blank", 
      #   p.mat = p.mat, 
         # hide correlation coefficient on the principal diagonal
         diag=FALSE)
dev.off()


####### covariance matrix for M #################################################
#### it is about the relationshp among VARIABLES (columns) of thematrix "mydata" 

# M <- cbind(mydata[,2:18], mydata[,22:23])  ### leave out CTDR, POP and DWEL  
# M <- as.matrix(sapply(M, as.numeric))
M <- mydata[,-1]

# Standardise variables (rescale data based on the meand and Standard deviation)
M_standard <- data.frame(scale(M))

#### covariance matrix #################
M_cov <- cov(M, use='pairwise') 
M_cov_norm <- cov(M_standard, use='pairwise')  #### standardized
write.csv(M_cov, file = "Covar_Matrix.csv", row.names=TRUE)
write.csv(M_cov_norm, file = "Covar_Matrix_Standard.csv", row.names=TRUE)


# levelplot(M_cov)

### Covariance Matrix non normalized
# pdf(paste("Covariance_Matrix_1",".pdf",sep=""), width = 15, height = 15)
jpeg('Covariance_Matrix.jpg',
     quality = 100, bg = "white", res = 200, width = 10, height = 7, units = "in")
par(mar=c(4, 10, 9, 2) + 0.3)
oldpar <- par(las=1)


heatmap <- pheatmap(M_cov, col=bluered(256), cluster_cols=TRUE, 
                    fontsize_row=15,fontsize_col=15, border_color=NA, margins=c(7,7))


dev.off()



### Covariance Matrix normalized
# pdf(paste("Covariance_Matrix_1",".pdf",sep=""), width = 15, height = 15)
jpeg('Covariance_Matrix_normalized.jpg',
     quality = 100, bg = "white", res = 200, width = 10, height = 7, units = "in")
par(mar=c(4, 10, 9, 2) + 0.3)
oldpar <- par(las=1)


heatmap <- pheatmap(M_cov_norm, col=bluered(256), cluster_cols=TRUE, 
                    fontsize_row=15,fontsize_col=15, border_color=NA, margins=c(7,7))


dev.off()



### Covariance Matrix 2
pdf(paste("Covariance_Matrix_a",".pdf",sep=""), width = 17, height = 15)
x  <- as.matrix(M_cov_norm)
rc <- rainbow(nrow(x), start = 0, end = .3)
cc <- rainbow(ncol(x), start = 0, end = .3)
hv <- heatmap(x, col = cm.colors(256), scale = "column",
              RowSideColors = rc, ColSideColors = cc, margins = c(15,12),
              cexRow = 3, cexCol = 3)
dev.off()



# k <- ncol(M) #number of variables
# n <- nrow(M) #number of subjects
# 
# #create means for each column
# M_mean <- matrix(data=1, nrow=n) %*% cbind(mean(M$STOPPED, na.rm=TRUE),
#                                            mean(M$IDENTIFIED, na.rm=TRUE),
#                                            mean(M$ILLEGAL, na.rm=TRUE),
#                                            mean(M$MIXE, na.rm=TRUE),
#                                            mean(M$NEWS, na.rm=TRUE),
#                                            mean(M$GLAS, na.rm=TRUE),
#                                            mean(M$STCN, na.rm=TRUE),
#                                            mean(M$ALCN,na.rm =TRUE),
#                                            mean(M$PLBT, na.rm=TRUE),
#                                            mean(M$HOUSE, na.rm=TRUE),
#                                            mean(M$RECY, na.rm=TRUE),
#                                            mean(M$DRY, na.rm=TRUE),
#                                            mean(M$DISP, na.rm=TRUE),
#                                            mean(M$LA_W, na.rm=TRUE),
#                                            mean(M$CPI, na.rm=TRUE),
#                                            mean(M$RPI, na.rm=TRUE),
#                                            mean(M$GDP, na.rm=TRUE),
#                                            #mean(M$POP, na.rm=TRUE),
#                                            #mean(M$DWEL, na.rm=TRUE),
#                                            #mean(M$CTDR, na.rm=TRUE),
#                                            mean(M$WATE, na.rm=TRUE),
#                                            mean(M$WATP, na.rm=TRUE))
# 
# 
#  M_mean <- as.matrix(sapply(M_mean, as.numeric))
#                                            
# #creates a difference matrix
# D <- M - M_mean
# D <- as.matrix(sapply(D, as.numeric))
# 
# 
# #creates the covariance matrix
# C <- ((n-1)^-1) * t(D) %*% D


########## Multilienar Regression ############################################
##############################################################################


 fit <- lm(STOPPED ~ MIXE + NEWS + GLAS + STCN + ALCN + PLBT + HOUSE + RECY +
           DRY + DISP + LA_W + CPI + RPI + GDP + WATE + WATP,
           data = mydata, na.action= na.exclude)  ### leave out CTDR, POP and DWEL  


# fit <- lm(IDENTIFIED ~ MIXE + NEWS + GLAS + STCN + ALCN + PLBT + HOUSE + RECY +
#             DRY + DISP + LA_W + CPI + RPI + GDP + WATE + WATP, 
#           data = mydata)  ### leave out CTDR, POP and DWEL  (few observations)

# fit <- lm(ILLEGAL ~ MIXE + NEWS + GLAS + STCN + ALCN + PLBT + HOUSE + RECY +
#             DRY + DISP + LA_W + CPI + RPI + GDP + WATE + WATP, 
#           data = mydata)  ### leave out CTDR, POP and DWEL  (few observations)


names <- c("RECYCLE_MIXEDPAPER_£", "RECYCLE_NEWSPAPER_£", "RECYCLE_MIXGLASS_£",
           "RECYCLE_STEELCANS_£", "RECYCLE_ALUMCANS_£", "RECYCLE_PLASTICBOTTLE_£",
           "HOUSEWASTE_T", "RECYCLE_P","DRYRECYCLE_P", "DISPOSAL_P","LA_COLLECTWASTE_T",
           "CPI_NOTNATIONAL_INDICIES_BASE2005", "RPI_NOTNATIONAL_INDICIES_BASE1987",
           "GDP_CURRENT_MARKETPRICE", "E : Water Supply; Sewerage, Waste Management and Remediation Activities (Index) (Seasonally adjusted) (2012 index year)",
           "Water Supply; Sewerage,Waste Management & Remediation Act (period-period growth (Seasonally adjusted) index 2012")
names <- as.data.frame(names)

summary(fit) # show results
signif(summary(fit)$r.squared, 5) ### R2
 

### Diagnostic plots ##### 
Coeff <- as.data.frame(coefficients(fit)) # model coefficients (intercef and slope from the regression)

 confint(fit, level=0.95) # CIs for model parameters 
 fitted(fit) # predicted values
 residuals(fit) # residuals
 anova(fit) # anova table 
 Covariance_Matrix <-  vcov(fit) # covariance matrix for model parameters 
 influence(fit) # regression diagnostics

# Calculate Relative Importance for Each Predictor
# Function to calculate relative importance metrics for linear models,
# calc.relimp calculates several relative importance metrics for the linear model.
# The recommended metrics are lmg (R^2 partitioned by averaging over orders,
# like in Lindemann, Merenda and Gold (1980, p.119ff)) and pmvd
# (a newly proposed metric by Feldman (2005) that is provided
# in the non-US version of the package only). For completeness and comparison
# purposes, several other metrics are also on offer (cf. e.g. Darlington (1968)).

library(relaimpo)

# pdf(paste("Relative_Importance_Site_Stopped",".pdf",sep=""), width = 20, height = 16.53)
# pdf(paste("Relative_Importance_Site_Identified",".pdf",sep=""), width = 20, height = 16.53)
pdf(paste("Relative_Importance_Site_Illegal",".pdf",sep=""), width = 20, height = 16.53)

# IMPO_Sites_Stopped <- calc.relimp(fit, type="lmg", rela=TRUE)
# IMPO_Sites_Identified <- calc.relimp(fit, type="lmg", rela=TRUE)
IMPO_Sites_Illegal <- calc.relimp(fit, type="lmg", rela=TRUE)

# plot(IMPO_Sites_Stopped, sort = FALSE,
#     main ="Relative Importance of Variables vs Sites Stopped")

# plot(IMPO_Sites_Identified, sort = FALSE,
#      main ="Relative Importance of Variables vs Sites Identified")

plot(IMPO_Sites_Illegal, sort = FALSE,
     main ="Relative Importance of Variables vs Sites Illegal")

dev.off()



# R2_IMPO <-  as.data.frame(IMPO_Sites_Stopped$lmg)*100  ### % of R2
# R2_IMPO <-  as.data.frame(IMPO_Sites_Identified$lmg)*100  ### % of R2
R2_IMPO <-  as.data.frame(IMPO_Sites_Illegal$lmg)*100  ### % of R2

# Ranking_IMPO <- as.data.frame(IMPO_Sites_Stopped$lmg.rank)  ### Rank of R2
# Ranking_IMPO <- as.data.frame(IMPO_Sites_Identified$lmg.rank)  ### Rank of R2
Ranking_IMPO <- as.data.frame(IMPO_Sites_Illegal$lmg.rank)  ### Rank of R2


IMPO_Stats <- cbind(names, Coeff[-1,], R2_IMPO, Ranking_IMPO)
IMPO_Stats <- cbind(Row.Names = rownames(IMPO_Stats), IMPO_Stats)
colnames(IMPO_Stats) <- c("Codes", "names", "Coeff. (slope)", "Importance_R2", "Ranking_R2")

# write.csv(IMPO_Stats, file = "Sites_Stopped_Importance_Statisitcs.csv", row.names=TRUE)
# write.csv(IMPO_Stats, file = "Sites_Identified_Importance_Statisitcs.csv", row.names=TRUE)
write.csv(IMPO_Stats, file = "Sites_Illegal_Importance_Statisitcs.csv", row.names=TRUE)

