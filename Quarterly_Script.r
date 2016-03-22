# install.packages("zoo")
# install.packages("ggplot2")


# library(ggplot2)
# library(zoo)

## Set working directory
# setwd("C:/Users/cj5/Desktop/Waste_TS/R")
setwd("C:/RWEM_processed_calvin/R")


#### STEP 1: IMPORT DATA - Set Zeros as Null Values 
mydata <- read.csv("Quarterly_Data_v2.csv", na.string = 0,  stringsAsFactors = FALSE)



# Set Index Class for Quarterly Data
library(zoo)
mydata$YYYY.QQ <- as.yearqtr(mydata$YYYY.QQ)
# Add date field for time-series
mydata$ts_date <- as.POSIXct(mydata$YYYY.QQ, tz = "Etc/GMT-1")



#### STEP 2: Produce time-series and regression plots
require(gridExtra)
require(ggplot2)

# Return the names of df columns (Remove first and last fields)
variables <- colnames(mydata)
variables <- variables[2:(length(variables)-1)]


# Time-Series Plots
	p1_name <- paste(variables[1])  ### Site stopped
	p2_name <- paste(variables[2])  ### Site identified
	p3_name <- paste(variables[3])  ### Illegal exported waste
	
	### Site stopped	
	
	p1 <- ggplot(mydata, aes_string("ts_date", p1_name)) + 
				geom_line(colour="black", size=1.5) +  	# Black lines
				geom_point(size=3.5, colour="red") +    # Red dots
				theme_bw() +							# Change background theme to white with grey grids
				xlab("") + ylab(p1_name) +
				ggtitle(p1_name)
	
	### Site identified

	p2 <- ggplot(mydata, aes_string("ts_date", p2_name)) + 
				geom_line(colour="steelblue", size=1.5) +  	# Blue lines
				geom_point(size=3.5, colour="black") +      # Black dots
				theme_bw() +								# Change background theme to white with grey grids
				xlab("") + ylab(p2_name) +
				ggtitle(p2_name)
	
	### Illegal exported waste
	
	p3 <- ggplot(mydata, aes_string("ts_date", p3_name)) + 
	  geom_line(colour="red", size=1.5) +  	# Blue lines
	  geom_point(size=3.5, colour="black") +      # Black dots
	  theme_bw() +								# Change background theme to white with grey grids
	  xlab("") + ylab(p3_name) +
	  ggtitle(p3_name)

	Site_Stopped <- p1
	Site_Identified <- p2
	Ill_Exp_Waste <- p3
	
	## Linear regression plot
	eq <- paste(paste("mydata$",p1_name,sep=""),"~", paste("mydata$",p2_name,sep=""))
	reg <- lm(eq)																# R basic regression to obtain performance information
	
	fit_lm <- ggplot(mydata, aes_string(p1_name, p2_name)) +
	  geom_point(shape=1, size=5) +    							# Use hollow circles
	  theme_bw() +												# Change background theme to white with grey grid
	  geom_smooth(method=lm, colour="red", size=1) +   			# Add linear regression line (includes 95% confidence)
	  xlab(p1_name) + ylab(p2_name) +
	  ggtitle(expression(bold("Waste Measurements"))) +
	  
	  # Add regression performance information ot title
	  labs(title = paste("Adj R2 = ",signif(summary(reg)$adj.r.squared, 5),
	                     "  Intercept =",signif(reg$coef[[1]],5 ),
	                     "  Slope =",signif(reg$coef[[2]], 5),
	                     "  P =",signif(summary(reg)$coef[2,4], 5)))	
	

	##  Export Plots
	pdf(paste(p2_name,".pdf",sep=""), width = 23.39, height = 16.53)		# Output set to A2 sheet dimensions								
	print(grid.arrange(arrangeGrob(p1, p2), fit_lm, ncol=2))				# require(gridExtra)
	dev.off()	


	
#### new Stuff from Federico ######	
	
	
	index <- c(1:22)
	# LIST_ALL <- list(variables)
	
	for (i in 1:22) {
	 
	p_name <- paste(variables[index[i]])  ### Site stopped
	
	mypath <- file.path("C:","RWEM_processed_calvin","R","plots",
	                    paste(p_name, index[i], ".jpg", sep = ""))

	p <- ggplot(mydata, aes_string("ts_date", p_name)) + 
	  geom_line(colour="black", size=1.5) +  	# Black lines
	  geom_point(size=3.5, colour="red") +    # Red dots
	  theme_bw() +							# Change background theme to white with grey grids
	  xlab("") + ylab(p_name) +
	  ggtitle(p_name)
	
	ggsave(mypath, p)
	}
