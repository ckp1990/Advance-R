#prevelance analysi of lecocytozoon
setwd("I:/onedrive/Leuco_analysis/Processed_datas") # set the woking direcotory in need 
Leco_data<-read.csv("UK_Leuco_Chandan_FI Edits.csv",header = T) #reading the data
## a total of 1737 observaton
##attachnng the wrandling libraries. 
library(dplyr)
library(reshape2)
library(tidyr)
##step 1: need to remove the bird species with less than 5 ind
indivdual_count<-group_by(Leco_data,species)%>%
  summarise(count=n())%>%
  as.data.frame()
speceis_retatined<-indivdual_count$species[indivdual_count$count>5]
##selection in the main data base. 
Leco_data<-subset(Leco_data,is.na(match(Leco_data$species,speceis_retatined))==F)
## total of 1579 will be left. 
final_summary_of_data<-Leco_data%>%summarise(n_bird=n_distinct(species),
                      n_bird_genera=n_distinct(genus),
                      n_family=n_distinct(family),
                      n_para=n_distinct(L.SEQ_changed))
##ploting the summary
library(ggplot2)
ggplot(data = melt(final_summary_of_data),aes(x=variable,y=value))+
  geom_bar(stat='identity',fill="Khaki4",color="black")+theme_classic()+
  labs(y="Unique data points")+
  scale_x_discrete(labels=c("Species (Birds)","Genera (Birds)",
                            "Family(Birds)","Lineage (Parasite)"), name="")+
  theme(plot.title = element_text(hjust = 0.5,face = "bold",size = 16))+
  ggtitle("Global summary of data point")

## input dataframe for prevelence. 
prevelance_input_df<-Leco_data[,c(1,2,4,8,9,10,11,13,14,15,16,17,18)]
##the above  mentioned column number may change if different data set
## to begin with 
#only "Sample.No.","Common.Name",Leco_infection" "year" ,"elevation"  ,alt"            "L.SEQ_changed"  "species"       
#"genus","family","order", "class were selected.
## creating the binary format for lecuo infection
Infection_binary<-rep(NA,nrow(prevelance_input_df))

Infection_binary[which(prevelance_input_df$Leco_infection=="Yes")]<-1
Infection_binary[which(prevelance_input_df$Leco_infection=="No")]<-0
prevelance_input_df$Infection_binary<-Infection_binary

## summary of prevelance al different levels. 

Simple_linear_model<-glm(Infection_binary~year+status+species*elevation,
                         data = prevelance_input_df,
                         family = binomial(link = "coglog"))
summary(Simple_linear_model)
as.data.frame(drop1(Simple_linear_model,test = "Chi"))


##linear mixed effect model for prevelance. 
##loading the pacages. 
require(ggplot2)
require(GGally)
require(reshape2)
require(lme4)
require(compiler)
require(parallel)
require(boot)
require(lattice)
##pair wise plots
ggpairs(prevelance_input_df[, c("year", "alt")])
##no Strong corr value 
Explor_prevelance <- melt(prevelance_input_df[, c("Infection_binary", "elevation","status","year")],
						id.vars="Infection_binary")
ggplot(Explor_prevelance, aes(factor(Infection_binary), y = value, fill=factor(Infection_binary))) +
	geom_boxplot() +
	facet_wrap(~variable, scales="free_y")+
	theme_classic()
##fitting the model
##deleting the LM
prevelance_input_df<-prevelance_input_df[prevelance_input_df$status!="LM",]
prevelance_input_df$log_el<-prevelance_input_df$elevation%>%log()
prevelance_input_df$year<-as.factor(prevelance_input_df$year)
Prev_model<- glmer(Infection_binary~log_el+year+status+(1|species)+
                     (1|genus)+(1|family), prevelance_input_df, family = binomial(link = "logit"),
                   control = glmerControl(optimizer = "bobyqa"))

#making plot for random effects for speceis. 
dotplot(ranef(Prev_model, which = "family", condVar = TRUE))
dotplot(ranef(Prev_model, which = "genus", condVar = TRUE))
dotplot(ranef(Prev_model, which = "species", condVar = TRUE))

# print the mod results without correlations among fixed effects
print(Prev_model, corr = FALSE)
summary(Prev_model)
drop1(Prev_model,test="Chi")
#

# print the mod results without correlations among fixed effects
print(Prev_model, corr = FALSE)
summary(Prev_model)
##getting SE and tablular formate
se <- sqrt(diag(vcov(Prev_model)))
# table of estimates with 95% CI
(tab <- cbind(Est = fixef(Prev_model), LL = fixef(Prev_model) - 1.96 * se, UL = fixef(Prev_model) + 1.96 *
								se))

#multilevel bootstrap to smaple and predict the values for ploting 

sampler <- function(dat, clustervar, replace = TRUE, reps = 1) {
	cid <- unique(dat[, clustervar[1]])
	ncid <- length(cid)
	recid <- sample(cid, size = ncid * reps, replace = TRUE)
	if (replace==T) {
		rid <- lapply(seq_along(recid), function(i) {
			cbind(NewID = i, RowID = sample(which(dat[, clustervar] == recid[i]),
																			size = length(which(dat[, clustervar] == recid[i])), replace = TRUE))
		})
	} else {
		rid <- lapply(seq_along(recid), function(i) {
			cbind(NewID = i, RowID = which(dat[, clustervar] == recid[i]))
		})
	}
	dat <- as.data.frame(do.call(rbind, rid))
	dat$Replicate <- factor(cut(dat$NewID, breaks = c(1, ncid * 1:reps), include.lowest = TRUE,
															labels = FALSE))
	dat$NewID <- factor(dat$NewID)
	return(dat)
}
##sampling with the function
set.seed(20)
temp_data<-sampler(prevelance_input_df,"log_el",reps=100)
bigdata <- cbind(temp_data, prevelance_input_df[temp_data$RowID, ])
#Next we refit the model on the resampled data. First we store the estimates from our original model, which we will use as start values for the bootstrap models.
f <- fixef(Prev_model)
r <- getME(Prev_model, "theta")

cl <- makeCluster(3)
clusterExport(cl, c("bigdata", "f", "r"))
clusterEvalQ(cl, require(lme4))

myboot <- function(i) {
	object <- try(glmer(Infection_binary ~ log_el+year+status + (1 | species) + (1 | genus) + (1 | family), data = bigdata, subset = Replicate == i, family = binomial(link="logit"),
											nAGQ = 1, start = list(fixef = f, theta = r)),silent = TRUE)
	if (class(object) == "try-error")
		return(object)
	c(fixef(object), getME(object, "theta"))
}

start <- proc.time()
res <- parLapplyLB(cl, X = levels(bigdata$Replicate), fun = myboot)
end <- proc.time()

# shut down the cluster
stopCluster(cl)
#Now that we have the data, the local cluster, and the fitting function setup, we are ready to actually do the bootstrapping.

#the bootstraps result are there. We can see how many models actually conserved. 
success <- sapply(res, is.numeric)
mean(success)
#Next we convert the list of bootstrap results into a matrix, and then calculate the 2.5th and 97.5th percentiles for each parameter.
# combine successful results
bigres <- do.call(cbind, res[success])

# calculate 2.5th and 97.5th percentiles for 95% CI
(ci <- t(apply(bigres, 1, quantile, probs = c(0.025, 0.975))))

# All results
finaltable <- cbind(Est = c(f, r), SE = c(se, NA), BootMean = rowMeans(bigres),
										ci)
# round and print
round(finaltable, 3)
##ploting the graphs with eleveation as x and predicted value on y
# the catergories I am using is different year and since the status do not say any thingh about the results . We are not dealing with that. 

tmpdat <- prevelance_input_df[, c( "log_el","year","status" , "species", "genus","family")]

jvalues <- with(prevelance_input_df, seq(from = min(log_el), to = max(log_el), length.out = 100))

# calculate predicted probabilities and store in a list
pp <- lapply(jvalues, function(j) {
	tmpdat$log_el <- j
	predict(Prev_model, newdata = tmpdat, type = "response")
})
sapply(pp[c(1, 20, 40, 60, 80, 100)], mean) ## checking is all right

# get the means with lower and upper quartiles
plotdat <- t(sapply(pp, function(x) {
	c(Mean = mean(x), quantile(x, c(0.25, 0.75)))
}))

# add in LengthofStay values and convert to data frame
plotdat <- as.data.frame(cbind(plotdat, jvalues))

# better names and show the first few rows
colnames(plotdat) <- c("PredictedProbability", "Lower", "Upper","elevation")

# plot average marginal predicted probabilities
ggplot(plotdat, aes(x =elevation, y = PredictedProbability)) + geom_ribbon(aes(ymin = Lower,ymax = Upper),alpha=0.5) + geom_line(size = 2) + ylim(c(0, 1))+
	theme_classic()+
	labs(x="Elevation (m) in log scale",y="Probabity of prevalence")
range(plotdat$PredictedProbability)
##for different years. 
# calculate predicted probabilities and store in a list
biprobs <- lapply(levels(prevelance_input_df$year), function(year) {
	tmpdat$year[] <- year
	lapply(jvalues, function(j) {
		tmpdat$log_el <- j
		predict(Prev_model, newdata = tmpdat, type = "response")
	})
})

# get means and quartiles for all jvalues for each level of CancerStage
plotdat2 <- lapply(biprobs, function(X) {
	temp <- t(sapply(X, function(x) {
		c(Mean_value=mean(x), quantile(x, c(.25, .75)))
	}))
	temp <- as.data.frame(cbind(temp, jvalues))
	colnames(temp) <- c("PredictedProbability", "Lower", "Upper", "elevation")
	return(temp)
})

# collapse to one data frame
plotdat2 <- do.call(rbind, plotdat2)

# add year to df
plotdat2$year <- factor(rep(levels(prevelance_input_df$year), each = length(jvalues)))

# show first few rows
head(plotdat2)
##ploting the data 

# graph it
ggplot(plotdat2, aes(x = log(elevation), y = PredictedProbability)) +
	geom_ribbon(aes(ymin = Lower, ymax = Upper,fill=year), alpha = .15) +
	geom_line(aes(colour = year), size = 2) +
	ylim(c(0, 1)) + facet_wrap(~  year)+
	labs(x="Elevation (m) in log scale",y="Probabity of prevalence")+
	theme_classic()

##working with Intensity.

intensity_df<-subset(Leco_data,Leco_data$Leco_infection=="Yes")
intensity_df<-intensity_df[,c(1,2,4,7,8,9,10,11,13,14,15,16,17,18)]
##the above  mentioned column number may change if different data set
## to begin with 
#only "Sample.No.","Common.Name",Leco_itensity" "year" ,"elevation"  ,alt"            "L.SEQ_changed"  "species"       
#"genus","family","order", "class were selected.
## creating the binary format for lecuo infection
intensity_df$Leco_intensity
## summary of prevelance al different levels. 
##linear mixed effect model for prevelance. 
##loading the pacages. 
##deleting the LM
intensity_df<-intensity_df[intensity_df$status!="LM",]
##fitting the model
intensity_df$log_el<-intensity_df$elevation%>%log()
intensity_df$year<-as.factor(intensity_df$year)
##chaning intensity to numberical values
number_intensiy<-data.frame(Leco_intensity=as.character(intensity_df$Leco_intensity)%>%unique(),value=c(0,3,2,1))
#linking to intensity df
intensity_df<-left_join(intensity_df,number_intensiy,by="Leco_intensity")
intensity_model<- lmer(value~alt+year+status+(1|species),data = intensity_df,REML = F)
coefs<-data.frame(coef(summary(intensity_model)))
coefs$p.Z<-round(2*(1-pnorm(abs(coefs$t.value))),3)

summary(intensity_model)
library(lmerTest)
library("emmeans", lib.loc="~/R/win-library/3.5")
#making plot for random effects for speceis. 
dotplot(ranef(intensity_model, which = "species", condVar = TRUE))
anova(intensity_model)
ls_means(intensity_model)
difflsmeans(intensity_model,test.effs="alt+year")
tukey_test_value<-emmeans(intensity_model,list(pairwise~alt+year),adjust="tukey")
tukey_test_value$`pairwise differences of alt, year`
library("effects", lib.loc="~/R/win-library/3.5")
ploting_intensity<-as.data.frame(effect(c("year","alt"),intensity_model))

post_test_result<-data.frame(alt=c("low","medium","high"),
														 text_value=c("*","*","#"))

ploting_intensity<-left_join(ploting_intensity,post_test_result,by="alt")
ploting_intensity$alt<-ordered(as.factor(ploting_intensity$alt),levels=c("low","medium","high"))

intenisty_plot<-ggplot(ploting_intensity,aes(x=alt,y=fit))
dodge<-position_dodge(width = 0.9)
error_bars<-aes(ymax=fit+se,ymin=fit-se)
intenisty_plot+geom_bar(stat = "identity",position = dodge,fill="khaki4")+
	geom_errorbar(error_bars,position = dodge,width=0.5)+
	theme_classic()+geom_text(aes(y=fit+se+0.1,x=alt,label=text_value),size=4)+
	labs(x="Elevation",y="Intenisty of leucocytozoon intfection")+
	facet_wrap(~year)
