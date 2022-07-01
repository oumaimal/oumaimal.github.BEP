library(remulate)
library(remstats)
library(remify)
library(relevent)

# library(magrittr)
# library(tidyverse)

fas_data <- read.csv('C:/Users/20193694/Documents/Quartile 4/BEP/Code/reduced_riskset_sentiment.csv')
names(fas_data) <- c('time', 'actor1', 'actor2', 'tweet', 'weight')

without_sentiment <- read.csv('C:/Users/20193694/Documents/Quartile 4/BEP/Code/reduced_riskset.csv')
names(without_sentiment) <- c('time', 'actor1', 'actor2', 'tweet')


library(Rcpp)

# Cpp function to compute lambda
cppFunction("
	
	arma::mat compute_lambda(arma::cube statistics, arma::colvec coef) {
	
		arma::mat lambda(statistics.n_rows, statistics.n_cols, arma::fill::zeros);
		
		for(arma::uword i = 0; i < coef.n_elem; ++i) {
			arma::mat statsslice = statistics.slice(i);
			arma::mat lambdaslice = coef(i)*statsslice;
			lambda += lambdaslice;
		}
		
		lambda = exp(lambda);
		return lambda;
	
	}", depends = "RcppArmadillo")

# Function to compute gof
compute_gof <- function(fit, statistics, evls, interactions) {
  # Stap 1: Compute lambda
  lambda <- compute_lambda(statistics, fit$coef)
  
  # Stap 2: Compute event probabilities
  eventprob <- t(apply(lambda, 1, function(x) {
    x/sum(x)
  }))
  
  # Stap 3: Inverse event ranks
  eventranks <- t(apply(eventprob, 1, function(x) {
    rank(x, ties.method = "max")/length(x)
  }))
  
  # Stap 4: Gof
  gof <- vector()
  for(j in 1:length(unique(interactions))) {
    top5 <- which(eventranks[which(interactions == unique(interactions)[j])[1],] >= 0.99)
    dyads <- evls[which(interactions == unique(interactions)[j]),1]
    gof[j] <- mean(dyads %in% top5)>0
  }
  
  mean(gof)
}








### MODEL 0
effects.stats <- ~ 1 

statsObject <- remstats::remstats(edgelist = fas_data, 
                                  tie_effects = effects.stats, 
                                  attributes = attributesObject)

fit <- relevent::rem(eventlist = statsObject$evls, statslist = statsObject$statistics,  estimator = "MLE", timing = "interval")
summary(fit)


temp <- eventseq[!duplicated(fas_data$time),]
cr <- sapply(temp$numofactors, function(x) {
  1-(0.99^choose(x, 2))
})
gof0 <- mean(cr)
gof0




### MODEL 1
effects.stats <- ~ 1 +
  remstats::inertia(scaling = 'std') + 
  remstats::reciprocity(scaling = 'std') + 
  remstats::totaldegreeSender(scaling = 'std') + 
  remstats::totaldegreeReceiver(scaling = 'std') 


statsObject <- remstats::remstats(edgelist = fas_data, 
                                  tie_effects = effects.stats)

fit <- relevent::rem(eventlist = statsObject$evls, statslist = statsObject$statistics,  estimator = "MLE", timing = "interval")
summary(fit)


statistics <- statsObject$statistics
evls <- statsObject$evls

gof1 <- compute_gof(fit, statistics, evls, fas_data$time)
gof1




### MODEL 1 - Trying to add inertia weighted by sentiment
effects.stats <- ~ 1 +
  remstats::inertia(scaling = 'std') + 
  remstats::reciprocity(scaling = 'std') + 
  remstats::totaldegreeSender(scaling = 'std') + 
  remstats::totaldegreeReceiver(scaling = 'std') 


statsObject <- remstats::remstats(edgelist = fas_data, 
                                  tie_effects = effects.stats)
sentiment_inertia <- statsObject$statistics[,,2]

effects.stats <- ~ 1 +
  remstats::inertia(scaling = 'std') + 
  remstats::reciprocity(scaling = 'std') + 
  remstats::totaldegreeSender(scaling = 'std') + 
  remstats::totaldegreeReceiver(scaling = 'std') 

statsObject <- remstats::remstats(edgelist = without_sentiment, 
                                  tie_effects = effects.stats)
stats <- statsObject$statistics
stats <- abind::abind(stats, sentiment_inertia, along = 3)
dimnames(stats)[[3]] <- c("baseline", "inertia", "reciprocity", "totaldegreeSender", "totaldegreeReceiver", "sentiment_inertia")

fit <- relevent::rem(eventlist = statsObject$evls, statslist = stats,  estimator = "MLE", timing = "interval")
summary(fit)


statistics <- stats
evls <- statsObject$evls

gof1_sent <- compute_gof(fit, statistics, evls, fas_data$time)
gof1_sent






### MODEL 2
effects.stats <- ~ 1 +
  remstats::inertia(scaling = 'std') + 
  remstats::reciprocity(scaling = 'std') + 
  remstats::isp(scaling = 'std') +
  remstats::osp(scaling = 'std')


statsObject <- remstats::remstats(edgelist = fas_data, 
                                  tie_effects = effects.stats)

fit <- relevent::rem(eventlist = statsObject$evls, statslist = statsObject$statistics,  estimator = "MLE", timing = "interval")
summary(fit)


statistics <- statsObject$statistics
evls <- statsObject$evls

gof2 <- compute_gof(fit, statistics, evls, fas_data$time)
gof2




### MODEL 2 - Trying to add inertia weighted by sentiment
effects.stats <- ~ 1 +
  remstats::inertia(scaling = 'std') + 
  remstats::reciprocity(scaling = 'std') + 
  remstats::isp(scaling = 'std') +
  remstats::osp(scaling = 'std')


statsObject <- remstats::remstats(edgelist = fas_data, 
                                  tie_effects = effects.stats)
sentiment_inertia <- statsObject$statistics[,,2]

effects.stats <- ~ 1 +
  remstats::inertia(scaling = 'std') + 
  remstats::reciprocity(scaling = 'std') + 
  remstats::isp(scaling = 'std') +
  remstats::osp(scaling = 'std')

statsObject <- remstats::remstats(edgelist = without_sentiment, 
                                  tie_effects = effects.stats)
stats <- statsObject$statistics
stats <- abind::abind(stats, sentiment_inertia, along = 3)
dimnames(stats)[[3]] <- c("baseline", "inertia", "reciprocity", "isp", "osp", "sentiment_inertia")

fit <- relevent::rem(eventlist = statsObject$evls, statslist = stats,  estimator = "MLE", timing = "interval")
summary(fit)


statistics <- stats
evls <- statsObject$evls

gof2_sent <- compute_gof(fit, statistics, evls, fas_data$time)
gof2_sent





### MODEL 3
effects.stats <- ~ 1 +
  remstats::inertia(scaling = 'std') + 
  remstats::reciprocity(scaling = 'std') + 
  remstats::outdegreeSender(scaling = 'std') + 
  remstats::indegreeReceiver(scaling = 'std')  


statsObject <- remstats::remstats(edgelist = fas_data, 
                                  tie_effects = effects.stats)

fit <- relevent::rem(eventlist = statsObject$evls, statslist = statsObject$statistics,  estimator = "MLE", timing = "interval")
summary(fit)


statistics <- statsObject$statistics
evls <- statsObject$evls

gof3 <- compute_gof(fit, statistics, evls, fas_data$time)
gof3




### MODEL 3 - Trying to add inertia weighted by sentiment
effects.stats <- ~ 1 +
  remstats::inertia(scaling = 'std') + 
  remstats::reciprocity(scaling = 'std') + 
  remstats::outdegreeSender(scaling = 'std') + 
  remstats::indegreeReceiver(scaling = 'std') 


statsObject <- remstats::remstats(edgelist = fas_data, 
                                  tie_effects = effects.stats)
sentiment_inertia <- statsObject$statistics[,,2]

effects.stats <- ~ 1 +
  remstats::inertia(scaling = 'std') + 
  remstats::reciprocity(scaling = 'std') + 
  remstats::outdegreeSender(scaling = 'std') + 
  remstats::indegreeReceiver(scaling = 'std') 

statsObject <- remstats::remstats(edgelist = without_sentiment, 
                                  tie_effects = effects.stats)
stats <- statsObject$statistics
stats <- abind::abind(stats, sentiment_inertia, along = 3)
dimnames(stats)[[3]] <- c("baseline", "inertia", "reciprocity", "outdegreeSender", "indegreeReceiver", "sentiment_inertia")

fit <- relevent::rem(eventlist = statsObject$evls, statslist = stats,  estimator = "MLE", timing = "interval")
summary(fit)


statistics <- stats
evls <- statsObject$evls

gof3_sent <- compute_gof(fit, statistics, evls, fas_data$time)
gof3_sent





### MODEL 4
effects.stats <- ~ 1 +
  remstats::inertia(scaling = 'std') + 
  remstats::reciprocity(scaling = 'std') + 
  remstats::totaldegreeSender(scaling = 'std') + 
  remstats::totaldegreeReceiver(scaling = 'std') + 
  remstats::isp(scaling = 'std') +
  remstats::osp(scaling = 'std') 


statsObject <- remstats::remstats(edgelist = fas_data, 
                                  tie_effects = effects.stats)

fit <- relevent::rem(eventlist = statsObject$evls, statslist = statsObject$statistics,  estimator = "MLE", timing = "interval")
summary(fit)


statistics <- statsObject$statistics
evls <- statsObject$evls

gof4 <- compute_gof(fit, statistics, evls, fas_data$time)
gof4





### MODEL 4 - Trying to add inertia weighted by sentiment
effects.stats <- ~ 1 +
  remstats::inertia(scaling = 'std') + 
  remstats::reciprocity(scaling = 'std') + 
  remstats::totaldegreeSender(scaling = 'std') + 
  remstats::totaldegreeReceiver(scaling = 'std') + 
  remstats::isp(scaling = 'std') +
  remstats::osp(scaling = 'std') 


statsObject <- remstats::remstats(edgelist = fas_data, 
                                  tie_effects = effects.stats)
sentiment_inertia <- statsObject$statistics[,,2]

effects.stats <- ~ 1 +
  remstats::inertia(scaling = 'std') + 
  remstats::reciprocity(scaling = 'std') + 
  remstats::totaldegreeSender(scaling = 'std') + 
  remstats::totaldegreeReceiver(scaling = 'std') + 
  remstats::isp(scaling = 'std') +
  remstats::osp(scaling = 'std') 

statsObject <- remstats::remstats(edgelist = without_sentiment, 
                                  tie_effects = effects.stats)
stats <- statsObject$statistics
stats <- abind::abind(stats, sentiment_inertia, along = 3)
dimnames(stats)[[3]] <- c("baseline", "inertia", "reciprocity", "totaldegreeSender", "totaldegreeReceiver", "isp", "osp", "sentiment_inertia")

fit <- relevent::rem(eventlist = statsObject$evls, statslist = stats,  estimator = "MLE", timing = "interval")
summary(fit)


statistics <- stats
evls <- statsObject$evls

gof4_sent <- compute_gof(fit, statistics, evls, fas_data$time)
gof4_sent



gof_all <- cbind(c(gof1, gof2, gof3, gof4), c(gof1_sent, gof2_sent, gof3_sent, gof4_sent))
colnames(gof_all) <- c("gof", "gof_sentiment")
rownames(gof_all) <- c("Model1", "Model2", "Model3", "Model4")
gof_all








?remstats
