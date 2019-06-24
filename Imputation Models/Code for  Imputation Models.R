#creating a dataset for the other apply functions to use as direction for number of iterations
j <- data.frame(c(seq(0,2000000, 1)))

#defining functions that perform different types of imputation.
#in all cases we are comparing random samples of 7 values (test) of a normal distribution. 
#In each case we impute three more values and examine on average what the p-values of a two sample comparison look like.

    #rnorm10 takes the mean and sd of the original sample and attempts to emulate the original distribution as best as possible.
    #Using rnorm, rnorm10 draws 10 sets of 3 values using the mean and sd of the original distribution.
    #Of these 10 options it then imputes the three values that change the mean and sd of the original distribution of values the least.
    #(mean and sd are weighted equally for this)
    f_modelimp_rnorm10 <- function(x) {
      set.seed(x)
      if(x %% 10000 == 0.0000000000000000) {
        print(x/10000)
      }
    
      set.seed(x)
      test1 <- c(rnorm(7))
      preextendstats <- c(mean(test1), sd(test1))
      
      winnercandidates <- c(rep(0, 4))
      
      for(i in 1:10) {
        set.seed(x+(i*nrow(j)))
        testextend <- c(rnorm(mean = mean(test1), sd = sd(test1), 3))
        test1a <- c(test1, testextend)
        
        postextendstats <- c(mean(test1a), sd(test1a))
        impdiff <- (postextendstats-preextendstats)
        impdiff <- (sum(abs(impdiff)))
        
        candidate <- c(impdiff, testextend)
        winnercandidates <- cbind(winnercandidates, candidate)
        
        
      }
      
      winnercandidates <- data.frame(winnercandidates)
      winnercandidates <- winnercandidates[, 2:ncol(winnercandidates)]
      winnercol <- which.min(winnercandidates[1, ])
      winner <- winnercandidates[2:nrow(winnercandidates), winnercol]
      
      test1 <- c(test1, winner)
      
      set.seed(-x)
      test2 <- c(rnorm(7))
      preextendstats <- c(mean(test2), sd(test2))
      
      winnercandidates <- c(rep(0, 4))
      
      for(i in 1:10) {
        set.seed(-x-(i*nrow(j)))
        testextend <- c(rnorm(mean = mean(test2), sd = sd(test2), 3))
        test2a <- c(test2, testextend)
        
        postextendstats <- c(mean(test2a), sd(test2a))
        impdiff <- (postextendstats-preextendstats)
        impdiff <- (sum(abs(impdiff)))
        
        candidate <- c(impdiff, testextend)
        winnercandidates <- cbind(winnercandidates, candidate)
        
        
      }
      
      winnercandidates <- data.frame(winnercandidates)
      winnercandidates <- winnercandidates[, 2:ncol(winnercandidates)]
      winnercol <- which.min(winnercandidates[1, ])
      winner <- winnercandidates[2:nrow(winnercandidates), winnercol]
      
      test2 <- c(test2, winner)
      
      t_testresults <- t.test(test1, test2, var.equal = TRUE)
      return(c(t_testresults["p.value"]))
    }
    
    
    # this "perfect" imputation just draws 10 values from a normal distribution. 
    # The three "imputed" values are just more values drawn from the original distribution.
    # While it is impossible to achieve this form of imputation in reality, it is a good baseline to compare other methods to.
    f_modelimp_perf <- function(x) {
      set.seed(x)
      if(x %% 10000 == 0.0000000000000000) {
        print(x/10000)
      }
      
      set.seed(x)
      test1 <- c(rnorm(10))
      set.seed(-x)
      test2 <- c(rnorm(10))
      
      t_testresults <- t.test(test1, test2, var.equal = TRUE)
      return(c(t_testresults["p.value"]))
    }
    
    # this method of imputation imputes the three values with the mean of the sample.
    f_modelimp_mean <- function(x) {
      if(x %% 10000 == 0.0000000000000000) {
        print(x/10000)
      }
      
      set.seed(x)
      test1 <- c(rnorm(7))
      testextend <- c(rep(mean(test1), 3))
      test1 <- c(test1, testextend)
      
      set.seed(-x)
      test2 <- c(rnorm(7))
      testextend <- c(rep(mean(test2), 3))
      test2 <- c(test2, testextend)
      
      t_testresults <- t.test(test1, test2, var.equal = TRUE)
      return(c(t_testresults["p.value"]))
    }
    
    # this method of imputation uses the mean and sd of the sample to approximate the original distribution.
    # It draws three values from the approximated distribution and uses those as the imputed values.
    # this method is very similar to rnorm10, rnorm10 essentially just does this process repeatedly and selects for random values that don't alter the original sample's mean and sd too much.
    f_modelimp_rnorm <- function(x) {
      set.seed(x)
      if(x %% 10000 == 0.0000000000000000) {
        print(x/10000)
      }
      
      set.seed(x)
      test1 <- c(rnorm(7))
      set.seed(x + (nrow(j)))
      testextend <- c(rnorm(mean = mean(test1), sd = sd(test1), 3))
      test1 <- c(test1, testextend)
      
      set.seed(-x)
      test2 <- c(rnorm(7))
      set.seed(-x - (nrow(j)))
      testextend <- c(rnorm(mean = mean(test2), sd = sd(test2), 3))
      test2 <- c(test2, testextend)
      
      t_testresults <- t.test(test1, test2, var.equal = TRUE)
      return(c(t_testresults["p.value"]))  
      
    }
    
    # this model of hot deck imputation assumes that there is another sample available that perfectly represents the sample that has missing values.
    # in reality, we draw 10 values from rnorm and then sample a few random values from it. This is our "foreign sample" that perfectly represents the sample with NAs.
    # then we shorten the "perfect sample" by three values (we "lose"  those values, they become NAs we need to impute). This is now our sample that has NAs.
    # we now have a sampling of values from a "perfect sample" that represents the reality of our sample with NAs. These become the imputed values for our sample.
    # In a real hot-deck imputation we would need to look within our own dataset for this perfect sample that very well reflects the sample that has NAs.
    # The assumption that we'd find a perfect sample to represent every one of our samples with NAs is a very large one, 
    # especially because strong correlation does not necessarily indicate identical distribution between samples. 
    # However, I don't know how else to model hot deck imputation than in this ideal scenario.
    f_modelimp_idealhotdeck <- function(x) {
      set.seed(x)
      if(x %% 10000 == 0.0000000000000000) {
        print(x/10000)
      }
      
      set.seed(x)
      test1 <- c(rnorm(10))
      set.seed(x + (nrow(j)))
      testextend <- c(sample(test1, 3, replace = TRUE))
      test1 <- test1[1:7]
      test1 <- c(test1, testextend)
      
      set.seed(-x)
      test2 <- c(rnorm(10))
      set.seed(x + (nrow(j)))
      testextend <- c(sample(test2, 3, replace = TRUE))
      test2 <- test2[1:7]
      test2 <- c(test2, testextend)
      
      t_testresults <- t.test(test1, test2, var.equal = TRUE)
      return(c(t_testresults["p.value"]))
    }
    
    # this method of imputation imputes the three values with the median of the sample.
    f_modelimp_median <- function(x) {
      set.seed(x)
      if(x %% 10000 == 0.0000000000000000) {
        print(x/10000)
      }
      
      test1 <- c(rnorm(7))
      testextend <- c(rep(mean(test1), 3))
      test1 <- c(test1, testextend)
      
      test2 <- c(rnorm(7))
      testextend <- c(rep(mean(test2), 3))
      test2 <- c(test2, testextend)
      
      t_testresults <- t.test(test1, test2, var.equal = TRUE)
      return(c(t_testresults["p.value"]))
    }

    
# runnning our functions through 2 million iterations to see what their distribution of p values end up looking like 
# (we want to see which option preserves the structure of the data most and creates the least false positives)
idealhotdeckimputation <- apply(j, 1, f_modelimp_idealhotdeck)
meanimputation <- apply(j, 1, f_modelimp_mean)
medianimputation <- apply(j, 1, f_modelimp_median)
perfectimputation <- apply(j, 1, f_modelimp_perf)
rnorm10imputation <- apply(j, 1, f_modelimp_rnorm10)
rnormimputation <- apply(j, 1, f_modelimp_rnorm)

# formatting
idealhotdeckimputation <- unlist(idealhotdeckimputation)
meanimputation <- unlist(meanimputation)
medianimputation <- unlist(medianimputation)
perfectimputation <- unlist(perfectimputation)
rnorm10imputation <- unlist(rnorm10imputation)
rnormimputation <- unlist(rnormimputation)

# generating two values that will help compare the best and worst options
meanmode <- (ecdf(meanimputation) (0.005))*length((meanimputation))
perfectaverage <- (nrow(j)/200)

# histograms for a visual comparison that show off the merits of some methods against others
hist(idealhotdeckimputation, 200, ylim = c(0,meanmode))
  abline(h = meanmode)
  abline(h = perfectaverage)

hist(meanimputation, 200, ylim = c(0,meanmode))
  abline(h = meanmode)
  abline(h = perfectaverage)

hist(medianimputation, 200, ylim = c(0,meanmode))
  abline(h = meanmode)
  abline(h = perfectaverage)

hist(perfectimputation, 200, ylim = c(0,meanmode))
  abline(h = meanmode)
  abline(h = perfectaverage)

hist(rnorm10imputation, 200, ylim = c(0,meanmode))
  abline(h = meanmode)
  abline(h = perfectaverage)

hist(rnormimputation, 200, ylim = c(0,meanmode))
  abline(h = meanmode)
  abline(h = perfectaverage)

#removing uneccesary data
remove(meanmode, perfectaverage, j)