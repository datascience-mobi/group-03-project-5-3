# this code is intended to illustrate how differing imputation methods change the structure of the data they are used on
# as measures for data structure we take mean and sd (via t test) as well as order of ranked values (via wilcox rank sum test)
# each time we randomly generate 10 values from a normal distribution, "lose" three, impute three more and compare the pre imputation and post imputation sets with our two tests.
# the average p value in each case is considered a measure for similarity of mean, sd and order of ranked values, which in turn each describe the similarity of the imputed and non imputed data structure
# additionally we calculate the average difference between the imputed and real data points when they are each ordered from smallest to largest, which indicates how close the imputed data points are to the original values
# lastly we check how well the mean and sd of the new distribution emulates the real distribution (mean = 0, sd= 1) by calculating the absolute difference between the mean and 0 and the sd and 1.


# creating a data frame for other apply functions to ppiggyback on
j <- data.frame(c(seq(1,5e5, 1)))
jcount <- (nrow(j)/100)

# the same rnorm10 method as described in "code for imputation models", it calculates a t-test, wilcoxon rank sum test and mean difference in values between the imputed vector and original vecctor, as well as how well the imputed distribution mimics the real one.
f_modelimp_rnorm10 <- function(x) {
  if(x %% jcount == 0) {
    print(paste0(x/jcount, " % completed"))
  }
  
  set.seed(x)
  test1 <- c(rnorm(10))
  test2 <- test1[1:7]
  preextendstats <- c(mean(test2), sd(test2))
  
  winnercandidates <- c()
  
  for(i in 1:10) {
    set.seed(x+(i*nrow(j)))
    testextend <- c(rnorm(mean = mean(test2), sd = sd(test2), 3))
    test2a <- c(test2, testextend)
    
    postextendstats <- c(mean(test2a), sd(test2a))
    impdiff <- (postextendstats-preextendstats)
    impdiff <- (sum(abs(impdiff)))
    
    candidate <- c(impdiff, testextend)
    winnercandidates <- cbind(winnercandidates, candidate)
    
    
  }
  
  winnercandidates <- data.frame(winnercandidates)
  winnercol <- which.min(winnercandidates[1, ])
  winner <- winnercandidates[2:nrow(winnercandidates), winnercol]
  
  test2 <- c(test2, winner)
  
  t_testresults <- t.test(test1, test2, var.equal = TRUE)
  w_testresults <- wilcox.test(test1, test2)
  impdiff <- mean(abs(sort(test1[8:10])-sort(winner)))
  impshift_m <- abs(0 - mean(test2))
  impshift_s <- abs(1 - sd(test2))
  return(c(t_testresults["p.value"], w_testresults["p.value"], impdiff, impshift_m, impshift_s))
}


# the mean imputation as described in "code for imputation models", it calculates a t-test, wilcoxon rank sum test and mean difference in values between the imputed vector and original vecctor, as well as how well the imputed distribution mimics the real one.
f_modelimp_mean <- function(x) {
  if(x %% jcount == 0) {
    print(paste0(x/jcount, " % completed"))
  }
  
  set.seed(x)
  test1 <- c(rnorm(10))
  test2 <- test1[1:7]
  testextend <- c(rep(mean(test2), 3))
  test2 <- c(test2, testextend)
  
  t_testresults <- t.test(test1, test2, var.equal = TRUE)
  w_testresults <- wilcox.test(test1, test2)
  impdiff <- mean(abs(sort(test1[8:10])-sort(testextend)))
  impshift_m <- abs(0 - mean(test2))
  impshift_s <- abs(1 - sd(test2))
  return(c(t_testresults["p.value"], w_testresults["p.value"], impdiff, impshift_m, impshift_s))
}

#meanimputation <- data.frame(matrix(unlist(apply(j, 1, f_modelimp_mean)), ncol = 3, byrow = TRUE))
#rnorm10imputation <- data.frame(matrix(unlist(apply(j, 1, f_modelimp_rnorm10)), ncol = 3, byrow = TRUE))

meanimputation <- apply(j, 1, f_modelimp_mean)
meanimputation <- unlist(meanimputation)
meanimputation <- matrix(meanimputation, ncol = 5, byrow = TRUE)
meanimputation <- data.frame(meanimputation)

rnorm10imputation <- apply(j, 1, f_modelimp_rnorm10)
rnorm10imputation <- unlist(rnorm10imputation)
rnorm10imputation <- matrix(rnorm10imputation, ncol = 5, byrow = TRUE)
rnorm10imputation <- data.frame(rnorm10imputation)

#calculating mean values
meant <- mean(meanimputation[, 1])
meanw <- mean(meanimputation[, 2])
meani <- mean(meanimputation[, 3])
meanm <- mean(meanimputation[, 4])
means <- mean(meanimputation[, 5])

rnorm10t <- mean(rnorm10imputation[, 1])
rnorm10w <- mean(rnorm10imputation[, 2])
rnorm10i <- mean(rnorm10imputation[, 3])
rnorm10m <- mean(rnorm10imputation[, 4])
rnorm10s <- mean(rnorm10imputation[, 5])

#defining the perfect case scenario, in which we compare two identical vectors, which means both p values are always 1 and the difference is always 0.
#as for how well the real sample represents the real distribution, a function is defined for that 

perfectt <- 1
perfectw <- 1
perfecti <- 0

f_modelimp_perfect <- function(x) {
  if(x %% jcount == 0) {
    print(paste0(x/jcount, " % completed"))
  }
  
  set.seed(x)
  test1 <- rnorm(10)
  impshift_m <- abs(0-mean(test1))
  impshift_s <- abs(1- sd(test1))
  return(c(impshift_m, impshift_s))
}

perfectimputation <- apply(j, 1, f_modelimp_perfect)
perfectimputation <- unlist(perfectimputation)
perfectimputation <- matrix(perfectimputation, ncol = 2, byrow = TRUE)
perfectimputation <- data.frame(perfectimputation)

perfectm <- mean(perfectimputation[, 1])
perfects <- mean(perfectimputation[, 2])

#formatting and presenting results
average_p_value_depending_on_imputation_t <- data.frame(c(perfectt, rnorm10t, meant))
average_p_value_depending_on_imputation_w <- data.frame(c(perfectw, rnorm10w, meanw))
average_p_value_depending_on_imputation_i <- data.frame(c(perfecti, rnorm10i, meani))
average_p_value_depending_on_imputation_m <- data.frame(c(perfectm, rnorm10m, meanm))
average_p_value_depending_on_imputation_s <- data.frame(c(perfects, rnorm10s, means))

average_p_value_depending_on_imputation <- data.frame(cbind(average_p_value_depending_on_imputation_t, average_p_value_depending_on_imputation_w, average_p_value_depending_on_imputation_i, average_p_value_depending_on_imputation_m, average_p_value_depending_on_imputation_s))

colnames(average_p_value_depending_on_imputation) <- c("T-Test", "Wilcoxon", "Shift of imp. data points", "Shift of the Mean", "Shift of SD")
rownames(average_p_value_depending_on_imputation) <- c("perfect", "rnorm10", "mean")

average_p_value_depending_on_imputation <- t(average_p_value_depending_on_imputation)

View(average_p_value_depending_on_imputation)

#removing unecessary data
remove(j, meani, meant, meanw, perfecti, perfectt, perfectw, rnorm10i, rnorm10t, rnorm10w, average_p_value_depending_on_imputation_i, average_p_value_depending_on_imputation_t, average_p_value_depending_on_imputation_w, average_p_value_depending_on_imputation_m, average_p_value_depending_on_imputation_s, jcount, meanm, means, perfectm, perfects, rnorm10m, rnorm10s)