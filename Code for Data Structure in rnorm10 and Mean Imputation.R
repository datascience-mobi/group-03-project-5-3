# this code is iintended to illustrate how differing imputation methods change the structure of the data they are used on
# as measures for data structure we take mean and sd (via t test) as well as order of ranked values (via wilcox rank sum test)
# each time we randomly generate 10 values from a normal distribution, "lose" three, impute three more and compare the pre imputation and post imputation sets with our two tests.
# the average p value in each case is considered a measure for similarity of mean, sd and order of ranked values.
# which in turn each describe the similarity of the imputed and non imputed data structure


# creating a data frame for other apply functions to ppiggyback on
j <- data.frame(c(seq(0,2000000, 1)))

# the rnorm10 as described in "code for imputation models", once for t and once for wilcox
f_modelimp_rnorm10 <- function(x) {
  if(x %% 1000 == 0.0000000000000000) {
    print(x/10000)
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
  return(c(t_testresults["p.value"], w_testresults["p.value"], impdiff))
}


# the mean imputation as described in "code for imputation models", once for t and once for wilcox
f_modelimp_mean <- function(x) {
  if(x %% 10000 == 0.0000000000000000) {
    print(x/1000)
  }
  
  set.seed(x)
  test1 <- c(rnorm(10))
  test2 <- test1[1:7]
  testextend <- c(rep(mean(test2), 3))
  test2 <- c(test2, testextend)
  
  t_testresults <- t.test(test1, test2, var.equal = TRUE)
  w_testresults <- wilcox.test(test1, test2)
  impdiff <- mean(abs(sort(test1[8:10])-sort(testextend)))
  return(c(t_testresults["p.value"], w_testresults["p.value"], impdiff))
}

meanimputation <- data.frame(matrix(unlist(apply(j, 1, f_modelimp_mean)), ncol = 3))
rnorm10imputation <- data.frame(matrix(unlist(apply(j, 1, f_modelimp_rnorm10)), ncol = 3))

#calculating mean p value
meant <- mean(meanimputation[, 1])
meanw <- mean(meanimputation[, 2])
meani <- mean(meanimputation[, 3])
rnorm10t <- mean(rnorm10imputation[, 1])
rnorm10w <- mean(rnorm10imputation[, 2])
rnorm10i <- mean(rnorm10imputation[, 3])

#defining the perfect case scenario, in which we compare two identical vectors, so p is always 1
perfectt <- 1
perfectw <- 1
perfecti <- 0

#formatting and pressenting results
average_p_value_depending_on_imputation_t <- data.frame(c(perfectt, rnorm10t, meant))
average_p_value_depending_on_imputation_w <- data.frame(c(perfectw, rnorm10w, meanw))
average_p_value_depending_on_imputation_i <- data.frame(c(perfecti, rnorm10i, meani))
average_p_value_depending_on_imputation <- data.frame(cbind(average_p_value_depending_on_imputation_t, average_p_value_depending_on_imputation_w, average_p_value_depending_on_imputation_i))

colnames(average_p_value_depending_on_imputation) <- c("T-Test", "Wilcoxon", "Average difference of imputed data points to real data points")
rownames(average_p_value_depending_on_imputation) <- c("perfect", "rnorm10", "mean")

View(average_p_value_depending_on_imputation)

#removing unecessary data
remove(j, meani, meant, meanw, perfecti, perfectt, perfectw, rnorm10i, rnorm10t, rnorm10w, average_p_value_depending_on_imputation_i, average_p_value_depending_on_imputation_t, average_p_value_depending_on_imputation_w)