rm(list=ls())
source("BKAGL_train.R")
source("BKAGL_test.R")

#input the training set and test set, 
TRAINING <- ??
TEST <- ??
#group features in each dimension

#set the hyper-parameters
parameters <- list()
parameters$alpha_lambda <- 1
parameters$beta_lambda <- 1
parameters$alpha_eta <- 1
parameters$beta_eta <- 1
parameters$alpha_gamma <- 1
parameters$beta_gamma <- 1
#the above three groups of parameters can be tuned to control the sparsity of the sample, group and dimension weights respectively
parameters$alpha_omega <- 1
parameters$beta_omega <- 1

#set the number of iterations
parameters$iteration <- 1000
  
#set the margin parameter
parameters$margin <- 1

#set the seed
parameters$seed <- 100
  
#set the standard deviation of intermediate representations
parameters$sigma_g <- 1

#calculate the kernel matrix set of the training set
Ktrain <- ??
ytrain <- ??

#training the model
state <- BKAGL_train(Ktrain, ytrain, parameters)

#display the learning parameters of dimensions and groups
print(state$ec$mu[-1])
print(state$b)

#calculate the kernel matrix set of the test set
Ktest <- ??

#predict the class labels of the test set
prediction <- BKAGL_test(Ktest, state)

#display the predicted probabilities
print(prediction$p)