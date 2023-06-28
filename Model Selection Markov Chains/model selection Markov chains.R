# Settings ####

rm(list = ls())
library(tidyverse)
library(plotly)
source("C:/code/MDL/Model Selection and the Principle of Minimum Description Length/StandardUniversalCodeIntegers.R")

maxLogOrder <- 6
maxLogOrderTry <- 5
nObservations <- 1e3
nSimulations <- 1e3
d <- 4 # the precision: theta is discretised with width 2^-d


# Simulate the Markov chains ####
GetTransitionMatrix <- function(maxLogOrder = 0){
  
  if(maxLogOrder == 0)
    transitionMatrix <- tibble(1)
  else{
    transitionMatrix <- crossing(0:1) %>% set_names(1)
    if(maxLogOrder > 1){
      for(iLogOrder in 2:maxLogOrder){
        transitionMatrix <- crossing(0:1, transitionMatrix) %>% set_names(1:iLogOrder)
      }  
    }
  }
  transitionMatrix <-
    transitionMatrix %>%
    set_names(rev(colnames(transitionMatrix))) %>% 
    mutate(from = paste0(!!!syms( rev( colnames(transitionMatrix) ) ) ),
           total = 0,
           `to 0 (#)` = 0,
           `to 1 (#)` = 0,
           `to 0 (%)` = 0,
           `to 1 (%)` = 0)
  if(maxLogOrder == 0){
    transitionMatrix <-
      transitionMatrix %>% select(-`1`, -from) 
  }
  return(transitionMatrix)
}

transitionProbabilities <- 
  GetTransitionMatrix(maxLogOrder) %>% 
  mutate(
    `to 0 (Pr)` = runif(n()), 
    `to 1 (Pr)` = 1 - `to 0 (Pr)`) %>%
  arrange(from)
ixFromStateStart <- sample.int(nrow(transitionProbabilities), 1)
timeSeriesSimulatedRaw <- 
  transitionProbabilities %>%
  slice(ixFromStateStart) %>%
  pull(from) %>%
  str_split('') %>%
  unlist() %>%
  as.numeric()
fromStates <- transitionProbabilities %>% select(1:maxLogOrder) %>% as.matrix()
for(iSimulation in 1:(nSimulations - maxLogOrder)){
  fromState <- tail(timeSeriesSimulatedRaw, n = maxLogOrder)
  ixFromState <- which( rowSums( sweep(fromStates, 2, fromState, '==') ) == maxLogOrder)
  probabilitiesNextState <-
    transitionProbabilities %>% 
    slice(ixFromState) %>%
    select(`to 0 (Pr)`, `to 1 (Pr)`) %>% 
    as.numeric()
  timeSeriesSimulatedRaw <-
    c(timeSeriesSimulatedRaw,
      sample(0:(length(probabilitiesNextState) - 1), 1, prob = probabilitiesNextState))
}



# Calculate the Markov Chain maximum likelihood for each possible order. ####

lag_multiple <- function(x, n_vec){
  map(n_vec, lag, x = x) %>% 
    set_names(paste('lag', n_vec)) %>% 
    as_tibble()
}
timeSeriesSimulated <- 
  tibble(observed = timeSeriesSimulatedRaw) %>%
  mutate(lag_multiple(observed, 1:maxLogOrderTry),
         target = lead(observed, n = 1))
probabilitiesTry <- seq(0, 1, by = 2^(-d)) 
nSuccess <- sum(timeSeriesSimulated$target, na.rm = TRUE)
thetaDiscretised <- seq(2^(-d), 1 - 2^(-d), by = 2^(-d))
thetaDiscretisedMin <- min(thetaDiscretised)
thetaDiscretisedMax <- max(thetaDiscretised)
logLikelihood <- tibble(order = 0:maxLogOrderTry, `log-likelihood` = NA, `code length order` = NA)
for(iOrder in 0:maxLogOrderTry){
  
  print(iOrder)
  
  # Determine the observed transition matrix.
  transitionMatrix <- 
    GetTransitionMatrix(iOrder)
  if(iOrder == 0){
    transitionMatrix <- 
      transitionMatrix %>%
      mutate(total = nObservations, `to 1 (#)` = nSuccess)
  }else{
    for(iRoute in 1:nrow(transitionMatrix)){
      
      route <- transitionMatrix %>% slice(iRoute) %>% select(all_of(1:iOrder)) %>% as.numeric()
      ixSelect <- (iOrder + 1):(nObservations - 1)
      timeSeriesSimulatedRemaining <- timeSeriesSimulated %>% slice(ixSelect)
      for(iLag in 1:length(route)){
        timeSeriesSimulatedRemaining <- 
          timeSeriesSimulatedRemaining %>%
          filter(.data[[paste('lag', iLag)]] == route[iLag])
      }
      transitionMatrix[iRoute, 'total'] <- nrow(timeSeriesSimulatedRemaining)
      transitionMatrix[iRoute, 'to 1 (#)'] <- timeSeriesSimulatedRemaining %>% pull(target) %>% sum
    }
  }
  transitionMatrix <- 
    transitionMatrix %>% 
    mutate(`to 0 (#)` = total - `to 1 (#)`,
           `to 0 (%)` = `to 0 (#)` / total,
           `to 1 (%)` = `to 1 (#)` / total)
  
  # The optimal theta is enclosed in a discretised grid surrounding the maximum likelihood theta. For example, for the zero-th order Markov
  # model, i.e. the Bernoulli model, suppose the maximum likelihood theta is 0.65 and the optimal theta needs to be chosen among {0.2, 0.4,
  # 0.6, 0.8}. Then the discretised value of theta that is closest to the maximum likelihood theta is optimal; in this case 0.6. For a
  # higher dimensional vector of theta the optimal theta is an enclosed hyper dimensional cube inscribing the maximum likelihood theta.
  thetaMaximumLikelihood <- transitionMatrix$`to 1 (%)`
  
  # Determine for each element of the maximum likelihood theta the discretised values of theta that are most close to it.
  isInitialise <- TRUE
  for(iTheta in 1:length(thetaMaximumLikelihood)){
    thetaMaximumLikelihoodHere <- thetaMaximumLikelihood[iTheta]
    
    if(thetaMaximumLikelihoodHere < thetaDiscretisedMin){
      thetaOptimalGridVertex <- thetaDiscretisedMin
    }else if(thetaMaximumLikelihoodHere > thetaDiscretisedMax){
      thetaOptimalGridVertex <- thetaDiscretisedMax
    }else{
      thetaOptimalGridVertex <- 
        c(
          tail( thetaDiscretised[thetaDiscretised <= thetaMaximumLikelihoodHere], n = 1 ),
          head( thetaDiscretised[thetaMaximumLikelihoodHere < thetaDiscretised], n = 1 )
        )
    }
    
    if(isInitialise){
      thetaOptimalGrid <- tibble(`theta 1` = thetaOptimalGridVertex)
      isInitialise <- FALSE
    }else{
      thetaOptimalGrid <- crossing(thetaOptimalGrid, thetaOptimalGridVertex) %>% set_names(paste('theta', 1:iTheta))
    }
  }
  logLikelihoodOptimalGrid <- 
    thetaOptimalGrid %>%
    mutate(`log-likelihood` = 0)
  for(iTheta in ncol(thetaOptimalGrid)){
    columnName <- paste('theta', iTheta)
    logLikelihoodOptimalGrid <-
      logLikelihoodOptimalGrid %>%
      mutate(`log-likelihood` = 
               `log-likelihood` + 
               transitionMatrix$`to 1 (#)` * log2(.data[[columnName]]) + 
               transitionMatrix$`to 0 (#)` * log2(1 - .data[[columnName]]))
  }
  thetaOptimal <- 
    logLikelihoodOptimalGrid %>%
    slice_max(`log-likelihood`) %>%
    slice_head(n = 1) # any row would suffice as they are all minima
  if(iOrder == 0){
    codeLengthOrder <- 1
  }else{
    codeLengthOrder <- StandardUniversalCodeIntegers(iOrder)
  }
  
  logLikelihood <- 
    logLikelihood %>%
    mutate(`log-likelihood` = if_else(order == iOrder, thetaOptimal$`log-likelihood`, `log-likelihood`),
           `code length order` = ifelse(order == iOrder, codeLengthOrder, `code length order`))
  
  print(logLikelihood)
}

ggplotly(
  logLikelihood %>% 
    pivot_longer(-theta, names_to = 'order', values_to = 'log-likelihood') %>% 
    ggplot(aes(x = theta, y = `log-likelihood`, colour = order)) + 
    geom_line()
  )