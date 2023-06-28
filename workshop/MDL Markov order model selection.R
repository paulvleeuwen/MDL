# Settings ####

rm(list = ls())
library(tidyverse)
library(plotly)

maxMarkovOrder <- 20

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

# Load the raw data. ####

load('DJI_workshop.RData')



# Determine the code length for each Markov order. ####

markovOrder <- 0:maxMarkovOrder
codeLength <- tibble(order = markovOrder, `L(H)` = NA, `L(D|H)` = NA, `L(H, D)` = NA)
nSuccessR <- sum(DJI_workshop$R > 0, na.rm = TRUE)
nSuccessV <- sum(DJI_workshop$V > 0, na.rm = TRUE)
nSuccess <- nSuccessV


lag_multiple <- function(x, n_vec){
  map(n_vec, lag, x = x) %>% 
    set_names(paste('lag', n_vec)) %>% 
    as_tibble()
}
timeSeries <- 
  DJI_workshop %>%
  mutate(observed = if_else(R > lag(R), 1, 0)) %>%
  mutate(lag_multiple(observed, 1:maxMarkovOrder)) %>%
  select(observed, starts_with('lag ')) %>%
  slice_tail(n = -1 - maxMarkovOrder) # omit the first entry for the lag of the observed variable is used
nObservations <- nrow(timeSeries)

for(iOrder in 1:length(markovOrder)){
  
  markovOrderHere <- markovOrder[iOrder]
  codeLengthH <- 1 + 2 * log2(markovOrderHere) + markovOrderHere * log2(nObservations + 1)
  codeLength$`L(H)`[iOrder] <- codeLengthH
  
  # Determine the observed transition matrix.
  transitionMatrix <- 
    GetTransitionMatrix(markovOrderHere)
  if(markovOrderHere == 0){
    transitionMatrix <- 
      transitionMatrix %>%
      mutate(total = nObservations, `to 1 (#)` = sum(timeSeries %>% head(-1) %>% pull(observed)))
  }else{
    
    ixSelect <- (1 + maxMarkovOrder):(nObservations - 1)
    timeSeriesRemaining <- timeSeries %>% select(observed, any_of(paste('lag', 1:markovOrderHere)))
    timeSeriesRemainingMatrix <- timeSeriesRemaining %>% select(-observed) %>% as.matrix()
    
    for(iRoute in 1:nrow(transitionMatrix)){
      
      route <- transitionMatrix %>% slice(iRoute) %>% select(all_of(1:markovOrderHere)) %>% as.numeric()
      timeSeriesRoute <- matrix(route, nrow = 1)[rep(1, nrow(timeSeriesRemaining)),]
      isRoute <- rowSums(timeSeriesRemainingMatrix == timeSeriesRoute) == markovOrderHere
      ixRoute <- which(isRoute)
      transitionMatrix[iRoute, 'total'] <- length(ixRoute)
      transitionMatrix[iRoute, 'to 1 (#)'] <- sum(timeSeriesRemaining$observed[ixRoute])
      #transitionMatrix[iRoute, 'to 0 (#)'] <- length(ixRoute) - sum(timeSeriesRemaining$observed[ixRoute])
      
      # 
      # isRemaining <- rep(TRUE, nrow(timeSeriesRemaining))
      # for(iLag in 1:length(route)){
      #   isRemaining <- isRemaining & timeSeriesRemainingMatrix[,paste('lag', iLag)] == route[iLag]
      #   if(!any(isRemaining)){
      #     next
      #   }
      # }
      # 
      # isUseFilter <- FALSE
      # if(isUseFilter){
      #   for(iLag in 1:length(route)){
      #     timeSeriesRemaining <- 
      #       timeSeriesRemaining %>%
      #       filter(.data[[paste('lag', iLag)]] == route[iLag])
      #     if(nrow(timeSeriesRemaining) == 0){
      #       next
      #     }
      #   }
      #   
      #   transitionMatrix[iRoute, 'total'] <- nrow(timeSeriesRemaining)
      #   transitionMatrix[iRoute, 'to 1 (#)'] <- timeSeriesRemaining %>% pull(observed) %>% sum
      # }
      # 
      # transitionMatrix[iRoute, 'total'] <- sum(isRemaining)
      # transitionMatrix[iRoute, 'to 1 (#)'] <- sum(timeSeries$observed[isRemaining])
    }
  }
  transitionMatrix <- 
    transitionMatrix %>% 
    mutate(`to 0 (#)` = total - `to 1 (#)`,
           `to 0 (%)` = `to 0 (#)` / total,
           `to 1 (%)` = `to 1 (#)` / total)
  maximumLogLikelihood <- 
    sum(transitionMatrix$`to 0 (#)` * log2(transitionMatrix$`to 0 (%)`), na.rm = TRUE) +
    sum(transitionMatrix$`to 1 (#)` * log2(transitionMatrix$`to 1 (%)`), na.rm = TRUE)
  codeLengthDGivenH <- -maximumLogLikelihood
  codeLength$`L(D|H)`[iOrder] <- codeLengthDGivenH
  
  print(transitionMatrix)
  
  print(codeLength %>% slice(1:iOrder))
  
  stopifnot(sum(transitionMatrix$total) == nrow(timeSeries))
  stopifnot(sum(transitionMatrix$`to 0 (#)`) + sum(transitionMatrix$`to 1 (#)`) == nrow(timeSeries))
  
}
codeLength <-
  codeLength %>%
  mutate(`L(H, D)` = `L(H)` + `L(D|H)`)

print(
  ggplotly(
    codeLength %>% 
      ggplot(aes(x = order, y = `L(H, D)`)) + 
      geom_line() + 
      geom_point(data = codeLength %>% slice_min(`L(H, D)`),
                 aes(x = order, y = `L(H, D)`),
                 colour = 'red') + 
      ggtitle('Code length for Markov models of different orders')
  )
)
