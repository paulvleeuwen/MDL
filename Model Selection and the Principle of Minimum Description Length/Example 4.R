rm(list = ls())

library(openxlsx)
library(tidyverse)
library(lubridate)
library(msm)

dateFormat <- '%m/%d/%Y'

# Retrieved from https://www.spglobal.com/spdji/en/web-data-downloads/reports/dja-performance-report-daily.xls?force_download=true.
dataRawFileName <- 'C:/code/MDL/Model Selection and the Principle of Minimum Description Length/dja-performance-report-daily.xlsx'
dataRaw <- 
  read.xlsx(dataRawFileName, startRow = 6, sep.names = " ") %>% 
  tibble()
stopifnot(
  all( 
    as.character( as.Date(dataRaw$`Effective Date`, dateFormat), dateFormat ) == dataRaw$`Effective Date`, 
    na.rm = TRUE)
)

dataCleaned <-
  dataRaw %>%
  select(`Effective Date`, `Close Value`) %>%
  na.omit() %>%
  mutate(`Effective Date` = as.Date(`Effective Date`, dateFormat),
         P = log(`Close Value`),
         R = P - lag(P),
         V = NA)

# The paper mentions 6430 trading days while the number of trading days here is found to be 6535 (in line with the number trading days found
# in https://www.timeanddate.com/date/workdays.html?d1=01&m1=07&y1=1962&d2=31&m2=05&y2=1988&ti=on&).
dataSelection <-
  dataCleaned %>%
  mutate(
    `R indicator` = ifelse(R > lag(R), 1, 0),
    `V indicator` = ifelse(V > lag(V), 1, 0)
  ) %>%
  filter(as.Date('07/01/1962', dateFormat) <= `Effective Date` & `Effective Date` < as.Date('07/01/1988', dateFormat))
dataSelection$V[1] <- var(dataSelection$P)
for(iRow in 2:nrow(dataSelection)){
  dataSelection$V[iRow] <- 0.9 * dataSelection$V[iRow - 1] + 0.1 * dataSelection$R[iRow]^2
}
dataSelection <-
  dataSelection %>%
    mutate(
      `R indicator` = ifelse(R > lag(R), 1, 0),
      `V indicator` = ifelse(V > lag(V), 1, 0)
    ) %>%
    slice_tail(n = -1)

# Reproduce Figure 2.
print(
  dataSelection %>% 
    mutate(`row number` = row_number(), R_diff = R - lag(R)) %>% 
    filter(5400 < `row number`) %>% 
    ggplot(aes(x = `row number`, y = R_diff)) + 
    geom_line()
)
print(
  dataSelection %>% 
    mutate(`row number` = row_number(), V_diff = V - lag(V)) %>% 
    filter(5400 < `row number`) %>% 
    ggplot(aes(x = `row number`, y = pmin(0.004, V_diff))) + 
    geom_line()
)
  
# Calculate the number of bits using Equation (7).
nObservations <- nrow(dataSelection)
R_k <- sum(dataSelection$`R indicator`)
nBitsR <- log2(nObservations) - R_k * log2(R_k / nObservations) - (nObservations - R_k) * log2(1 - R_k / nObservations)
V_k <- sum(dataSelection$`V indicator`)
nBitsV <- log2(nObservations) - V_k * log2(V_k / nObservations) - (nObservations - V_k) * log2(1 - V_k / nObservations)

# Calculate the number of bits using a first-order Markov model.
transitionMatrixNumbers_R <- statetable.msm(`R indicator`, data = dataSelection)
toState_R <- matrix( rowSums(transitionMatrixNumbers_R), ncol = 1 )
transitionMatrix_R <- transitionMatrixNumbers_R / rowSums(transitionMatrixNumbers_R)
# The total number of bits is (i) the number of bits to encode the transition probabilities (minus the number of columns because they need
# to sum up to 1), which is the same for each transition probability plus (ii) the sum of the transition probabilities times the number of
# occurrences.
nBitsRMarkov <- 2 * log2(nObservations) - sum( transitionMatrixNumbers_R * log2(transitionMatrix_R) )
writeLines( sprintf( 'The Markov model for the return series yields an improvement of %d bits.', floor(nBitsR - nBitsRMarkov) ) )

transitionMatrixNumbers_V <- statetable.msm(`V indicator`, data = dataSelection)
toState_V <- matrix( rowSums(transitionMatrixNumbers_V), ncol = 1 )
transitionMatrix_V <- transitionMatrixNumbers_V / rowSums(transitionMatrixNumbers_V)
nBitsVMarkov <- 2 * log2(nObservations) - sum( transitionMatrixNumbers_V * log2(transitionMatrix_V) )
writeLines( sprintf( 'The Markov model for the volatility series yields an improvement of %d bits.', floor(nBitsV - nBitsVMarkov) ) )