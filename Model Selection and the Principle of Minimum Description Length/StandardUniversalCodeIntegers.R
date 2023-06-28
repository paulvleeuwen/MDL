StandardUniversalCodeIntegers <- function(n){
  logPrevious <- log2(n)
  codeLength <- 0
  while(logPrevious > 0){
    codeLength <- codeLength + logPrevious
    logPrevious <- log2(logPrevious)
  }
  codeLength <- codeLength + log2(2.865)
  stopifnot(codeLength <= 2 * log2(n) + 2)
  return(codeLength)
}