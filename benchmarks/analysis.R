load('results.RData')
oldlevels <- levels(results$progname)
for (i in 1:length(oldlevels)){
  if (grepl('devel',oldlevels[i])) 
    oldlevels[i] <- 'R-devel'
  else if (grepl('R-Array-Hash',oldlevels[i])) 
    oldlevels[i] <- 'R-Array-Hash'
  else if (grepl('Revo',oldlevels[i])) 
    oldlevels[i] <- 'Revo-R'
  else if (grepl('vertica',oldlevels[i])) 
    oldlevels[i] <- 'HP-R'
}
levels(results$progname) <- oldlevels

results2 <- data.frame()
for (i in levels(results$progname)){
for (j in levels(results$progname)){
for (k in unique(results$hashsize)){
for (l in unique(results$datafile)){
  # Don't compare progname with itself
  if (i==j) next

  x1 <- subset(results,progname==i & hashsize==k & datafile==l)$runtime
  x1_mu <- mean(x1)
  x2 <- subset(results,progname==j & hashsize==k & datafile==l)$runtime
  x2_mu <- mean(x2)

  # Only calculate speed up, not slow down
  if (x1_mu >= x2_mu) next

  # Coefficient of variance
  x1_cov <- (1-1/length(x1)) * (sd(x1)/x1_mu)
  x2_cov <- (1-1/length(x2)) * (sd(x2)/x2_mu)

  wt <- wilcox.test(x1,x2)
  se<-sqrt(var(log(x2))+var(log(x1)))/sqrt(length(x2))
  
  ci <- 1-exp(mean(log(x1))-mean(log(x2)) + c(-1,1)*1.96*se)
  
  speedup <-  1-exp(mean(log(x1))-mean(log(x2)))
  
  results2 <- rbind(results2,
                data.frame(
                  prog1=i,prog2=j,hashsize=k,datafile=l,x1_mu=x1_mu,x2_mu=x2_mu,
                  x1_cov = x1_cov, x2_cov = x2_cov,
                  p.value = wt$p.value, cilo=ci[1], cihi=ci[2], speedup=speedup)
              )

} } } }
