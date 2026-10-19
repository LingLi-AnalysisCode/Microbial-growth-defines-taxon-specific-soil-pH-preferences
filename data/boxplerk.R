##非参数组间差异的Kruskal-Wallis检验
boxplerk <-  function(X,
                      Y,
                      main = NULL,
                      xlab = NULL,
                      ylab = NULL,
                      bcol = "bisque",
                      p.adj = "none",
                      cexy = 1,
                      varwidth = TRUE,
                      las = 1,
                      paired = FALSE)
{
  aa <- levels(as.factor(Y))
  an <- as.character(c(1:length(aa)))
  tt1 <- matrix(nrow = length(aa), ncol = 7)    
  for (i in 1:length(aa))
  {
    temp <- X[Y == aa[i]]
    tt1[i, 1] <- mean(temp, na.rm = TRUE)
    tt1[i, 2] <- sd(temp, na.rm = TRUE) / sqrt(length(temp))
    tt1[i, 3] <- sd(temp, na.rm = TRUE)
    tt1[i, 4] <- min(temp, na.rm = TRUE)
    tt1[i, 5] <- max(temp, na.rm = TRUE)
    tt1[i, 6] <- median(temp, na.rm = TRUE)
    tt1[i, 7] <- length(temp)
  }
  
  tt1 <- as.data.frame(tt1)
  row.names(tt1) <- aa
  colnames(tt1) <- c("mean", "se", "sd", "min", "max", "median", "n")
  
  boxplot(
    X ~ Y,
    main = main,
    xlab = xlab,
    ylab = ylab,
    las = las,
    col = bcol,
    cex.axis = cexy,
    cex.lab = cexy,
    varwidth = varwidth
  )    
  require(agricolae)
  Yn <- factor(Y, labels = an)
  comp <- kruskal(X, Yn, p.adj = p.adj)
  sig <- "ns"
  
  if (paired == TRUE & length(aa) == 2)
  {
    coms <- wilcox.test(X ~ Yn, paired = TRUE)
    pp <- coms$p.value
  }    else
  {
    pp <- comp$statistics$p.chisq
  }    
  if(pp <= 0.1)
    sig <- "."
  if(pp <= 0.05)
    sig <- "*"
  if(pp <= 0.01)
    sig <- "**"
  if(pp <= 0.001)
    sig <- "***"
  
  gror <- comp$groups[order(rownames(comp$groups)), ]
  tt1$rank <- gror$X
  tt1$group <- gror$groups
  mtext(
    sig,
    side = 3,
    line = 0.5,
    adj = 0,
    cex = 2,
    font = 1
  )   
  if(pp <= 0.1)
    mtext(
      tt1$group,
      side = 3,
      at = c(1:length(aa)),
      line = 0.5,
      cex = 1,
      font = 4
    )
  
  list(comparison = tt1, p.value = pp)
  
}