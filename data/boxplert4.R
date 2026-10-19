boxplert <- function(X,
                     Y,
                     main = NULL,
                     xlab = NULL,
                     ylab = NULL,
                     bcol = "bisque",
                     p.adj = "BH",   # 认改成BH（推荐）
                     cexy = 1,
                     varwidth = TRUE,
                     las = 1,
                     paired = FALSE) 
{
  # 因子化
  Y <- as.factor(Y)
  aa <- levels(Y)
  an <- as.character(1:length(aa))
  
  # =========================
  # 1. 描述统计
  # =========================
  tt1 <- data.frame(
    mean   = tapply(X, Y, mean, na.rm = TRUE),
    se     = tapply(X, Y, function(x) sd(x, na.rm = TRUE)/sqrt(length(x))),
    sd     = tapply(X, Y, sd, na.rm = TRUE),
    min    = tapply(X, Y, min, na.rm = TRUE),
    max    = tapply(X, Y, max, na.rm = TRUE),
    median = tapply(X, Y, median, na.rm = TRUE),
    n      = tapply(X, Y, length)
  )
  
  # =========================
  # 2. 作图
  # =========================
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
  
  # =========================
  # 3. 建模
  # =========================
  Yn <- factor(Y, labels = an)
  model <- aov(X ~ Yn)
  
  # =========================
  # 4. 前提检验（重要！）
  # =========================
  shapiro_p  <- tryCatch(shapiro.test(residuals(model))$p.value, error=function(e) NA)
  bartlett_p <- tryCatch(bartlett.test(X ~ Yn)$p.value, error=function(e) NA)
  
  # =========================
  # 5. 显著性检验
  # =========================
  if (paired == TRUE & length(aa) == 2) {
    pp <- t.test(X ~ Yn, paired = TRUE)$p.value
  } else {
    pp <- summary(model)[[1]][["Pr(>F)"]][1]
  }
  
  # 显著性符号
  sig <- "ns"
  if (pp <= 0.1)  sig <- "."
  if (pp <= 0.05) sig <- "*"
  if (pp <= 0.01) sig <- "**"
  if (pp <= 0.001) sig <- "***"
  
  mtext(sig, side = 3, line = 0.5, adj = 0, cex = 2)
  
  # =========================
  # 6. 多重比较（仅显著时）
  # =========================
  tt1$group <- NA
  
  if (pp <= 0.05) {
    if (!requireNamespace("agricolae", quietly = TRUE)) {
      stop("Please install package 'agricolae'")
    }
    
    comp <- agricolae::LSD.test(
      model,
      "Yn",
      alpha = 0.05,
      p.adj = p.adj,
      group = TRUE
    )
    
    gror <- comp$groups[order(rownames(comp$groups)), , drop = FALSE]
    tt1$group <- gror$groups
    
    mtext(
      tt1$group,
      side = 3,
      at = 1:length(aa),
      line = 0.5,
      cex = 1,
      font = 4
    )
  }
  
  # =========================
  # 7. 输出结果
  # =========================
  return(list(
    comparison = tt1,
    p.value = pp,
    shapiro_p = shapiro_p,
    bartlett_p = bartlett_p
  ))
}