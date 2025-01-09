require(Yano)

## sSCC.rds is generated at figure04
obj <- readRDS("sSCC.rds")

DefaultAssay(obj) <- "var"

sel <- c("5:170837733CT=/+",  "5:170837733CT>C/+")
dat <- GetAssayData(obj, layer = "data")
dat <- dat[sel,]

x <- dat[1,]
y <- dat[2,]
z <- x+y

W <- obj[["pca_wm"]]
W <- as.sparse(W)

##############################################################################################################
##   Supplementary Figure 1A
##############################################################################################################

## DScore() can be found at Yano/dissimilarity.R

for (n in c(10, 20, 50, 100, 200, 500, 1000, 10000)) {
  s <- unlist(lapply(1:n, function(i) {  x0 <- sample(x); DScore(x0,z,W) }))
  png(file = paste0("perm_",n,".png"))
  par(cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2, mar = c(5, 5, 2, 2))
  br <- ifelse(n < 100, 12, 50)
  h <- hist(s, breaks = br, density = 10, col = "black", xlab = "permutated D scores", main = paste0("n = ",n))
  h
  mn <- mean(s)
  std <- sqrt(var(s))
  xfit <- seq(min(s), max(s), length = 40)
  yfit <- dnorm(xfit, mean = mn, sd = std)
  yfit <- yfit * diff(h$mids[1:2]) * length(s)
  abline(v=mn, lwd=2, col = "red")
  text(x=mn, y = max(h$counts), labels = paste0("mean = ", round(mn, digits = 2)), cex=2,col="red", adj = 0)
  lines(xfit, yfit, col = "blue", lwd = 2)
  dev.off()
}

##############################################################################################################
##   Supplementary Figure 1B
##############################################################################################################

png(file = "QQplot_n10000.png")
par(cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2, mar = c(5, 5, 2, 2))
qqnorm(s)
qqline(s)
dev.off()
