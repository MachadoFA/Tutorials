## ----plotting-defaults, include=FALSE-----------------------------------------
options(width=85)
blue <- "#206A9B"; orange <- "#BE541E"; gray <- "#707780"
par_default <- function() par(mar=c(4.3,4.4,2.4,1), las=1, bty="l")


## ----projection, fig.cap="The dashed segment meets the chosen axis at a right angle. The projected point has two coordinates, but its position along the axis is summarized by one score: s = 2."----
par_default()
y <- c(2,1); u <- c(3,4)/5
s <- sum(y*u); projected <- s*u
plot(NA, xlim=c(-0.5,2.8), ylim=c(-0.5,2.8), asp=1,
     xlab="Shape coordinate 1", ylab="Shape coordinate 2")
abline(h=0,v=0,col="gray90")
segments(-0.45*u[1],-0.45*u[2],3.2*u[1],3.2*u[2],col=blue,lwd=2)
arrows(0,0,u[1],u[2],length=.09,col=blue,lwd=3)
segments(y[1],y[2],projected[1],projected[2],lty=2,col=gray)
points(y[1],y[2],pch=19,cex=1.4,col=orange)
points(projected[1],projected[2],pch=15,cex=1.2,col=blue)
text(y[1],y[2],"Specimen (2, 1)",pos=4,cex=.82)
text(projected[1],projected[2],"Projection (1.2, 1.6)",pos=2,cex=.82)
text(.22,.67,"Unit vector u",pos=2,cex=.82,col=blue)
points(0,0,pch=3); text(0,0,"Mean",pos=2,cex=.8)


## ----simulate-----------------------------------------------------------------
set.seed(42)
n <- 120
x <- seq(-2,2,length.out=n)
beta_true <- c(0.025,0.055)
u_true <- beta_true/sqrt(sum(beta_true^2))
u_perp <- c(-u_true[2],u_true[1])

# Predictor effect + noise along that direction + larger perpendicular noise.
Y <- outer(x,beta_true) +
     outer(rnorm(n,sd=0.028),u_true) +
     outer(rnorm(n,sd=0.14),u_perp)
colnames(Y) <- c("shape_1","shape_2")
Yc <- scale(Y,center=TRUE,scale=FALSE)
xc <- x-mean(x)

# Do not standardize shape coordinates to unit variance.
pca <- prcomp(Y,center=TRUE,scale.=FALSE)
v1 <- pca$rotation[,1]
pc1_score <- drop(Yc %*% v1)
stopifnot(isTRUE(all.equal(pc1_score,unname(pca$x[,1]))))


## ----pca, fig.cap="Left: the shape cloud and the PC1 axis. Right: the same specimens reduced to their signed positions along PC1. The axis is chosen without using x."----
par(mfrow=c(1,2),mar=c(4.3,4.4,2.4,1),las=1,bty="l")
lim <- max(abs(Yc))*1.16
L <- seq(-lim,lim,length.out=100)
plot(Yc,asp=1,xlim=c(-lim,lim),ylim=c(-lim,lim),pch=19,
     col=adjustcolor(gray,.6),xlab="Shape coordinate 1",
     ylab="Shape coordinate 2",main="PC1 in shape space")
lines(L*v1[1],L*v1[2],col=blue,lwd=2)
ii <- c(10,30,60,90,110)
PP <- outer(pc1_score,v1)
segments(Yc[ii,1],Yc[ii,2],PP[ii,1],PP[ii,2],col=blue,lty=2)
points(PP[ii,1],PP[ii,2],pch=15,col=blue)
plot(pc1_score,rep(0,n),pch=19,col=adjustcolor(blue,.55),
     xlab="PC1 score",ylab="",yaxt="n",ylim=c(-1,1),
     main="Positions along PC1")
abline(v=0,lty=3,col=gray)


## ----regression---------------------------------------------------------------
fit <- lm(Y ~ x)
b <- as.numeric(coef(fit)["x",])
b_length <- sqrt(sum(b^2))
stopifnot(b_length > 0)
u_reg <- b/b_length
b


## ----slopes, fig.cap="Each regression has a slope. The pair of slopes gives the displacement in shape space associated with one unit of increase in x."----
par(mfrow=c(1,2),mar=c(4.3,4.4,2.4,1),las=1,bty="l")
for(j in 1:2){
  plot(x,Y[,j],pch=19,col=adjustcolor(gray,.6),xlab="Predictor x",
       ylab=paste("Shape coordinate",j),
       main=sprintf("Slope b%d = %.4f",j,b[j]))
  abline(coef(fit)[1,j],b[j],col=orange,lwd=2)
}


## ----axes, fig.cap="The same shape space, with two different axes. Color marks the predictor from low (blue) to high (orange). PC1 follows the dominant variation; the regression axis follows the estimated change with x."----
par_default()
cols <- colorRampPalette(c(blue,"#D4D4D4",orange))(n)
plot(Yc,asp=1,xlim=c(-lim,lim),ylim=c(-lim,lim),pch=19,col=cols,
     xlab="Shape coordinate 1",ylab="Shape coordinate 2")
lines(L*v1[1],L*v1[2],col=blue,lwd=2,lty=2)
lines(L*u_reg[1],L*u_reg[2],col=orange,lwd=2)
arrows(0,0,.20*u_reg[1],.20*u_reg[2],col=orange,lwd=2,length=.09)
legend("topleft",c("PC1 axis","Regression axis: increasing x"),
       col=c(blue,orange),lty=c(2,1),lwd=2,bty="n",cex=.85)


## ----scores-------------------------------------------------------------------
shape_score <- drop(Yc %*% u_reg)
projection <- outer(shape_score,u_reg)
fitted_centered <- outer(xc,b)
fitted_score <- xc*b_length


## ----score-plot, fig.cap="Left: observed shapes projected orthogonally onto the regression axis. Right: their scores plotted against x. The line is the score predicted by the multivariate regression, not a separate regression fitted to independent data."----
par(mfrow=c(1,2),mar=c(4.3,4.4,2.4,1),las=1,bty="l")
plot(Yc,asp=1,xlim=c(-lim,lim),ylim=c(-lim,lim),pch=19,col=cols,
     xlab="Shape coordinate 1",ylab="Shape coordinate 2",
     main="Project onto the axis")
lines(L*u_reg[1],L*u_reg[2],col=orange,lwd=2)
segments(Yc[ii,1],Yc[ii,2],projection[ii,1],projection[ii,2],
         lty=2,col=gray)
points(projection[ii,1],projection[ii,2],pch=15,col=orange)
plot(x,shape_score,pch=19,col=cols,xlab="Predictor x",
     ylab="Shape score",main="Plot scores against x")
lines(x,fitted_score,col=orange,lwd=2)


## ----observed-vs-fitted, fig.cap="Observed shapes give scores with scatter. Fitted shapes give scores exactly on a line. Use observed scores to visualize how specimens vary around the fitted relationship."----
par(mfrow=c(1,2),mar=c(4.3,4.4,2.4,1),las=1,bty="l")
yr <- range(c(shape_score,fitted_score))
plot(x,shape_score,pch=19,col=cols,ylim=yr,xlab="Predictor x",
     ylab="Observed shape score",main="Observed shapes")
lines(x,fitted_score,col=orange,lwd=2)
plot(x,fitted_score,pch=19,col=orange,ylim=yr,xlab="Predictor x",
     ylab="Fitted shape score",main="Fitted shapes only")


## ----variance-----------------------------------------------------------------
E <- residuals(fit)
R2_shape <- 1-sum(E^2)/sum(Yc^2)
R2_score <- summary(lm(shape_score ~ x))$r.squared
round(c(total_shape_R2=R2_shape,score_plot_R2=R2_score),3)


## ----reusable-function--------------------------------------------------------
regression_shape_scores <- function(Y,x) {
  Y <- as.matrix(Y)
  if(!is.numeric(Y) || !is.numeric(x) || nrow(Y)!=length(x))
    stop("Use numeric shape coordinates and a matching numeric predictor.")
  if(anyNA(Y) || anyNA(x) || any(!is.finite(Y)) || any(!is.finite(x)))
    stop("Align specimens and resolve missing/nonfinite values before fitting.")
  if(nrow(Y)<3 || ncol(Y)<1 || var(x)==0)
    stop("Need at least three observations, shape coordinates, and varying x.")
  Yc <- scale(Y,center=TRUE,scale=FALSE)
  fit <- lm(Y ~ x)
  b <- as.numeric(coef(fit)["x",])
  length_b <- sqrt(sum(b^2))
  if(!is.finite(length_b) || length_b==0)
    stop("A zero regression vector has no direction to project onto.")
  u <- b/length_b
  list(score=drop(Yc %*% u),direction=u,beta=b,
       fitted_score=(x-mean(x))*length_b,
       mean_shape=colMeans(Y),fit=fit,
       R2_shape=1-sum(residuals(fit)^2)/sum(Yc^2))
}
result <- regression_shape_scores(Y,x)
stopifnot(isTRUE(all.equal(result$score,shape_score)))


## ----full-pc-check------------------------------------------------------------
pc_fit <- lm(pca$x ~ x)
b_pc <- as.numeric(coef(pc_fit)["x",])
s_pc <- drop(pca$x %*% (b_pc/sqrt(sum(b_pc^2))))
stopifnot(max(abs(s_pc-shape_score)) < 1e-10)


## ----checks, include=FALSE----------------------------------------------------
stopifnot(abs(sum(u_reg^2)-1)<1e-12,
          max(abs(drop(crossprod(xc,E))))<1e-10,
          max(abs(shape_score-(fitted_score+drop(E %*% u_reg))))<1e-10,
          abs(coef(lm(shape_score ~ x))[2]-b_length)<1e-10,
          R2_score >= R2_shape-1e-12)

