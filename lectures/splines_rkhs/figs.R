### B-spline basis ###

Seq = function(a,b) {
  if (a<=b) return(a:b)
  else return(integer(0))
}

Id = function(n) {
  return(bandSparse(n, k=0, diag=list(rep(1,n))))
}

library(Matrix)
Dk = function(n, k, x=1:n/n) {
  I = Id(n)
  D = bandSparse(n, k=c(-1,0), diag=list(rep(-1,n-1), rep(1,n)))
  B = I
  for (i in Seq(1,k)) {
    wts = c(rep(1,i), i/(x[(i+1):n]-x[1:(n-i)]))
    M = bdiag(I[Seq(1,i-1),Seq(1,i-1)], D[1:(n-i+1),1:(n-i+1)])
    B = (wts * M) %*% B
  }
  return(B[-(1:k),])
}

tpf = function(z, k) {
  return(function(x) (x-z)^k * (x>z))
}

bs = function(z) {
  k = length(z) - 2
  D = Dk(k+2,k+1,z)
  funs = sapply(z, function(z0,k) return(tpf(z0,k)), k)
  return(Vectorize(function(x0) 
    return((-1)^(k+1) * (z[k+2]-z[1]) * (1/factorial(k+1)) *
             as.numeric(D %*% sapply(funs, function(f) f(x0))))))
}

n = 16
x = 1:n/(n+1)
kvec = 0:3
n0 = 1000
x0 = seq(0,1,length=n0)
Mmat = matrix(0,n0,length(kvec))

mar = c(3,3,2,1)
h = 4.5; w = 5.5
by = 3
i0 = c(7,5,4,2)

for (i in 1:length(kvec)) {
  k = kvec[i]
  z = x[seq(i0[i],by=by,length=k+2)]
  M = bs(z)
  Mmat[,i] = M(x0)
}

for (i in 1:length(kvec)) {
  k = kvec[i]
  z = x[seq(i0[i],by=by,length=k+2)]
  pdf(file=sprintf("bs%i.pdf",k), height=h, width=w)
  par(mar=mar)
  plot(x0, Mmat[,i], type="l", xlab="", ylab="",
       main=paste("Degree",k))
  abline(v=z, lwd=0.75, col=4, lty=2)
  graphics.off()
}

### Demmler-Reinsch basis ###

n = 50
x = 1:n/n
S = matrix(0,n,n)
lam = rep(0,n)
for (i in 1:n) {
  y = rep(0,n); y[i] = 1
  out = smooth.spline(x,y,df=10)
  S[,i] = out$y
  lam[i] = out$lambda
}

library(MASS)
lambda = mean(lam)
K = (ginv(S) - diag(n))/lambda
s = eigen(K)

mar = c(4.5,4.5,1,1)
h = 5; w = 5.5

m = 8
lgrid = c(1e-5,5e-5,1e-4,5e-4,1e-3,5e-3,0.01,0.05)
evals = matrix(0,n,length(lgrid))
for (i in 1:length(lgrid)) {
  evals[,i] = 1/(1+lgrid[i]*pmax(Re(s$values),0))
}

colorfun = colorRampPalette(c("purple4", "lightblue1"))

pdf(file="reinsch_basis.pdf", height=h, width=w)
par(mar=mar)
matplot(x, s$vectors[,n:(n-m+1)], type="l", lty=1,
        col=colorfun(m), xlab="x", ylab="Eigenvector")
graphics.off()

pdf(file="reinsch_weight.pdf", height=h, width=w)
par(mar=mar)
matplot(1:n, evals, type="l", lty=1,
        col=colorfun(length(lgrid)),
        xlab="Number", ylab="Weight")
graphics.off()

### Equivalent kernel ###

n = 500
x = 1:n/n
S = matrix(0,n,n)
lam = rep(0,n)
for (i in 1:n) {
  y = rep(0,n); y[i] = 1
  out = smooth.spline(x,y,df=10)
  S[,i] = out$y
  lam[i] = out$lambda
}

mar = c(2.5,2.5,2.25,1)
h = 5; w = 5.5
rownums = round(c(n/5,n/2,4*n/5))

pdf(file="ss_kernel1.pdf", height=h, width=w)
par(mar=mar)
matplot(x, t(S[rownums,]), type="o", pch=20, cex=0.75, lty=1,
        main="Inputs points on a grid", xlab="", ylab="")
legend("topright", lty=1, col=1:3, legend=paste("Row",rownums))
graphics.off()

set.seed(1)
x = sort(rnorm(n, mean=0.5, sd=0.1))
S = matrix(0,n,n)
lam = rep(0,n)
for (i in 1:n) {
  y = rep(0,n); y[i] = 1
  out = smooth.spline(x,y,df=10)
  S[,i] = out$y
  lam[i] = out$lambda
}

rownums = round(c(n/5,n/2,4*n/5))

pdf(file="ss_kernel2.pdf", height=h, width=w)
par(mar=mar)

matplot(x, t(S[rownums,]), type="o", pch=20, cex=0.75, lty=1,
        main="Inputs points ~ N(0.5,0.01)", xlab="", ylab="")
legend("topright", lty=1, col=1:3, legend=paste("Row",rownums))
graphics.off()
