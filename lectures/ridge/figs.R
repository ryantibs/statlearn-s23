sim = function(n, d, rfun=rnorm, ...) {
  X = matrix(rfun(n*d, ...), n, d); X = scale(X)
  sample(eigen(crossprod(X)/n, only.values=TRUE)$val)
}
    
mp_dens = function(s, gam, a, b) {
  1/(2*pi*gam*s) * sqrt((b-s)*(s-a))
}
         
n = 2000; d = 1000
gam = d/n
a = (1-sqrt(gam))^2
b = (1+sqrt(gam))^2
svals = seq(a, b, length=100)
evals1 = sim(n, d, rnorm)
evals2 = sim(n, d, rbinom, size=1, prob=0.5)

mar = c(4,4,2,1)
h = 4.5; w = 9.5
ylim = c(0, 0.9)

pdf(file="rmt.pdf", height=h, width=w)
par(mfrow=c(1, 2))
hist(evals1, breaks=30, col="pink", prob=TRUE,
     xlab="Eigenvalues", ylab="Density", ylim=ylim,
     main=sprintf("Gaussian: n = %i, d = %i", n, d)) 
lines(svals, mp_dens(svals, d/n, a, b), lwd=2)

hist(evals2, breaks=30, col="lightblue", prob=TRUE,
     xlab="Eigenvalues", ylab="Density", ylim=ylim,
     main=sprintf("Bernoulli: n = %i, d = %i", n, d)) 
lines(svals, mp_dens(svals, d/n, a, b), lwd=2)
graphics.off()
