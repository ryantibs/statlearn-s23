sim_once = function(n, d, rfun=rnorm, ...) {
  X = matrix(rfun(n*d, ...), n, d); X = scale(X)
  sample(eigen(crossprod(X)/n, only.values=TRUE)$val)
}

mp_dens = function(s, gam, a, b) {
  1/(2*pi*gam*s) * sqrt((b-s)*(s-a))
}

mp_dist = function(s, gam, a, b, tol=.Machine$double.eps^0.5) {
  integrate(mp_dens, lower=a, upper=s, abs.tol=tol, gam, a, b)$val
}

sim_loop = function(gam=0.5, rfun=rnorm, ..., rname=NULL, 
                    nmin=100, nmax=2000, nn=100,
                    t0=0.1, nbreaks=30, col="pink") {
  
  nvals = round(exp(seq(log(nmin), log(nmax), length=nn)))
  dvals = round(gam*nvals)
  errs = rep(NA,nn)
  
  a = (1-sqrt(gam))^2
  b = (1+sqrt(gam))^2
  svals = seq(a, b, length=100)
  f0 = sapply(svals, mp_dist, gam, a, b)
  par(mfrow=c(1,2), mar=c(4.5,4.5,2,2))
  pre = ifelse(is.null(rname), "", paste0(rname, ": "))
  
  for (i in 1:nn) {
    t1 = system.time({
      evals = sim_once(nvals[i], dvals[i], rfun, ...)
    })[3]
    f1 = sapply(svals, ecdf(evals))
    errs[i] = max(abs(f1-f0))
    beta = coef(lm(log10(errs) ~ log10(nvals)))
    
    hist(evals, breaks=nbreaks, col=col, prob=TRUE,
         xlab="Eigenvalues", ylab="Density", 
         main=sprintf("%sn = %i, d = %i", pre, nvals[i], dvals[i]))  
    lines(svals, mp_dens(svals, gam, a, b), lwd=2)
    
    plot(nvals, errs, log="xy", type="o", pch=20,
         xlab="n", ylab="sup_x |F_n(x) - F(x)|",
         main="Kolmogorov Distance")
    if (all(!is.na(beta))) {
      abline(coef=beta, col="blue", lty=2, lwd=2)
      legend("topright", col="blue", lty=2, lwd=2,
             legend=sprintf("Slope = %0.3f",beta[2]))
    }
    if (t1 <= t0) Sys.sleep(t0-t1)
  }
}

# Gaussian
sim_loop(gam=0.5, rfun=rnorm, rname="Gaussian")

# Binomial
sim_loop(gam=0.5, rfun=rbinom, size=1, prob=0.5, rname="Binomial")

# t, df=9
sim_loop(gam=0.5, rfun=rt, df=9, rname="t dist, df = 9")

# t, df=1
sim_loop(gam=0.5, rfun=rt, df=1, rname="t dist, df = 1")
