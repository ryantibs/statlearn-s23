### Curse of dimensionality ###

d = 1:10
x = seq(min(d),max(d),length=100)
eps = 1e-1

mar = c(4.5,4.5,1,1)
h = 5; w = 6

pdf(file="curse.pdf", height=h, width=w)
par(mar=mar)
plot(x, eps^(-(2+x)/2), type="l",
     xlab="d", ylab="eps^(-(2+d)/d)")
points(d, eps^(-(2+d)/2), pch=19, col=2)
graphics.off()

### Higher-order kernel ###

K = function(x) return(3/8*(3-5*x^2)*(abs(x)<=1))
x = seq(-1.5,1.5,length=100)

mar = c(4.5,4.5,1,1)
h = 5; w = 6

pdf(file="order4.pdf", height=h, width=w)
par(mar=mar)
plot(x, K(x), type="l", xlab="t", ylab="K(t)")
graphics.off()
