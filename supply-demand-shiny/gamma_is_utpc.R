z <- seq(-4,1,0.1)

x <- 1 - z


aa <- 2
th <- 1


y_gamma <- ((x)^(aa-1)) * exp(-x/th)
y_gamma <- y_gamma / max(y_gamma)

y_utpc <- (1-z) * exp(z)


plot(z,y_gamma, type = "l")
points(z,y_utpc, col = "red", lty = 2)
