# Create some arbitrary data
x <- sort(runif(10, 0, 100))
y <- rnorm(10, x^0.5, 0.05 * x^0.5)
plot(x, y) # what do the fake data look like

# interpolate the data
xout <- 1:100
a <- approx(x, y, xout)
lines(a$x, a$y) # this looks OK but...
# it's still distinct points, as can be seen here...
plot(a$x, a$y)

# We want a function, so use approxFun
f <- approxfun(x, y) # this allows us to derive a density for any depth
# this allows us to derive a density for any depth within the range of the data
x.data <- 1:100
y.test <- f(x.data)
plot(x.data, y.test) # this looks exactly the same as the previous plot...
# but the important thing is that our function can be called with any new data...
x.data <- seq(1, 100, 0.1) # for instance, this will now return 10 times greater resolution using the same function
y.test <- f(x.data)
plot(x.data, y.test)


