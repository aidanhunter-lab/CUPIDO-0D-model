library(ggplot2)


# Functions ---------------------------------------------------------------
# Describe faecal pellet decay rate declining with depth by using a Martin curve.
rD_fun <- function(a,b,z) a * z ^ {-b} # Martin curve (is just an exponential decay with depth)


# Parameters --------------------------------------------------------------
# We know decay rates at two specific depths from Morata
r13 <- 7 # degradation rate at chlorophyll max (13 m)
r90 <- 2 # degradation rate at 90 m

# Convert these values into the Martin curve parameters (see algebra below)
pow <- log(r90 / r13) / log(13 / 90) # power term
coef <- r13 * 13 ^ pow # coefficient multiplier
coef <- r90 * 90 ^ pow # or, equivalently

# Algebra -----------------------------------------------------------------
# rD = a * z ^ -b # Martin curve
# log(rD) = log(a * z ^ -b) = log(a) - b * log(z) # take logarithms to get 'b' down from the exponent
# log(rD_13) = log(a) - b * log(13) and log(rD_90) = log(a) - b * log(90) # substitute two known values (from Morata) of rD into the logged Martin curve to produce two equations
# log(a) = log(rD_13) + b * log(13) # use the 1st of these equations to get an expression for log(a)
# log(rD_90) = log(rD_13) + b * log(13) - b * log(90) # substitute the expressio for log(a) into the 2nd equation to get an expression in terms of b
# b = (log(rD_90) - log(rD_13)) / (log(13) - log(90)) = (log(rD_90 / rD_13) / log(13 / 90) # solve for b
# log(rD_13) = log(a) - b * log(13) => a = rD_13 * 13 ^ b # substitute b into one of the orginla equations to solve for a



# Evaluate and plot -------------------------------------------------------
# Define some depths (z > 0, curve is undefined for z = 0)
depths <- seq(1, 1000, length.out = 500)

rD <- rD_fun(coef,pow,depths) # evaluate the parameterised Martin curve at these depths to estimate degradation rate at depth

# view output
dat <- data.frame(depth = depths, rD = rD)
head(dat)

ggplot(data = dat) + 
  geom_point(aes(x = rD, y = depths)) + 
  xlab(expression(degradation ~ rate* ',' ~ r[D] ~ (day^{-1}))) + 
  ylab(expression(depth ~ (m))) + 
  scale_y_continuous(trans = 'reverse')






