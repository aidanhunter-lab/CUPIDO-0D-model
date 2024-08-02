# An example plot for visualising the theory and a useful output metric.

# Faecal pellets degrade through bacterial action as they sink. It would be good
# to determine the proportion of the faecal matter remaining as a function of
# depth. Make a plot to illustrate...

# Parameters -- these are dummy values used only for the plot.
r <- 0.25 # degradation rate [1/day] (assume constant across all depths)
v <- 200 # sink speed [m/day]
m0 <- 100 # initial mass of faecal pellet [ug] (sometimes 'u' is used to represent micro)

depth_min <- 0 # surface depth [m] 
depth_max <- 5000 # (arbitrary) maximum depth [m]
dd <- 1 # depth increment [m]
depth <- seq(depth_min, depth_max, dd)

t0 <- 0 # initial time (day)
tf <- 10 # final time (day)
dt <- 1/24 # time increment (day)
time <- seq(t0, tf, dt)
nt <- length(time) # number of discrete times

m <- m0 * exp(- r * time) 

plot(time, m / m0)

plot(m / m0, time)

plot(rev(m / m0), time)







