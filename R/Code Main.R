pkgs <- c("gsw", "tidyverse", "ncdf4", "lubridate", "dplyr") # required packages
for(i in 1:length(pkgs)){
  # Load packages, installing if necessary
  j <- library(pkgs[i], character.only = TRUE, logical.return = TRUE)
  if(!j){
    install.packages(pkgs[i])
    library(pkgs[i], character.only = TRUE)}}

data.list <- list() # 'list' is a reserved keyword R uses to create lists -- avoid using reserved keywords to define objects
# list <- list()
getwd()
data_files <- list.files(recursive = TRUE)
data_files <- data_files[grepl('.nc', data_files)]
for(i in 1:length(data_files)){
  nc_ds <- nc_open(data_files[i])
  dim_lat <- ncvar_get(nc_ds, "LATITUDE")
  dim_lon <- ncvar_get(nc_ds, "LONGITUDE")
  dim_temp <- ncvar_get(nc_ds, "TEMP")
  dim_psal <- ncvar_get(nc_ds, "PSAL")
  dim_pres <- ncvar_get(nc_ds, "PRES")
  dim_time <- ncvar_get(nc_ds, "TIME")
  
  dim_depth <- gsw_z_from_p(dim_pres, dim_lat) # calculate depth from pressure and latitude
  
  # origin <- as.Date('1950-01-01')
  origin <- '1950-01-01'
  sec <- dim_time * 24 * 60 * 60
  # sec <- dim_time[1] * 24 * 60 * 60
  dim_time <- as.POSIXct(sec, origin = origin) # convert stored 'time' into actual sample time
  dim_year <- as.numeric(strftime(dim_time, '%Y'))
  dim_month <- as.numeric(strftime(dim_time, '%m'))
  dim_day <- as.numeric(strftime(dim_time, '%d'))
  
  j <- which.min(apply(dim_depth, 2, function(z) min(z, na.rm = TRUE)))
  
  # name the variables for clarity, rather than using X1, X2 etc
  tesa <- data.frame(depth = dim_depth[,j], psal = dim_psal[,j], temp = dim_temp[,j], pres = dim_pres[,j])
  # tesa <- data.frame(cbind(dim_psal[,j], dim_temp[,j], dim_pres[,j]))
  tesa
  
  # Can combine new variables within a single call of 'mutate'
  tesa <- tesa %>% mutate(lon = dim_lon[j],
                          lat = dim_lat[j],
                          year = dim_year[j],
                          month = dim_month[j],
                          day = dim_day[j]
  )
  
  
  data.list[[i]] <- tesa
  
  par(mfrow=c(1,2))
  
  
  plot(depth~psal, data=data.list[[i]])
  plot(depth~temp, data=data.list[[i]])
  Title1 <- paste(dim_day[j], month.name[dim_month[j]], dim_year[j])
  mtext(Title1, side = 3, outer = TRUE, line = -1)
  Title2 <- bquote('(' * .(dim_lon[j]) * degree ~ E * ',' ~ .(dim_lat[j]) * degree ~ S * ')')
  mtext(Title2, side = 3, outer = TRUE, line = -2.5)
}

# Include a sample index in each list element
data.list <- lapply(1:length(data.list), function(z){
  x <- data.list[[z]]
  x$sample <- z
  x})


data.list
new.data <- do.call('rbind', data.list)
summary(na.omit(new.data$lat))



new.data$asal <- gsw_SA_from_SP(new.data$psal, p = new.data$pres,
                                longitude = new.data$lon, latitude = new.data$lat)
# Calculate seawater density
new.data$rho <- gsw_rho(SA = new.data$asal, CT = new.data$temp, p = new.data$pres)
new.data$rho <- new.data$rho / 1000 # convert units to g/cm^3
# we can convert back to a list, but it's probably easier to work with the data
# frame, using the sample index column to differentiate samples.
data.list <- lapply(unique(new.data$sample), function(z){
  new.data[new.data$sample == z,]})

par(mfrow = c(2,2))
for(i in unique(new.data$sample)){
  d <- new.data[new.data$sample == i,]
  Title <- paste(d$day[1], month.abb[d$month[1]], d$year[1])
  plot(d$rho, d$depth, xlab = expression(sewater ~ density ~ (g/cm^3)),
       ylab = 'depth (m)', main = Title)
}

# Average the density-at-depth measurements across samples within months. The depths
# may not line up, so first create smooth interpolating functions to find density
# at the same depths in every sample. Then create a smooth interpolating function
# that returns density-at-depth for each month (this is needed for our equations).
new.data$depth <- new.data$depth * 100 # convert units of depth to cm (our interpolating functions are used within the sinking speed equations which use g-cm-s units)
months <- unique(new.data$month)
nmonths <- length(months)
depth.range <- data.frame(min = rep(NA, nmonths), max = rep(NA, nmonths), month = months) # store range of sampled depths for each month
rho.interp <- setNames(vector(mode = 'list', length = nmonths), month.abb[months]) # store interpolating function for each month
for(i in 1:nmonths){
  d <- new.data[new.data$month == months[i],] # subset data for month i
  samples <- unique(d$sample) # sample indices for month i
  nsamples <- length(samples)
  dr <- sapply(samples, function(z){ # depth range of each sample
    x <- d[d$sample == z,]
    range(x$depth, na.rm = TRUE)})
  dr <- c(ceiling(min(dr[2,])), floor(max(dr[1,]))) # a depth range for month i  
  depth.range[depth.range$month == months[i], c('min','max')] <- dr
  depth.inc <- seq(dr[1], dr[2], -10) # uniform depth increments [cm] to feed into interpolating functions
  ndepth <- length(depth.inc)
  rho.matrix <- matrix(NA, ndepth, nsamples) # store interpolated densities for each sample
  for(j in 1:nsamples){
    dj <- d[d$sample == samples[j],] # subset data for sample j in month i
    f <- approxfun(dj$depth, dj$rho)
    rho.matrix[,j] <- f(depth.inc) # interpolate
  }
  rho <- rowMeans(rho.matrix, na.rm = TRUE) # average density (at depths depth.inc) across samples
  rho.interp[[i]] <- approxfun(depth.inc, rho, rule = 2) # create and store a single interpolating function for density-at-depth per month
}
# Now we can use rho.interp to return density at any depth (within the bounds of
# the data) for each month. 
# In terms of ocean properties, all we need is the range of sampled depths, 'depth.range',
# interpolating functions for density, 'rho.interp', and viscosity which we assume as constant.

# Now that we have these values and functions, we no longer need the data tables
# of pressure and temperature etc.

# Set up a numerical integration of the sinking & degradation equations.
# Choose a time step. The value matters -- there is a trade off between accuracy
# and time taken to solve the equations. A small time step will provide more
# accurate answers, but take longer (this can be 'dialed-in' later).
dt <- 10 * 60 # 10 minute time step [s]
# Set parameter values (in correct units)
L.init <- 0.2927 # initial length [cm]
D.init <- 0.0183 # initial diameter [cm]
gd <- 1.121 # faecal pellet density [g/cm^3]
rD <- 1.2 / 86400 # degradation rate [1/s]
vis <- 0.0189 # seawater viscosity []
coef <- 0.0790 # scaling coefficient from Komar
g <- 981 # acceleration due to gravity [cm/s^2]
p <- -1.644 # power term from Komar

# Get the initial conditions...

# I've moved some equations here so I can see what's going a bit better... move them back if you want.
# I've also slightly altered how they're written, but not made any major changes.
volume <- function(L,D){
  pi / 4 * L * D^2}
mass <- function(V, gd){
  V * gd}
shape <- function(L, D){
  L / D}

V.init <- volume(L.init, D.init) # initial volume [cm^3]
m.init <- mass(V.init, gd) # initial mass [g]
S <- shape(L.init, D.init) # shape parameter describing faecal pellets, assumed constant [unitless]

mass.degrade <- function(t, r = rD, m0 = m.init){
  m0 * exp(-r * t)}

speed <- function(m, rhoF, rhoW, visc = vis, 
                  grav = g, pow = p, const = coef, s = S){
  const / visc * {rhoF - rhoW} * grav * s^pow * 
    {{{4 * s^2} / {pi * rhoF}} * m}^{2/3}
}

depth.init <- sapply(1:nmonths, function(z) depth.range$min[z]) # initial depths [cm], one for each month, depend on depth range of samples
rhoW.init <- sapply(1:nmonths, function(z) rho.interp[[z]](depth.init[z])) # initial seawater density

sp.init <- speed(m.init, gd, rhoW.init) # initial speed [cm/s]

# Create list of data frames to store output for each month, and populate with initial values
output <- setNames(lapply(1:nmonths, function(z){
  data.frame(time = 0, depth = depth.init[z], volume = V.init, mass = m.init, speed = sp.init[z])}),
  month.abb[months])

# Integrate with a basic forward Euler method, for each month.
for(i in 1:nmonths){
  stopLoop <- FALSE # TRUE stops the while loop
  t <- 0 # time [s]
  j <- 0 # iteration counter
  while(!stopLoop){
    t <- t + dt # each loop iteration steps forward one time increment
    j <- j + 1
    distance <- output[[i]]$speed[j] * dt # distance [cm] traveled in time dt
    depth <- output[[i]]$depth[j] - distance # new depth
    m <- mass.degrade(t) # new mass
    rho.W <- rho.interp[[i]](depth) # seawater density
    sp <- speed(m, gd, rho.W) # updated speed
    V <- m / gd # new volume -- we don't need this to solve the equations, but it may be useful once plastics are considered...
    output[[i]][j+1,] <- data.frame(time = t, depth = depth, volume = V, mass = m, speed = sp)
    # Conditions for breaking out of the while loop...
    depth.range.exceeded <- output[[i]]$depth[j+1] < depth.range$max[i] # has pellet descended below depth limit of measurements?
    too.slow <- output[[i]]$speed[j+1] * 86400 < 1 # is pellet moving so slow that it's effectively stationary? (I've chosen 1 cm/day, but this can be changed).
    if(depth.range.exceeded | too.slow) stopLoop <- TRUE
    print(paste0(month.name[months[i]], ': j = ', j))
  }
}

# Convert units to something more readable
output <- lapply(1:length(output), function(z){
  x <- output[[z]]
  x$time <- x$time / 86400 # [day]
  x$depth <- x$depth / 100 # [m]
  x$volume <- x$volume * 1000 # [mm^3]
  x$mass <- x$mass * 1000 # [mg]
  x$speed <- x$speed * 864 # [m/day]
  x
})

# View results
i <- 1
head(output[[i]])
tail(output[[i]])

plot(output[[i]]$time, output[[i]]$depth, type = 'l',
     xlab = 'time (days)', ylab = 'depth (m)')
plot(output[[i]]$mass, output[[i]]$depth, type = 'l',
     xlab = expression(faecal ~ pellet ~ mass ~ (mg)), ylab = 'depth (m)')
plot(output[[i]]$speed, output[[i]]$depth, type = 'l',
     xlab = 'sinking speed (m/day)', ylab = 'depth (m)')
plot(output[[i]]$volume, output[[i]]$depth, type = 'l',
     xlab = expression(faecal ~ pellet ~ volume ~ (mm^3)), ylab = 'depth (m)')
mtext(month.name[months[i]], outer = TRUE, line = -1.5)


# MP volume calculations based on assumed diameters (values from Johnston et al, 2023)
# Fiber - L of 240 microns - assume D = 1/100 of L
fib <- volume(0.024, 0.00024)
fib

# Fragment - L of 20 microns - assume D = 1/4 of L
frag <- volume(0.002,0.0005)
frag

# 14,740 MP particles found in FP / 15 krill / day (Dawson et al, 2018)
# particle # is an underestimation, as 14740 particles counted during the 1st day
# of MP-free diet after 10 days of MP-contaminated diet
MpFP <- 14740
NumK <- 15

# particles / krill
PerK <- MpFP/NumK

# Krill egestion occurs ~1/3hrs
PerKEg <- PerK / 8

# Fibers account for ~82% of MP, fragments for ~18% (assume no beads or sheets)
fibK <- 0.82*PerKEg
fragK <- 0.18*PerKEg

# To get total volume of MPs in FP
tVperK <-  (fibK * fib) + (fragK * frag)
tVperK

# % of FP volume that is MPs
rat <- tVperK/V.init
rat * 100

# Polymer type densities (g/cm^3)
PA6 <- 1.14
PE <- 0.95
EVA <- 0.93
PBR <- 0.9

# Approximate ratio in fecal pellet into overall density of MPs
Ovrll <- (0.55*PA6) + (0.18*PE) + (0.17*EVA) + (0.1*PBR)
Ovrll

# Overall density of MP contaminated FP
(gd * (1-rat)) + (Ovrll * rat)

