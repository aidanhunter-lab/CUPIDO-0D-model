
volume_fun <- function(L, W = NULL, A = NULL, plastic.type = NULL,
                       fibre.width.default = 30){
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Calculate volume of plastic pieces given...
  # L = length [required]
  # W = width [wanted but not strictly necessary]
  # A = area (in length-width cross section) [wanted but not strictly necessary]
  # plastic.type = 'fragment' or 'fibre'
  # fibre.width.default = assumed width (μm) of fibres if not specified by W
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if(is.null(plastic.type)){
    plastic.type <- 'fragment'
    warning("Argument 'plastic.type' was not specified: assuming plastic.type = 'fragment'")}
  
  # Get plastic shape from plastic.type and the available dimensions
  switch(plastic.type,
         fibre = {
           shape <- 'cylinder'
           if(is.null(W)){
             W <- fibre.width.default
             warning(paste0("Fibre width, 'W', was not included as an argument so using the default W = ",
                            fibre.width.default, " μm. Note the units of this default and change argument 'fibre.width.default' if needed."))
           }
         },
         fragment = {
           if(is.null(W)) shape <- 'sphere' else{
             if(is.null(A)) shape <- 'ellipsoid' else shape <- 'ragged'}
         }
  )
  
  if(shape == 'ragged' & A > L * W) warning("Measured area, 'A', exceeds theoretical (rectangular) maximum, 'L x W'!")
  
  # Calculate volume from dimensions and shape
  V <- switch(shape,
              cylinder = pi / 4 * L * W^2,
              sphere = pi / 6 * L^3,
              ellipsoid = {
                r <- W / L
                H <- r * W # estimate height from ratio of width to length
                pi / 6 * L * W * H
              },
              ragged = {
                r <- W / L
                H <- r * W # estimate height from ratio of width to length
                R <- L * W # area assuming fragment is rectangular
                P <- A / R # proportion of rectangle filled by fragment
                P * H * A # volume given by height*area, assuming same proportion, P, is filled in the height dimension
              }
  )
  return(V)
}


