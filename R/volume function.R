
volume_fun <- function(L, W = NULL, A = NULL, shape = NULL, plastic.type = NULL,
                       fibre.width.default = 30, suppress.warnings = FALSE){
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Calculate volume of objects assumed to be ellipsoidal or cylindrical.
  # L = length [required]
  # W = width [wanted but not strictly necessary]
  # A = area (in length-width cross section) [wanted for plastic fragments but not strictly necessary]
  # shape = 'cylinder' or 'ellipsoid' -- may be left NULL only if plastic.type is specified
  # plastic.type = 'fragment' or 'fibre' -- this may be left as NULL provided shape is specified, or set to coerce the shape.
  # fibre.width.default = assumed width (μm) of fibres if not specified by W
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if(is.null(shape) && is.null(plastic.type)){
    warning("'shape' or 'plastic.type' must be specified")
    return(NA)}
  
  if(!is.null(plastic.type)){
    if(!plastic.type %in% c('fragment', 'fibre')){
      warning("'plastic.type' must be either 'fragment' or 'fibre'")
      return(NA)}
    # plastic.type (and available dimensions) determines shape
    switch(plastic.type,
           fibre = {
             shape_ <- 'cylinder'
             if(is.null(W)){
               W <- fibre.width.default
               if(!suppress.warnings) warning(
                 paste0("Fibre width, 'W', was not included as an argument so using the default W = ",
                        fibre.width.default, " μm. Note the units of this default and change argument 'fibre.width.default' if needed."))}
           },
           fragment = {
             if(is.null(W)) shape_ <- 'sphere' else{
               if(is.null(A)) shape_ <- 'ellipsoid' else shape_ <- 'ragged'}
           })
    if(!suppress.warnings && !is.null(shape) && shape != shape_) warning("Input 'shape' has been overwritten by choice of 'plastic.type'")
    shape <- shape_
  }
  
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


