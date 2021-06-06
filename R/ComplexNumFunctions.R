#' @name Magnitude
#' @title Get the magnitude of a complex number
#' @param x Input is a single complex number or a matrix of complex numbers.
Magnitude <- function (x)
{
  if (is.numeric(x) == FALSE & is.complex(x) == FALSE){
    stop('You need to input a damn number!', call. = FALSE)
  }
  if (is.complex(x) == FALSE){
    stop('This function requires a complex number input.', call. = FALSE)
  }
  magtemp <- abs(x)
}

#' @name AngleRads
#' @title Get the angle of a complex number (radians)
#' @param x Input is a single complex number or a matrix of complex numbers.
AngleRads <- function (x)
{
  if (is.numeric(x) == FALSE & is.complex(x) == FALSE){
    stop('You need to input a damn number!', call. = FALSE)
  }
  if (is.complex(x) == FALSE){
    stop('This function requires a complex number input.', call. = FALSE)
  }
  radtemp <- Arg(x)
}

#' @name AngleDegrees
#' @title Get the angle of a complex number (degrees)
#' @param x Input is a single complex number or a matrix of complex numbers.
AngleDegrees <- function (x)
{
  if (is.numeric(x) == FALSE & is.complex(x) == FALSE){
    stop('You need to input a damn number!', call. = FALSE)
  }
  if (is.complex(x) == FALSE){
    stop('This function requires a complex number input.', call. = FALSE)
  }
  degtemp <- Arg(x) * (180/pi)
}

#' @name DotPlotComplex
#' @title Plots complex number(s) in Cartesian space
#' @param x Input is a single complex number or a matrix of complex numbers.
DotPlotComplex <- function(x)
{
  if (is.numeric(x) == FALSE & is.complex(x) == FALSE){
    stop('You need to input a damn number!', call. = FALSE)
  }
  if (is.complex(x) == FALSE){
    stop('This function requires a complex number input.', call. = FALSE)
  }
  realTemp <- Re(x)
  imagTemp <- Im(x)
  CompPlot <- plot(realTemp,imagTemp)
}
