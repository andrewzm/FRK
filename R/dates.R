## Taken from package Hmisc
## Copyright notice from Hmisc:
## This program is free software; you can redistribute it and/or modify it
## under the terms of the GNU General Public License as published by the
## Free Software Foundation; either version 2, or (at your option) any
## later version.
##
## These functions are distributed in the hope that they will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## The text of the GNU General Public License, version 2, is available
## as http://www.gnu.org/copyleft or by writing to the Free Software
## Foundation, 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
##



yearDays <- function(time) {
  time <- as.POSIXlt(time)

  time$mon[] <- time$mday[] <- time$sec[] <- time$min <- time$hour <- 0
  time$year <- time$year + 1

  return(as.POSIXlt(as.POSIXct(time))$yday + 1)
}

monthDays <- function(time) {
  time <- as.POSIXlt(time)
  time$mday[] <- time$sec[] <- time$min <- time$hour <- 0
  time$mon <- time$mon + 1

  return(as.POSIXlt(as.POSIXct(time))$mday)
}

#' @export
round.POSIXt <- function(x, digits=c("secs", "mins", "hours", "days", "months", "years"))
  {
    ## this gets the default from the generic, as that has two args.
    if(is.numeric(digits) && digits == 0.0) digits <-"secs"
    units <- match.arg(digits)

    month.length <- monthDays(x)
    x <- as.POSIXlt(x)

    if(length(x$sec) > 0)
      switch(units,
             "secs"   = {x$sec <- x$sec + 0.5},
             "mins"   = {x$sec <- x$sec + 30},
             "hours"  = {x$sec <- 0; x$min <- x$min + 30},
             "days"   = {x$sec <- 0; x$min <- 0; x$hour <- x$hour + 12
                         isdst <- x$isdst <- -1},
             "months" = {x$sec <- 0; x$min <- 0; x$hour <- 0;
                         x$mday <- x$mday + trunc(monthDays(x)/2);
                         isdst <- x$isdst <- -1},
             "years"  = {x$sec <- 0; x$min <- 0; x$hour <- 0;
                         x$mday <- 0; x$mon <- x$mon + 6;
                         isdst <- x$isdst <- -1}
             )

    return(trunc(as.POSIXct(x), units=units))
  }

#' @export
trunc.POSIXt <- function(x, units=c("secs", "mins", "hours", "days", "months", "years"), ...) {
    units <- match.arg(units)

    x <- as.POSIXlt(x)

    isdst <- x$isdst
    if(length(x$sec) > 0)
      switch(units,
             "secs" = {x$sec <- trunc(x$sec)},
             "mins" = {x$sec <- 0},
             "hours"= {x$sec <- 0; x$min <- 0},
             "days" = {x$sec <- 0; x$min <- 0; x$hour <- 0; isdst <- x$isdst <- -1},
             "months" = {
               x$sec <- 0
               x$min <- 0
               x$hour <- 0
               x$mday <- 1
               isdst <- x$isdst <- -1
             },
             "years" = {
               x$sec <- 0
               x$min <- 0
               x$hour <- 0
               x$mday <- 1
               x$mon <- 0
               isdst <- x$isdst <- -1
             }
             )

    x <- as.POSIXlt(as.POSIXct(x))
    if(isdst == -1) {
      x$isdst <- -1
    }
    return(x)
  }

ceil <- function(x, units, ...) {
  UseMethod('ceil', x)
}

ceil.default <- function(x, units, ...) {
  ceiling(x)
}

ceil.POSIXt <- function(x, units=c("secs", "mins", "hours", "days", "months", "years"), ...) {
  units <- match.arg(units)

  x <- as.POSIXlt(x)

  isdst <- x$isdst
  if(length(x$sec) > 0 && x != trunc.POSIXt(x, units=units)) {
    switch(units,
           "secs" = {
             x$sec <- ceiling(x$sec)
           },
           "mins" = {
             x$sec <- 0
             x$min <- x$min + 1
           },
           "hours"= {x$sec <- 0; x$min <- 0; x$hour <- x$hour + 1},
           "days" = {
             x$sec <- 0
             x$min <- 0
             x$hour <- 0
             x$mday <- x$mday + 1
             isdst <- x$isdst <- -1
           },
           "months" = {
             x$sec <- 0
             x$min <- 0
             x$hour <- 0
             x$mday <- 1
             x$mon <- x$mon + 1
             isdst <- x$isdst <- -1
           },
           "years" = {
             x$sec <- 0
             x$min <- 0
             x$hour <- 0
             x$mday <- 1
             x$mon <- 0
             x$year <- x$year + 1
             isdst <- x$isdst <- -1
           }
           )

    x <- as.POSIXlt(as.POSIXct(x))
    if(isdst == -1) {
      x$isdst <- -1
    }
  }
  return(x)
}
