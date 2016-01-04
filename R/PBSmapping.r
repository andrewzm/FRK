#==============================================================================
# Copyright (C) 2003-2013  Fisheries and Oceans Canada
# Nanaimo, British Columbia
# This file is part of PBS Mapping.
#
# PBS Mapping is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# PBS Mapping is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with PBS Mapping; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
#==============================================================================
# PBS Mapping
#
# File:
#   PBSmapping.R
#
# Version:
#   See DESCRIPTION.
#
# Description:
#   Mapping software for R.
#
# Authors:
#   Specifications: Jon Schnute, Nicholas Boers, Rowan Haigh
#   Code:           Nicholas Boers (and others as documented)
#
# Change Log:
#    See PBSmapping/ChangeLog
#==============================================================================

PBSprint <- FALSE

#==============================================================================
.addAxis <- function(xlim, ylim, tckLab, tck, tckMinor, ...)
  # 'xlim'/'ylim': vectors of length 2
  # 'tckLab'/'tck'/'tckMinor': vectors of length 1 or 2
  #
  # Returns: NULL (invisible)
{
  # force 'tck*' to length 2
  tckLab   <- rep(tckLab,   length.out = 2);
  tck      <- rep(tck,      length.out = 2);
  tckMinor <- rep(tckMinor, length.out = 2);

  # reduce 'cex' and 'mex' to appropriate sizes for axis labels; these settings
  # should remain after the function termimates; .addLabels() will use the same
  # units for measuring lines

  # changing 'cex' causes plot(...) to change 'mai', so let's reset 'mai'
  # after setting 'cex'
  mai <- par()$mai;
  par(cex = par()$cex * 0.8);   # decrease font size
  par(mai = mai);
  par(mex = par()$cex);         # decrease spacing to match font size

  # 1 is the horizontal axis, 2 is the vertical axis
  lim <- list(xlim, ylim);
  rotate <- list(0, 90)
  for (i in 1:2) {
    if (((i == 1) && (par()$xaxt != "n")) ||
        ((i == 2) && (par()$yaxt != "n"))) {
      # create both major and minor ticks
      ticks <- pretty(lim[[i]]);
      ticksMinor <- pretty(c(0, diff(ticks)[1]));
      ticksMinor <- (sort(rep(ticks, length(ticksMinor))) +
                     rep(ticksMinor, length(ticks)));

      # filter them for anything on the extents
      ticks <- ticks[ticks > lim[[i]][1] & ticks < lim[[i]][2]];
      ticksMinor <-
        ticksMinor[ticksMinor > lim[[i]][1] & ticksMinor < lim[[i]][2]];

      if (!tckLab[i]) {
        tickLabels <- FALSE;
      } else {
        tickLabels <- as.character(ticks);
      }

      # plot the axis
      # major ticks
      axis(side = i, at = ticks, labels = tickLabels, tck = tck[i],
           srt = rotate[[i]]);
      # minor ticks
      axis(side = i, at = ticksMinor, labels = FALSE, tck = tckMinor[i]);
    }
  }

  invisible(NULL);
}

#==============================================================================
.addCorners <- function(polys, ptSummary)
{
  xlim <- range(polys$X)
  ylim <- range(polys$Y)

  corners <- list(tl=c(xlim[1], ylim[2]), tr=c(xlim[2], ylim[2]),
                  bl=c(xlim[1], ylim[1]), br=c(xlim[2], ylim[1]));
  # two "tests" for each corner:
  # each list element is the function to apply to X and then Y
  tests <- list(tl=c(min, max), tr=c(max, max),
                bl=c(min, min), br=c(max, min));

  for (c in names(corners)) {
    # do we need to add a corner point? check whether a point is already in the
    # corner
    if (nrow(polys[polys$X == (corners[[c]])[1] &
                   polys$Y == (corners[[c]])[2], ]) > 0)
      next

    # find candidate PIDs for the corner;
    # candidates touch a boundary on either side of the corner point
    polysA <- polys[polys$X == (corners[[c]])[1], ]
    if (nrow(polysA) > 0)
      # for those polys with an X on the boundary, which ones have a Y
      # that is the min/max? get their PIDs
      PIDsA <- polysA[polysA$Y == ((tests[[c]])[2])[[1]](polysA$Y), "PID"]

    polysB <- polys[polys$Y == (corners[[c]])[2], ]
    if (nrow(polysB) > 0)
      PIDsB <- polysB[polysB$X == ((tests[[c]])[1])[[1]](polysB$X), "PID"]

    if (nrow(polysA) == 0 || nrow(polysB) == 0) {
      warning(paste("Unable to close a corner (", c, ").", sep=""));
      next
    }

    # determine candidates
    cand <- intersect(PIDsA, PIDsB)

    # no candidates (this shouldn't happen...)
    if (length(cand) == 0)
      stop(paste("Unable to close a corner (", c, ") since no candidates exist.",
                 sep=""));

    # more than one candidate and we need to find the appropriate one
    if (length(cand) > 1) {
      # compute the distance from each candidate point to the corner
      pts <- data.frame (x=ptSummary[cand, "x"], y=ptSummary[cand, "y"])
      pts$cornerX <- (corners[[c]])[1];
      pts$cornerY <- (corners[[c]])[2];
      pts$dist <- sqrt((pts$x - pts$cornerX)^2 + (pts$y - pts$cornerY)^2)
      shortest <- which(pts$dist == min(pts$dist))
      if (length(shortest) != 1)
        stop(paste(
"Unable to determine the appropriate polygon to close corner ", c, ".", sep=""))
      cand <- cand[shortest]
    }

    # add corner point
    newPoly <- polys[polys$PID == cand, ]
    polys <- polys[polys$PID != cand, ]

    xydata <- data.frame(X=c(newPoly$X, corners[[c]][1]),
                         Y=c(newPoly$Y, corners[[c]][2]));

    newPoly <- calcConvexHull(xydata)
    newPoly$PID <- cand

    polys <- rbind(polys, newPoly);
  }

  polys <- polys[order(polys$PID), ]
  return (polys)
}

#==============================================================================

#.addBubblesLegend----------------------2012-03-01
#  Construct legend for the function addBubble.
#--------------------------------------------DC/NB
.addBubblesLegend <- function(radii.leg, usr.xdiff, usr.ydiff,
                              symbol.zero, symbol.fg, symbol.bg,
                              legend.pos, legend.breaks,
                              legend.type, legend.title, legend.cex, ...)
{
  # ratio of y to x: units-per-inch (Y) / units-per-inch (X)
  ratio.y.x = (usr.ydiff / par("pin")[2]) / (usr.xdiff / par("pin")[1])

  # calculate the height and width of the legend, which is essential
  # to calculating its position
  gap.x <- par("cxy")[1] * legend.cex / 2
  gap.y <- par("cxy")[2] * legend.cex / 2
  radii.leg.y <- radii.leg * ratio.y.x
  leg.tex.w <- strwidth(legend.breaks, units = "user") * legend.cex
  title.w = strwidth(legend.title)
  max.tex.w <- max(leg.tex.w)

  switch(legend.type,
         nested = {
           # height of the legend is the biggest bubble + title
           legend.height <- 2 * max(radii.leg.y) + 3 * gap.y
           legend.width <- 2 * max(radii.leg) + gap.x + max.tex.w },
         horiz = {
           # height of the legend is the biggest bubble + title + tag
           legend.height <- 2 * max(radii.leg.y) + 3 * gap.y
           legend.width <- 2 * sum(radii.leg) +
             (length(legend.breaks) - 1) * gap.x },
         vert = {
           legend.height <- 2 * sum(radii.leg.y) +
             (length(legend.breaks) - 1) * gap.y + 3 * gap.y
           legend.width <- 2 * max(radii.leg) + gap.x + max.tex.w }
         )

  # reflect an adjustment in X if the title is broader than the legend
  if (title.w > legend.width) {
    w.adj <- (title.w - legend.width) / 2
  } else {
    w.adj <- 0
  }

  # if we already have positions, keep them; otherwise, calculate
  # positions given the described corner
  if (class(legend.pos) == "numeric") {
    legend.loc <- legend.pos
  } else {
    corners <- c("bottomleft", "bottomright", "topleft", "topright")
    if (legend.pos %in% corners) {
      legend.loc <- switch(legend.pos,
        bottomleft =
          c(par("usr")[1] + 0.025 * usr.xdiff + w.adj,
            par("usr")[3] + 0.025 * usr.ydiff + legend.height),
        bottomright =
          c(par("usr")[2] - (0.025 * usr.xdiff + legend.width + w.adj),
            par("usr")[3] + 0.025 * usr.ydiff + legend.height),
        topleft =
          c(par("usr")[1] + 0.025 * usr.xdiff + w.adj,
            par("usr")[4] - 0.025 * usr.ydiff),
        topright =
          c(par("usr")[2] - (0.025 * usr.xdiff + legend.width + w.adj),
            par("usr")[4] - 0.025 * usr.ydiff));
    }
  }

  # the calling function should already have validated legend.type
  switch(legend.type,
         nested = {
           # legend.loc[1] specifies X for *center* of circle; we need
           # to shift it right by the largest radius
           legend.loc[1] <- legend.loc[1] + max(radii.leg)
           # the legend will be drawn from the bottom up; shift
           # legend.loc[2] so that it refers to the bottom
           legend.loc[2] <- legend.loc[2] - legend.height

           r <- rev(radii.leg)
           bb <- rev(legend.breaks)

           # position of the right edge of the text labels legend
           x.text.leg <- legend.loc[1] + max(r) + gap.x + max.tex.w

           # draw the circles, lines, and labels
           for (i in 1:length(r)) {
             symbols(legend.loc[1], legend.loc[2] + r[i] * ratio.y.x,
                     circles=r[i], inches=FALSE, add=TRUE, bg=symbol.bg[length(r)-i+1],
                     fg=symbol.fg)
             lines(c(legend.loc[1], legend.loc[1] + r[1] + gap.x),
                   rep(legend.loc[2] + 2 * r[i] * ratio.y.x, 2))
             text(x.text.leg, legend.loc[2] + 2 * r[i] * ratio.y.x, bb[i],
                  adj=c(1, .5), cex=legend.cex)
           }

           # add the title
           x.title.leg <- legend.loc[1] - max(radii.leg) + (legend.width / 2)
           text(x.title.leg, legend.loc[2]+legend.height, legend.title,
                adj=c(0.5,0.5), cex=legend.cex+0.2, col="black")

           # set positions for plotting of zero (later)
           zlab <- c(x.title.leg, legend.loc[2]+legend.height/4) },

         horiz = {
           # legend.loc[2] currently identifies the top of the legend
           legend.loc[2] <- legend.loc[2] + max(radii.leg.y) - legend.height

           # compute offsets for horizontal spacing
           offset <- vector()
           for (i in 1:length(radii.leg))
             offset[i] <- 2 * sum(radii.leg[1:i]) - radii.leg[i] +
               (i - 1) * gap.x

           # draw circles, labels
           symbols(legend.loc[1] + offset,
                   rep(legend.loc[2],length(radii.leg)),
                   circles = radii.leg, inches = FALSE, bg = symbol.bg,
                   fg = symbol.fg, add = TRUE)
           text(legend.loc[1] + offset, legend.loc[2] + radii.leg.y + gap.y,
                legend.breaks, adj = c(0.5, 0.5), cex = legend.cex)

           # add the title
           text(legend.loc[1] + legend.width / 2,
                legend.loc[2] + legend.height - max(radii.leg.y),
                legend.title, adj = c(0.5, 0.5), cex = legend.cex + 0.2,
                col = "black")

           # set positions for plotting of zero (later)
           zlab <- c(legend.loc[1], legend.loc[2] - legend.height / 8) },

         vert = {
           # part of the calculations have been made above to
           # calculate the height of the legend
           if (any(legend.pos == c("bottomleft","topleft")))
             legend.loc[1] <- legend.loc[1] + 0.05 * usr.xdiff

           # compute offsets for vertical spacing
           offset <- vector()
           for (i in 1:length(legend.breaks))
             offset[i] <- gap.y + 2 * sum(radii.leg.y[1:i]) - radii.leg.y[i] +
               i * gap.y

           # draw circles, labels
           symbols(rep(legend.loc[1], length(legend.breaks)),
                   legend.loc[2] - offset,
                   circles = radii.leg, bg = symbol.bg, fg = symbol.fg,
                   inches = FALSE, add = TRUE)
           x.text.leg <- legend.loc[1] + max(radii.leg) + gap.x + max.tex.w
           text(rep(x.text.leg, length(legend.breaks)), legend.loc[2] - offset,
                legend.breaks,
                cex = legend.cex, adj = c(1, 0.5), col="black")

           # add the title
           text(legend.loc[1] + legend.width / 2 - max(radii.leg),
                legend.loc[2], legend.title, adj = c(0.5, 0.5),
                cex = legend.cex + 0.2, col = "black")

           # set positions for plotting of zero (later)
           zlab <- c(legend.loc[1] + legend.width / 8, legend.loc[2]) }
         )

  # plot the zero if need be
  if (!is.logical(symbol.zero))
    legend(zlab[1], zlab[2], legend = "zero", pch = symbol.zero, xjust = 0,
           yjust = 1, bty = "n", cex = 0.8, x.intersp = 0.5)

  invisible()
}

#==============================================================================
.addFeature <- function(feature, data, polyProps, isEventData,
                        cex = NULL, col = NULL, font = NULL, pch = NULL, ...)
  # 'feature == "points"': add 'data' with points()
  # 'feature == "labels"': add 'data' with text()
  # 'data': if is PolyData, relaxed unique PID requirement (for addStipples())
  # 'isEventData': if TRUE, look for an EID column in 'data'
  # '...' contains arguments for either the 'text()' or 'points()' function
  #
  # Returns: PolyProps
{
  data <- .mat2df(data);
  if (isEventData) {
    type <- "e";
  } else {
    type <- "p";
  }

  # given 'feature', appropriate columns in polyProps
  if (feature == "points") {
    relevantProps <- c("cex", "col", "pch");
  } else {
    relevantProps <- c("cex", "col", "font");
  }

  # validate the polyProps argument
  polyProps <- .validatePolyProps(polyProps, parCols = relevantProps);
  if (is.character(polyProps))
    stop(paste("Invalid PolyData 'polyProps'.\n", polyProps, sep=""));

  # at this point, 'data' has:
  # (!isEventData): PID, (optional SID), X, Y, and label columns
  # (isEventData):  EID, X, Y, and label columns

  # defaults for 'polyProps'
  parValues <- list(cex = par("cex"), col = par("col"),
                    font = par("font"), pch = par("pch"));
  # don't replace existing values with defaults
  parValues <- parValues[setdiff(names(parValues), names(polyProps))];
  # keep only applicable values
  parValues <- parValues[intersect(names(parValues), relevantProps)];
  # override defaults with values passed as arguments
  if (!is.null(cex))  parValues[["cex"]]  <- cex;
  if (!is.null(col))  parValues[["col"]]  <- col;
  if (!is.null(font)) parValues[["font"]] <- font;
  if (!is.null(pch))  parValues[["pch"]]  <- pch;

  # make basic 'polyProps' if necessary
  if (is.null(polyProps)) {
    if (isEventData) {
      polyProps <- data.frame(EID = unique(data$EID));
    } else {
      polyProps <- data.frame(PID = unique(data$PID));
    }
  }
  # merge SIDs into PolyProps if necessary
  if (!isEventData
      && is.element("SID", names(data))
      && !is.element("SID", names(polyProps))) {
    # grab SIDs from 'polys' for all PIDs in 'polyProps'
    p <- data[is.element(data$PID, unique(polyProps$PID)), c("PID", "SID")];
    # by filtering things first (above), we speed up the paste() (below)
    # add those SIDs to 'polyProps'
    polyProps <- merge(polyProps,
                       p[!duplicated(paste(p$PID, p$SID)), c("PID", "SID")],
                       by="PID");
  }

  # flesh out 'polyProps'
  if (length(parValues) > 1)
    polyProps <- .addProps(type = type, polyProps = polyProps, parValues);
  polyPropsReturn <- polyProps;

  # reduce data to IDs found in 'polyProps'
  if (isEventData) {
    data <- data.frame(data[is.element(data$EID, unique(polyProps$EID)), ]);
  } else {
    if (is.element("SID", names(polyProps))) {
      # 'polyProps' may have an SID field and 'polys' may not
      if (!is.element("SID", names(data)))
        stop("Since 'polyProps' contains an SID column, 'data' must as well.\n");
      data <- data[is.element(paste(data$PID, data$SID),
                              unique(paste(polyProps$PID, polyProps$SID))), ];
    } else {
      data <- data[is.element(data$PID, unique(polyProps$PID)), ];
    }
  }

  # flatten relevant columns
  propColumns <- intersect(names(polyProps), relevantProps);
  # paste the columns together
  exprStr <- paste("paste(",
                   paste(paste('polyProps[, "', propColumns, '"]', sep=""),
                         collapse=", "),
                   ");", sep="");
  polyProps$props <- eval(parse(text=exprStr))

  # merge that 'props' column into data
  if (isEventData) {
    data <- merge(data, polyProps[, c("EID", "props")], by="EID");
  } else {
    # if SID exists in 'polyProps' it exists in 'data' (and vice-versa)
    if (is.element("SID", names(polyProps))) {
      data <- merge(data,
                    polyProps[, c("PID", "SID", "props")], by=c("PID", "SID"));
    } else {
      data <- merge(data, polyProps[, c("PID", "props")], by="PID");
    }
  }

  # split data on the 'props' column
  data <- split(data, data$props);
  # keep polyProps for getting at individual properties
  polyProps <- polyProps[!duplicated(polyProps$props), ];

  # plot each set of properties on its own
  for (c in names(data)) {
    d <- (data[[c]]);
    p <- (as.list(polyProps[polyProps$props == c, propColumns]));

    if (feature == "labels") {
      text(x = d$X, y = d$Y, labels = as.character(d$label),
           cex = p$cex, col = p$col, font = p$font, ...);
    } else {
      points (x = d$X, y = d$Y,
              cex = p$cex, col = p$col, pch = p$pch, ...);
    }
  }

  return (polyPropsReturn);
}

#==============================================================================
.addProps <- function(type, polyProps, ...)
  # cycles properties in '...' by PID
  # 'type = "e"': EventData
  # 'type = "p"': PolyData
  # 'polyProps': a data frame
  # '...': parameters in the form (col = 1, cex = 1, ...),
  #   each of which can be a list object
  #
  # If '...' contains parameters that already exist in 'polyProps', these
  # parameters will overwrite the existing columns in 'polyProps' -- avoid
  # this behaviour by not including properties in '...' that already exist
  # in 'polyProps'.
  #
  # Returns: polyProps
{
  # clean up param -- in case it contains one or more lists
  param <- list(...);
  newParam <- list();
  for (i in 1:length(param)) {
    if (is.list(param[[i]])) {
      newParam <- c(newParam, param[[i]]);
    } else {
      newParam <- c(newParam, param[i]);
    }
  }
  param <- newParam;

  if (type == "e") {
    polyProps$IDX <- polyProps$EID;
  } else if (type == "p") {
    # IDX is PID rather than paste(PID, SID) because we cycle by PID
    polyProps$IDX <- polyProps$PID;
  } else {
    stop (
"Unknown 'type'.  Must be either \"e\" or \"p\".");
  }

  uIDX <- unique(polyProps$IDX);

  for (c in names(param)) {
    # add it if it isn't NULL
    if (!is.null(param[[c]])) {
      # remove column if already exists
      if (is.element(c, names(polyProps))) {
        # use data.frame() to ensure it isn't reduced to a vector
        polyProps <- data.frame(polyProps[, !is.element(names(polyProps), c)]);
        # if reduced to 1 column, it loses the name; reset it
        if (ncol(polyProps) == 1) {
          if (type == "p") names(polyProps) <- "PID";
          if (type == "e") names(polyProps) <- "EID";
        }
      }

      # build new structure
      newColumn <- data.frame(IDX=uIDX);
      newColumn[, c] <- rep(param[[c]], length.out = length(uIDX));

      # merge new structure
      polyProps <- merge(polyProps, newColumn, by = "IDX");
    }
  }

  # remove IDX column; use data.frame() to ensure it isn't reduced to a vector
  polyProps <- data.frame(polyProps[, !is.element(names(polyProps), "IDX")]);
  # if reduced to 1 column, it loses the name; reset it
  if (ncol(polyProps) == 1) {
    if (type == "p") names(polyProps) <- "PID";
    if (type == "e") names(polyProps) <- "EID";
  }

  return (polyProps);
}

#==============================================================================
.calcDist <- function(polys)
  # Assumes 'polys' contains valid PolySet with 'projection' attribute
  # containing LL/UTM/1.
  #
  # Returns: distance vector (distances between each point)
{
  # calculate distance for UTM/1:1
  if (!is.null(attr(polys, "projection"))
      && !is.na(attr(polys, "projection"))
      && ((attr(polys, "projection") == "UTM")
          || (attr(polys, "projection") == 1))) {
    len <- nrow(polys)

    D <- c(sqrt((polys$X[1:(len-1)] - polys$X[2:len])^2
                 + (polys$Y[1:(len-1)] - polys$Y[2:len])^2),
            0);

  }
  # calculate distance for LL
  else if (!is.null(attr(polys, "projection"))
           && !is.na(attr(polys, "projection"))
           && (attr(polys, "projection") == "LL")) {
    # Equatorial radius 6,378.14 km
    # Polar radius 6,356.78 km
    # Mean radius 6,371.3 km
    # Sources:
    #   http://en.wikipedia.org/wiki/Earth
    #   http://en.wikipedia.org/wiki/Earth_radius
    R <- 6371.3;

    # degrees to radians
    polys[, c("X", "Y")] <- polys[, c("X", "Y")] * pi / 180.0
    len <- nrow(polys)

    # Source:
    #   http://www.census.gov/cgi-bin/geo/gisfaq?Q5.1
    # Algorithm originally in pseudocode.
    s0 <- 1:(len - 1);  # 's' for 's'hift
    s1 <- 2:len;        # (s0 + 1;)

    dlon <- polys$X[s1] - polys$X[s0];   # dlon = lon2 - lon1
    dlat <- polys$Y[s1] - polys$Y[s0];   # dlat = lat2 - lat1
    cosPolysY <- cos(polys$Y);
    # a = (sin(dlat/2))^2 + cos(lat1) * cos(lat2) * (sin(dlon/2))^2
    a <- (sin(dlat / 2))^2 + cosPolysY[s0] * cosPolysY[s1] *
      (sin(dlon / 2))^2;
    # c = 2 * arcsin(min(1, sqrt(a)))
    a <- sqrt(a);
    a[a > 1] <- 1;
    cc <- 2 * asin(a);
    # d = R * c
    D <- c(R * cc, 0);
  }
  # unknown projection
  else {
    stop(paste(
"Invalid projection attribute.  Supported projections include \"LL\",",
"\"UTM\", and 1.\n"));
  }

  return (D);
}

#==============================================================================
.calcOrientation <- function(polys)
  # Assumes 'polys' contains a valid PolySet.
  #
  # Returns:
  #   data frame (invisible) with 'orientation' column (-1 when
  #     counter-clockwise; 0 when N/A; +1 when clockwise)
  #   OR: NULL (invisible) (if no rows in output)
{
  inRows <- nrow(polys);
  # Memory requirement for output are lower than nrow(polys); in fact, they
  # are length(unique(paste(polys$PID, polys$SID))).
  # The extra time to compute the real memory requirements makes doing so
  # not worth it.
  outCapacity <- nrow(polys);

  # create the data structure that the C functions expect
  # Using $<col> notation seems faster than [, "col"] notation.
  if (!is.element("SID", names(polys))) {
    inID <- c(polys$PID, integer(length = inRows), polys$POS);
  } else {
    inID <- c(polys$PID, polys$SID, polys$POS);
  }
  inXY <- c(polys$X, polys$Y);

  # call the C function
  results <- .C("calcOrientation",
                inID = as.integer(inID),
                inXY = as.double(inXY),
                inVerts = as.integer(inRows),
                outID = integer(2 * outCapacity),
                outOrientation = double(outCapacity),
                outRows = as.integer(outCapacity),
                outStatus = integer(1),
                PACKAGE = "PBSmapping");
  # note: outRows is set to how much space is allocated -- the C function
  #       should consider this

  if (results$outStatus == 1) {
    stop(
"Insufficient physical memory for processing.\n");
  }
  if (results$outStatus == 2) {
    stop(paste(
"Insufficient memory allocated for output.  Please upgrade to the latest",
"version of the software, and if that does not fix this problem, please",
"file a bug report.\n",
               sep = "\n"));
  }

  # determine the number of rows in the result
  outRows <- as.vector(results$outRows);

  # extract the data from the C function results
  if (outRows > 0) {
    d <- data.frame(PID = results$outID[1:outRows],
                    SID = results$outID[(outCapacity+1):(outCapacity+outRows)],
                    orientation = results$outOrientation[1:outRows]);

    if (!is.element("SID", names(polys)))
      d$SID <- NULL;

    invisible(d);
  } else {
    invisible(NULL);
  }
}

#==============================================================================
.checkClipLimits <- function(limits)
{
  # Makes sure that the X & Y limits are within the bounds of the GSHHS databases
  if (limits[1] > limits[2])
    stop("xlim[1] is larger than xlim[2]")
  if (limits[3] > limits[4])
    stop("ylim[1] is larger than ylim[2]")
  if (limits[1] > 360 || limits[2] < -20)
    stop("xlim are outside of the range of c(-20,360)")
  if (limits[3] > 90 || limits[4] < -90)
    stop("ylim are outside of the range of c(-90,90)")
}

#==============================================================================
.checkProjection <- function(projectionPlot, projectionPoly)
{
  if (is.null(projectionPlot)) {
    projMapStr <- "NULL";
  } else {
    projMapStr <- as.character(projectionPlot);
  }
  if (is.null(projectionPoly)) {
    projPolyStr <- "NULL";
  } else {
    projPolyStr <- as.character(projectionPoly);
  }
  msg <- paste(
"The data's 'projection' attribute (", projPolyStr, ") differs from the\n",
"projection of the plot region (", projMapStr, ").\n", sep="");

  if (xor(is.null(projectionPlot), is.null(projectionPoly))) {
    warning(msg);
  } else if ((!is.null(projectionPlot) && !is.null(projectionPoly)) &&
             (xor(is.na(projectionPlot), is.na(projectionPoly)))) {
    warning(msg);
  } else if (!is.null(projectionPlot) && !is.null(projectionPoly) &&
             !is.na(projectionPlot) && !is.na(projectionPoly) &&
             (projectionPlot != projectionPoly)) {
    warning(msg);
  }
}

#==============================================================================
.checkRDeps <- function(caller = "unspecified", requires = NULL)
{
  if (is.null(version$language) || (version$language != "R")) {
    stop (paste ("
The function '", caller, "' requires several dependencies available only in R.\n",
"Please try again from within R.\n", sep=""));
  }

  err <- NULL;
  for (pkg in requires) {
    if (!require(pkg, character.only = TRUE)) {
      err <- append (err, pkg);
    }
  }
  if (!is.null (err)) {
    err <- paste (err, collapse="', '");
    stop (paste ("
The function '", caller, "' requires the package(s) '", err, "'.\n",
"Please install the package(s) and try again.\n", sep=""));
  }
}

#==============================================================================
.clip <- function(polys, xlim, ylim, isPolygons, keepExtra)
  # Does not validate 'polys'; called directly from addPolys since addPolys()
  # already validated data.
  #
  # Maintains extra attributes.
  #
  # Returns:
  #   data frame with 'oldPOS' column
  #   OR: NULL (if no rows in output)
{
  # if 'keepExtra', create PolyData from PolySet; after processing, merge these
  # data back into 'polys'
  if (keepExtra)
    pdata <- extractPolyData(polys);

  # save the attributes of the data frame
  attrNames <- setdiff(names(attributes(polys)),
                       c("names", "row.names", "class"));
  attrValues <- attributes(polys)[attrNames];

  inRows <- nrow(polys);
  outCapacity <- as.integer(2 * inRows);

  # create the data structures that the C function expects
  if (!is.element("SID", names(polys))) {
    inID <- c(polys$PID, integer(length = inRows), polys$POS);
  } else {
    inID <- c(polys$PID, polys$SID, polys$POS);
  }
  inXY <- c(polys$X, polys$Y);
  limits <- c(xlim, ylim);

  # call the C function
  results <- .C("clip",
                inID = as.integer(inID),
                inXY = as.double(inXY),
                inVerts = as.integer(inRows),
                polygons = as.integer(isPolygons),
                limits = as.double(limits),
                outID = integer(4 * outCapacity),
                outXY = double(2 * outCapacity),
                outRows = as.integer(outCapacity),
                outStatus = integer(1),
                PACKAGE = "FRK");
  # note: outRows is set to how much space is allocated -- the C function
  #       should take this into consideration

  if (results$outStatus == 1) {
    stop(
"Insufficient physical memory for processing.\n");
  }
  if (results$outStatus == 2) {
    stop(paste(
"Insufficient memory allocated for output.  Please upgrade to the latest",
"version of the software, and if that does not fix this problem, please",
"file a bug report.\n",
               sep = "\n"));
  }

  # determine the number of rows in the result
  outRows <- as.vector(results$outRows);

  # extract the data from the C function results
  if (outRows > 0) {
    # after creating a data frame several different ways, the following seemed
    # to be quickest
    d <- data.frame(PID = results$outID[1:outRows],
                    SID = results$outID[(outCapacity+1):(outCapacity+outRows)],
                    POS = results$outID[(2*outCapacity+1):(2*outCapacity+outRows)],
                    oldPOS = results$outID[(3*outCapacity+1):(3*outCapacity+outRows)],
                    X = results$outXY[1:outRows],
                    Y = results$outXY[(outCapacity+1):(outCapacity+outRows)]);

    # remove "SID" column if not in input
    if (!is.element("SID", names(polys)))
      d$SID <- NULL;

    # change oldPOSs of -1 to NA to meet specs.
    d$oldPOS[d$oldPOS == -1] <- NA;

    # if keepExtra, merge the PolyData back in...
    # 'merge' loses attributes!  Merge before restoring attributes
    if (keepExtra)
      d <- merge(x = d, y = pdata, all.x = TRUE,
                 by = intersect(c("PID", "SID"), names(d)));

    # restore the attributes
    attributes(d) <- c(attributes(d), attrValues);

    # a final test, to ensure the PolySet contains more than just
    # degenerate edges
    if (all(d$X == xlim[1]) || all(d$X == xlim[2]) ||
        all(d$Y == ylim[1]) || all(d$Y == ylim[2]))
      return(NULL);

    return(d);
  } else {
    return(NULL);
  }
}

#==============================================================================
.closestPoint <- function(pts, pt)
  # "pts": data frame with columns "X" and "Y"
  # "pt":  data frame with columns "X" and "Y", and ONE row
  #
  # Returns a vector of length "pts" where T indicates that the point is
  # closest to "pt".  Returns several Ts when several points equidistant.
{
  pts$Xorig <- pt$X
  pts$Yorig <- pt$Y

  pts$dist <- (pts$Xorig - pts$X) ^ 2 + (pts$Yorig - pts$Y) ^ 2

  return (pts$dist == min(pts$dist))
}

#==============================================================================
.createFastIDdig <- function(polysA, polysB = NULL, cols)
  # Determines the maximum number of digits in the second column of a data
  # frame.  If given two data frames ('polysA' and 'polysB', determines that
  # maximum between the two data frames.
  #
  # Assumptions:
  #   both 'cols' contain integers
  #
  # Arguments:
  #   'polysA': first data frame
  #   'polysB': second data frame, which may be missing one or both 'cols';
  #     if missing one or more 'cols', it's assumed they'll be derived from
  #     those in 'polysA'
  #   'cols': vector of length 2, listing columns to use
  #
  # Note:
  #   If it returns 0, then the first and second columns will not fit into
  #   a double, and they must be 'paste'd toegether.
  #
  # Returns:
  #   number of digits in SECOND column (max. of 'polysA' and 'polysB')
  #   OR NULL if the columns don't exist in the data frame(s)
{
  if ((length(cols) == 2)
      && all(is.element(cols, names(polysA)))) {
    digitsL <- floor(log10(max(polysA[[cols[1]]])) + 1);
    digitsR <- floor(log10(max(polysA[[cols[2]]])) + 1);

    # check 'polysB' as well
    if (!is.null(polysB)) {
      if (is.element(cols[1], names(polysB)))
        digitsL <- max(digitsL, floor(log10(max(polysB[[cols[1]]])) + 1));
      if (is.element(cols[2], names(polysB)))
        digitsR <- max(digitsR, floor(log10(max(polysB[[cols[2]]])) + 1));
    }

    # 'double' has 15 digits of precision (decimal), according to my
    # 'C Pocket Reference' by Prinz, P. and U. Kirch-Prinz
    if ((digitsL + digitsR) <= 15) {
      return (digitsR);
    } else {
      return (0);
    }
  } else {
    return (NULL);
  }
}

#==============================================================================
.createGridIDs <- function(d, addSID, byrow)
  # Create IDs for a grid according to the addSID and byrow arguments.
  #
  # Arguments:
  #  'd': PolySet (grid) created using addSID = T and byrow = T
  #  'addSID': if TRUE, include an SID column
  #  'byrow': if TRUE, increment PID along X
{
  if (addSID && !byrow) {
    # swap
    tmp <- d$PID;
    d$PID <- d$SID;
    d$SID <- tmp;
  } else if (!addSID && byrow) {
    d$PID <- (d$SID - 1) * (length(unique(d$X)) - 1) + d$PID;
    d$SID <- NULL;
  } else if (!addSID && !byrow) {
    d$PID <- (d$PID - 1) * (length(unique(d$Y)) - 1) + d$SID;
    d$SID <- NULL;
  }

  return (d);
}

#==============================================================================
.createIDs <- function(x, cols, fastIDdig = NULL)
  # Creates an IDs (or IDX) column from its input.
  #
  # Arguments:
  #   'x': data frame with one or more columns
  #   'cols': columns to use when creating the index; OK if individual
  #     columns are missing
  #   'fastIDdig': (optional) maximum number of digits in the second column, if
  #     only two integer columns, and the maximum number of digits between them
  #     is less than 15; often the output from '.createFastIDdig'; a user would
  #     want to pass in 'fastIDdig' when creating (some) matching indices for two
  #     different data frames
  #
  # Note:
  #   If 'fastIDdig' equals NULL, executes '.createFastIDdig' to create it (but
  #   only if two columns).
  #
  # Returns:
  #   index column if everything OK
  #   NULL if error
{
  # use 'is.element' instead of 'intersect' because 'intersect' gives no
  # guarantee of order (as far as I can see)
  presentCols <- cols[is.element(cols, names(x))]

  if (length(presentCols) == 1) {
    return (x[[presentCols]]);
  } else if (length(presentCols) == 2) {
    # if a fastIDdig wasn't passed in, try to create one
    if (is.null(fastIDdig)) {
      fastIDdig <- .createFastIDdig(polysA=x, polysB=NULL, cols=presentCols);
    }

    # if called the function above, 'fastIDdig' won't equal NULL (but might
    # equal 0)
    # if the user passed it in, it should only equal NULL if there is
    # only one ID column, in which case we shouldn't be in this branch
    # of the 'if'
    if (fastIDdig > 0) {
      return (as.double(x[[presentCols[1]]]
                        + (x[[presentCols[2]]] / 10^fastIDdig)));
    } else {
      return (paste(x[[presentCols[1]]], x[[presentCols[2]]], sep = "-"));
    }
  } else if (length(presentCols) > 2) {
    exprStr <- paste("paste(",
                     paste(paste("x$", presentCols, sep=""),
                           collapse=", "),
                     ");", sep="");
    return (eval(parse(text=exprStr)));
  }

  return (NULL);
}

#==============================================================================
.expandEdges <- function(polys, pts, xlim, ylim)
  # expects all arguments to be specified
{
  polyRange <- c(range(polys$X), range(polys$Y))
  ptsRange <- c(range(pts$X), range(pts$Y))

  # set "toFix" based on points outside the PolySet
  toFixPts <- c(ptsRange[1] < polyRange[1],
                ptsRange[2] > polyRange[2],
                ptsRange[3] < polyRange[3],
                ptsRange[4] > polyRange[4])
  # update "toFix" to extend PolySet to limits outside of its range
  toFixLim <- c(signif(xlim[1], 5) < signif(polyRange[1], 5),
                signif(xlim[2], 5) > signif(polyRange[2], 5),
                signif(ylim[1], 5) < signif(polyRange[3], 5),
                signif(ylim[2], 5) > signif(polyRange[4], 5))
  # the use of "signif" reduces the number of unnecessary edges in
  # "toFixLim" due to rouding error; digits = 6 seemed insufficient
  # to compensate for rounding error

  toFix <- toFixPts | toFixLim

  # if nothing to fix, return
  if (!any(toFix))
    return (polys)

  # side: 1 = left, 2 = right, 3 = bottom, 4 = top
  for (side in which(toFix)) {
    # a point that causes expansion
    if (side == 1) {
      PID <- which(pts$X < polyRange[side])
    } else if (side == 2) {
      PID <- which(pts$X > polyRange[side])
    } else if (side == 3) {
      PID <- which(pts$Y < polyRange[side])
    } else if (side == 4) {
      PID <- which(pts$Y > polyRange[side])
    } else {
      stop ("Internal error: unrecognized value of \"size\" in point check.")
    }

    # no points caused the expansion, so it must be a limit that causes it
    if (length(PID) == 0) {
      # must determine PID using another method
      if (side == 1 || side == 2) {
        PID <- which(.closestPoint(pts, data.frame(X=xlim[side], Y=mean(ylim))))
      } else if (side == 3 || side == 4) {
        PID <- which(.closestPoint(pts, data.frame(X=mean(xlim), Y=ylim[side - 2])))
      } else {
        stop ("Internal error: unrecognized value of \"size\" in limit check.")
      }
    }

    # sanity check
    if (length(PID) != 1)
      stop ("Internal error: unable to determine appropriate PID for expansion.")

    newPoly <- polys[polys$PID == PID, ]
    polys <- polys[polys$PID != PID, ]

    if (side == 1 || side == 2) {
      newXY <- data.frame(X=c(newPoly$X, rep(xlim[side], 2)),
                          Y=c(newPoly$Y, range(newPoly$Y)))
    } else if (side == 3 || side == 4) {
      newXY <- data.frame(X=c(newPoly$X, range(newPoly$X)),
                          Y=c(newPoly$Y, rep(ylim[side - 2], 2)))
    } else {
      stop ("Internal error: unrecognized value of \"size\" in data frame setup.")
    }

    newPoly <- calcConvexHull(newXY)
    newPoly$PID <- PID

    polys <- rbind(polys, newPoly)
  }

  polys <- polys[order(polys$PID), ]

  return (polys);
}

#==============================================================================
.fixGSHHSWorld <- function (world) {
  # store desired limits
  xlim <- range(world$X)
  ylim <- range(world$Y)
  ylim[1] <- -90

  # determine PID of Antarctica: we'll use it to extract the current Antarctica,
  # which we'll grow west/east, then clip, then merge with the other polygons
  event <- data.frame(EID = 1, X = 85, Y = -72)
  event <- findPolys(event, world)
  pid <- event$PID[1]

  # extract current Antarctica
  curAnt <- world[is.element(world$PID,pid), ]

  # ensure the points are ordered from left to right
  if (curAnt$X[1] > curAnt$X[nrow(curAnt)])
    curAnt <- curAnt[nrow(curAnt):1, ]

  # create a copy to the left and a copy to the right; add corners
  left <- curAnt[c(1,1:nrow(curAnt)), ]
  left$X <- left$X - 360
  left$Y[1] <- -90
  right <- curAnt[c(1:nrow(curAnt), nrow(curAnt)), ]
  right$X <- right$X + 360
  right$Y[nrow(right)] <- -90

  # merge to create a new, very wide Antarctica
  curAnt <- rbind(left, curAnt, right)
  curAnt$POS <- 1:nrow(curAnt)

  # clip
  curAnt <- clipPolys(curAnt, xlim, ylim)
  curAnt$oldPOS <- NULL

  # merge into the existing world
  world <- rbind(world[world$PID < pid, ], curAnt, world[world$PID > pid, ])
  row.names(world) <- 1:length(row.names(world))

  invisible(world)
}

#=============================================================================
.getBasename <- function (fn, ext)
  # If appropriate, remove the extension from 'fn' to obtain the shapefile
  # name without the extension.  When testing for file existance, use the
  # extension 'ext'.
{
  # if appending .shp does not give a valid file...
  if (!file.exists (paste(fn, ".", ext, sep=""))) {
    # attempt to remove extension and try again
    fn <- sub ("\\..{3}$", "", fn);

    if (!file.exists (paste(fn, ".", ext, sep=""))) {
      stop (paste ("Cannot find the file \"", fn, ".", ext, "\".", sep=""));
    }
  }

  return (fn);
}

#=============================================================================
.getGridPars <- function (polys, fullValidation = TRUE)
{
  res <- list();
  res$x <- sort(unique(polys$X));
  res$y <- sort(unique(polys$Y));

  lenx <- length(res$x);
  leny <- length(res$y);
  if ((lenx < 2) || (leny < 2)) {
    return (NULL);
  }

  # determine "addSID"
  if (is.element("SID", names(polys))) {
    res$addSID <- TRUE;
  } else {
    res$addSID <- FALSE;
  }

  # determine "byrow"
  if (lenx == 2 && leny == 2) {
    # special case: only one polygon; byrow does not matter
    res$byrow <- TRUE;
  } else if (lenx == 2) {
    # special case: the second polygon in polys is _above_ the
    # first
    if (polys[1, "PID"] + 1 == polys[5, "PID"]) {
      res$byrow <- FALSE;
    } else {
      res$byrow <- TRUE;
    }
  } else {
    # the second polygon in polys is to the _right_ of the
    # first
    if (polys[1, "PID"] + 1 == polys[5, "PID"]) {
      res$byrow <- TRUE;
    } else {
      res$byrow <- FALSE;
    }
  }

  res$projection <- attr(polys, "projection");
  res$zone <- attr(polys, "zone");

  if (fullValidation) {
    t <- makeGrid(x = res$x, y = res$y, byrow = res$byrow,
                  addSID = res$addSID, projection = res$projection,
                  zone = res$zone);
    if (is.character (all.equal (polys, t))) {
      return (NULL);
    }
  }

  return (res);
}

#==============================================================================
.initPlotRegion <- function(projection, xlim, ylim, plt)
  # Initialize the plot region, accounting for the aspect ratio.
  #
  # 'projection': projection
  #   "UTM"   => UTM        <-- as per specs.
  #   "LL"    => Geographic <-- as per specs.
  #   NA      => none       <-- for plotLines, plotPolys
  #   numeric => specifies aspect ratio
  # 'xlim': x-limits of the plot
  # 'ylim': y-limits of the plot
  # 'plt': requested plot size
  #
  # Returns: NULL (invisible)
  #
  # Notes:
  # In a new plot, mar/mai are 'in-sync' with pin/plt, meaning that changing
  # mar/mai causes pin/plt to change.  Manually setting pin/plt causes mar/mai
  # to fall 'out-of-sync', meaning that changing mar/mai will have no effect
  # on pin/plt.  In this function, do not cause pin/plt to fall out-of-sync
  # (i.e., never directly set pin/plt). However, if they are already
  # out-of-sync, honour the pin setting (which is consistent with what R/S-PLUS
  # does).
  #
  # Summary:
  # Read plt; whether in-sync or out-of-sync, this value reflects the
  # maximum plot region.  Calculate desired plot region, assume in-sync, and
  # set new region indirectly via mar/mai. If new region doesn't match the
  # expected region, it is out-of-sync, so set new region directly via pin/plt.
  #
  # If the user allows pin/plt to fall out-of-sync and sets plt to NULL, the
  # function will still produce the correct aspect ratio, but may not honour
  # the mai/mar par() parameters.
{
  # validate 'xlim' and 'ylim'
  if ((missing(xlim) || is.null(xlim)) || !is.vector(xlim) ||
      (length(xlim) != 2) ||
      (missing(ylim) || is.null(ylim)) || !is.vector(ylim) ||
      (length(ylim) != 2)) {
    stop(
"xlim and/or ylim is missing or invalid.\n");
  }
  if (diff(ylim) == 0 || diff(xlim) == 0) {
    stop(
"xlim/ylim must specify a region with area greater than 0.\n");
  }

  # since all high-level plot functions call this function, we can
  # advance to the next frame here, unless par(new=T) is specified.
  if(!par()$new)
    frame();

  # if 'plt' is NULL, the user doesn't want to use the default, so read its
  # current value from 'par' parameters
  if (is.null(plt))
    plt <- par("plt");

  # 'projection' cannot equal NULL
  if (is.null(projection))
    stop(
"'projection' argument must not equal NULL.\n");

  # if projection isn't NA, adjust for the aspect ratio; in this section
  # we must indirectly set plt/pin via mai/mar if possible; if that doesn't
  # work, we must set them directly
  if (!is.na(projection))
    {
      # use an additional scaling factor if LL
      xyRatio <- ifelse(projection == "LL", cos((mean(ylim) * pi) / 180), 1);

      # determine plot region aspect ratio; calculate indirectly using
      # fin and plt
      aspPlotRegion <-
        (par()$fin[1]*diff(plt[1:2])) / (par()$fin[2]*diff(plt[3:4]));
      if (is.infinite(aspPlotRegion))
        stop("Plot region must have an area greater than 0.\n");

      # determine desired aspect ratio; xyRatio != 1 only when proj == "LL"
      aspPolySet <- diff(xlim) / diff(ylim) * xyRatio;
      if (is.numeric(projection))
        aspPolySet <- aspPolySet / projection;
      if (is.infinite(aspPolySet) || (aspPolySet == 0))
        stop(paste(
"Either 'projection' is 0 or 'xlim'/'ylim' specify a region with an area",
"of 0.\n", sep="\n"));

      # adjust 'plt' to honour the aspect ratio
      if (aspPlotRegion < aspPolySet) {
        # shrink the Y direction
        pinX <- par()$fin[1]*diff(plt[1:2]);
        pinY <- pinX * aspPolySet^-1;
        toMove <- (diff(plt[3:4]) - pinY/par()$fin[2]) / 2;
        plt[3:4] <- plt[3:4] + c(toMove, -toMove);
      } else if (aspPlotRegion > aspPolySet) {
        # shrink the X direction
        pinY <- par()$fin[2]*diff(plt[3:4]);
        pinX <- pinY * aspPolySet;
        toMove <- (diff(plt[1:2]) - pinX/par()$fin[1]) / 2;
        plt[1:2] <- plt[1:2] + c(toMove, -toMove);
      }

      # set the plot region indirectly via mai/mar
      # mai: (bottom, left, top, right)
      # fin: (width, height)
      # plt: (left, right, bottom, top)
      # If out-of-sync, may issue a warning re: plot specified in inches
      # too large.
      par(mai=c(par()$fin[2]*plt[3], par()$fin[1]*plt[1],
            par()$fin[2]*(1-plt[4]), par()$fin[1]*(1-plt[2])));

      # verify that it set OK; if not, set it directly via plt
      parPlt <- signif(as.double(par()$plt), digits=5);
      locPlt <- signif(as.double(plt), digits=5);
      if (any(parPlt != locPlt))
        par(plt = plt);
    }
  # if projection is NA, continue at this point with no changes to
  # the plot region

  # must set 'usr' again because changing the other par() parameters changes it
  par(usr = c(xlim, ylim));

  invisible(NULL);
}

#==============================================================================
.insertNAs <- function(polys, idx) {
  # extract the polys of interest
  sel.polys <- polys[is.element(names(polys), idx)]
  # NOTE: the above code is _much_ faster in S-PLUS than an alternative...
  #  sel.polys <- polys[as.character(idx)]
  # ... which I tried first
  if (is.null(unlist(sel.polys)))
    return (NA)

  # lapply's are slower in S-PLUS than R; the next line is a bottle-neck
  # in S-PLUS (while having minimal impact on R runtime)
  lenPolys <- lapply(sel.polys, "length")
  nPolys <- length(lenPolys)

  # create a vector T, F, ..., T
  TFT <- rep(c(TRUE, FALSE), length.out=(nPolys * 2) - 1)

  # create a vector with the number of vertices/polygon, with each separated
  # by a 1
  reps <- rep(1, len=(nPolys * 2) - 1)
  reps[TFT] <- unlist(lenPolys)

  # create a vector that can select all of the NAs that seperate polygons
  # in the "new.polys" vector
  NAs <- rep(!TFT, times=reps)

  new.polys <- vector(length=(length(unlist(sel.polys)) + nPolys - 1))
  new.polys[NAs] <- NA
  new.polys[!NAs] <- unlist(sel.polys)

  return (new.polys)
}

#==============================================================================
.mat2df <- function(data) {
  # Convert matrices to data.frames, preserving certain attributes.
  # Attributes to ignore: dim, dimnames, class
  pbsClass <- intersect(attributes(data)$class, c("EventData", "LocationSet",
              "PolyData", "PolySet"));
  attrNames <- names(attributes(data));
  if(all(is.element(attrNames, c("dim", "dimnames", "class")) == TRUE)) {
    addValues <- NULL;
  } else {
    addNames <- setdiff(attrNames, c("dim", "dimnames", "class"));
    addValues <- attributes(data)[addNames];
  }
  data <- data.frame(unclass(data));
  if(!is.null(addValues)) {
    attributes(data) <- c(attributes(data), addValues);
  }
  if(length(pbsClass) > 0) {
    attr(data, "class") <- c(pbsClass, "data.frame");
  }
  return(data);
}

#==============================================================================
.rollupPolys <- function(polys, rollupMode, exteriorCCW, closedPolys,
                         addRetrace)
  # rollupMode: method for rolling up the PolySet; essentially controls when
  #   to introduce a new PID (or PID/SID)
  #    1 = roll-up to the PID level (only PIDs in the result)
  #    2 = roll-up to the outer contour level (only outer contours in the
  #        result)
  #    3 = do not roll-up
  #
  # exteriorCCW: modify vertices orientation (CW/CCW)?
  #   -1 = don't modify
  #    0 = exterior should be CW
  #   +1 = exterior should be CCW
  #
  # closedPolys: whether the last and first vertices should be the same
  #   -1 = don't modify
  #    0 = ensure polygons do not close
  #   +1 = close the polygons
  #
  # addRetrace: determines whether it adds retrace lines to the first vertex
  # of the parent after outputting a child
  #    0 = don't add
  #    1 = add
  #
  # Note: does not validate the PolySet
  #
  # Returns: PolySet (invisible) or NULL (invisible)
{
  # save the attributes for the data frame (.validatePolySet returns a data
  # frame)
  attrNames <- setdiff(names(attributes(polys)),
                       c("names", "row.names", "class"));
  attrValues <- attributes(polys)[attrNames];

  inRows <- nrow(polys);
  outCapacity <- 2 * inRows;

  # create the data structure that the C function expects
  if (!is.element("SID", names(polys))) {
    inID <- c(polys$PID, integer(length = inRows));
  } else {
    inID <- c(polys$PID, polys$SID);
  }
  inXY <- c(polys$X, polys$Y);

  results <- .C("rollupPolys",
                inID = as.integer(inID),
                inPOS = as.double(polys$POS),
                inXY = as.double(inXY),
                inVerts = as.integer(inRows),
                outID = integer(3 * outCapacity),
                outXY = double(2 * outCapacity),
                outRows = as.integer(outCapacity),
                rollupMode = as.integer(rollupMode),
                exteriorCCW = as.integer(exteriorCCW),
                closedPolys = as.integer(closedPolys),
                addRetrace = as.integer(addRetrace),
                outStatus = integer(1),
                PACKAGE = "PBSmapping");
  # note: outRows is set to how much space is allocated -- the C function
  #       should take this into consideration

  if (results$outStatus == 1) {
    stop(
"Insufficient physical memory for processing.\n");
  }
  if (results$outStatus == 2) {
    stop(paste(
"Insufficient memory allocated for output.  Please upgrade to the latest",
"version of the software, and if that does not fix this problem, please",
"file a bug report.\n",
               sep = "\n"));
  }
  if (results$outStatus == 3) {
    stop(paste(
"Unable to rollup the polygons, as one or more children did not have a",
"parent.\n"));
  }

  # determine the number of rows in the result
  outRows <- as.vector(results$outRows);

  # extract the data from the C function results
  if (outRows > 0) {
    d <- data.frame(PID = results$outID[1:outRows],
                    SID = results$outID[(outCapacity+1):(outCapacity+outRows)],
                    POS = results$outID[(2*outCapacity+1):(2*outCapacity+outRows)],
                    X = results$outXY[1:outRows],
                    Y = results$outXY[(outCapacity+1):(outCapacity+outRows)]);

    if (!is.element("SID", names(polys)) || rollupMode == 1)
      d$SID <- NULL;

    # restore the attributes
    attributes(d) <- c(attributes(d), attrValues);

    invisible(d);
  } else {
    invisible(NULL);
  }
}

#==============================================================================
# Default for 'projection' must be logical type.
# For 'plotPoints()', 'cex'and 'pch' are part of '...' -- and they may equal
# NULL; when they equal NULL, don't try adding them to par()!
.plotMaps <- function(polys, xlim, ylim, projection, plt, polyProps,
                      border, lty, col, colHoles, density, angle, bg,
                      axes, tckLab, tck, tckMinor, isType, ...)
{
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # PART 1: deal with the '...' business; backup original 'par' values, and set
  #   values specified in '...' in 'par'
  # DO NOT backup 'new' (and some others)
  legalNames <- c("adj", "ann", "ask", "bg", "bty", "cex", "cex.axis",
                  "cex.lab", "cex.main", "cex.sub", "col", "col.axis",
                  "col.lab", "col.main", "col.sub", "crt", "csi", "err",
                  "exp", "fg", "font", "font.axis", "font.lab",
                  "font.main", "font.sub", "lab", "las", "lty",
                  "lwd", "mgp", "mkh", "pch", "smo", "srt", "tck",
                  "tcl", "tmag", "type", "xaxp", "xaxs", "xaxt", "xpd",
                  "yaxp", "yaxs", "yaxt");
  if (!is.null(version$language) && (version$language == "R")) {
    legalNames <- setdiff(legalNames, "csi");   # read-only in R
  }
  legalNames <- intersect(legalNames, names(par()));

  backupPar <- par(legalNames);
  on.exit(par(backupPar));

  dots <- list(...);

  extraArgs <- setdiff(names(dots), legalNames);
  extraArgs <- setdiff(extraArgs, c("main", "sub", "type", "xlab", "ylab"));
  if (length(extraArgs) > 0) {
    warning(paste(
"Ignored unrecognized argument '", paste(extraArgs, collapse="', '"), "'.\n",
                  sep = ""));
  }

  # special case: since they aren't explicit parameters, they end up being
  # part of '...'; but they may equal NULL (invalid values); remove them, since
  # they are explicitly passed into addPoints(), where they become part
  # of 'polyProps'
  if (isType == "points")
    legalNames <- setdiff(legalNames, c("cex", "pch"));

  par(dots[intersect(names(dots), legalNames)]);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # PART 2: set 'projection' and 'labelprojection' to the appropriate values
  #  Since we validate here, it's unnecessary to validate when this function
  #    calls 'addLines'/'addPolys'

  # although 'projection' cannot equal NA from caller, we may set it
  # to NA within this function (to mean unprojected)
  if (is.null(projection) || is.na(projection))
    projection <- FALSE;

  if (is.logical(projection)) {
    if (projection) {
      if (!is.null(attr(polys, "projection"))) {
        projection <- attr(polys, "projection");
      } else {
        projection <- 1;
        warning(
"'projection' set to 1:1 since unspecified 'projection' argument/attribute.\n");
      }
    } else {
      projection <- NA;
    }
  } else {
    # user specified 'projection': should be numeric, "LL", "UTM", and should
    # not conflict with 'polys'
    if (is.numeric(projection) || is.element(projection, c("LL", "UTM"))) {
      if (!is.null(attr(polys, "projection")) &&
          !is.element(projection, attr(polys, "projection"))) {
        projection <- attr(polys, "projection");
        warning(
"'projection' argument overwritten with PolySet's 'projection' attribute.\n");
      }
    } else {
      stop(paste(
"Either omit 'projection' argument or set it to a numeric value, \"LL\", or",
"\"UTM\".\n",
                 sep="\n"));
    }
  }
  # prevent duplicate warnings on call to 'addLines'/'addPolys'
  if (!is.null(polys))
    attr(polys, "projection") <- projection;
  labelProjection <- projection;

  # 'projection' type no longer 'logical'; now
  #   NA: no projection
  #   numeric: specified projection
  #   "LL": longitude-latitude
  #   "UTM": UTM

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # PART 3: get on with the routine

  # special case of NULL PolySet
  if (is.null(polys)) {
    if (is.null(xlim) || is.null(ylim) || is.null(projection)) {
      stop(
"To plot a NULL PolySet, pass 'xlim', 'ylim', and 'projection' arguments.\n");
    }
  }
  else {
    # as far as this function is concerned, 'polys' must have an X/Y, and
    # that's it; do a basic check
    polys <- .validateXYData(polys);
    if (is.character(polys)) stop(paste("Invalid PolySet.\n", polys, sep=""));
  }

  # detect limits if necessary
  if (is.null(xlim)) xlim <- range(polys$X);
  if (is.null(ylim)) ylim <- range(polys$Y);

  # check 'dots' for 'type'
  if (is.element("type", names(dots))) {
    if (dots$type == "n") {
      # set 'polys' to NULL since we've decided that type == "n" is equivalent
      # to passing a NULL PolySet
      polys <- NULL;
    } else {
      stop(
"Either omit 'type' argument or set it to \"n\".\n");
    }
  }

  # set 'col' parameter appropriately for adding labels/etc.
  if (!is.null(polys)) {
    par(col = 1);
  } else {
    if (length(col) > 1) {
      stop(paste(
"Either omit 'col' argument or set it to a single-element vector when 'polys'",
"equals NULL or 'type = \"n\"'.\n",
                 sep = "\n"));
    } else if (!is.null(col)) {
      par(col = col);
    }
  }

  # save settings in 'options'
  options(map.xlim = xlim);
  options(map.ylim = ylim);
  options(map.projection = projection);

  # create plot region
  .initPlotRegion(projection=projection, xlim=xlim, ylim=ylim, plt=plt);

  # plot background colour
  if (!is.null(bg))
    polygon(x = c(xlim[1], xlim[2], xlim[2], xlim[1]),
            y = c(ylim[1], ylim[1], ylim[2], ylim[2]),
            col = bg, border = 0);

  # plot PolySet 'polys'
  if (!is.null(polys)) {
    if (isType == "polygons") {
      # add polygons
      ret <- addPolys(polys, xlim = xlim, ylim = ylim, polyProps = polyProps,
                      border = border, lty = lty, col = col, colHoles = colHoles,
                      density = density, angle = angle);
    } else if (isType == "lines")  {
      # add lines
      ret <- addLines(polys, xlim = xlim, ylim = ylim, polyProps = polyProps,
                      lty = lty, col = col);
    } else if (isType == "points") {
      cex <- list(...)$cex;     # these were passed in '...'
      pch <- list(...)$pch;
      ret <- addPoints(polys, xlim = xlim, ylim = ylim, polyProps = polyProps,
                       cex = cex, col = col, pch = pch);
    } else {
      stop(
"Unrecognized 'isType'.\n");
    }
  } else {
    ret <- NULL;
  }

  if (axes) {
    .addAxis(xlim = xlim, ylim = ylim, tckLab = tckLab, tck = tck,
             tckMinor = tckMinor, ...);
  }
  else {
    options(map.xline = 1);
    options(map.yline = 1);
  }

  # labels must go after axis
  .addLabels(projection = labelProjection, ...);

  if (axes) {
    # since R won't plot outside the figure region, we don't want the
    # box to have a thicker line width
    box();
  }

  invisible(ret);
}

#==============================================================================
# .preparePolyProps:
#   Performs at least the following tasks:
#   1) creates 'polyProps' if it equals NULL
#   2) adds an SID column to 'polyProps' if one exists in 'polys'
#   3) removes from 'polyProps' any PIDs that do not exist in 'polys'
.preparePolyProps <- function(polysPID, polysSID, polyProps)
{
  # make 'polyProps' if necessary
  if (is.null(polyProps)) {
    polyProps <- data.frame(PID = unique(polysPID));
  } else {
    polyProps <- .validatePolyData(polyProps);
    if (is.character(polyProps))
      stop(paste("Invalid PolyData 'polyProps'.\n", polyProps, sep=""));
  }

  # if SIDs in 'polys' but not 'polyProps', add SIDs from 'polys' into
  # 'polyProps'
  if (!is.null(polysSID) && !is.element("SID", names(polyProps))) {
    # identify unique PIDs/SIDs within "polys"
    p <- data.frame(PID=polysPID, SID=polysSID);
    p <- p[!duplicated(.createIDs(p, cols = c("PID", "SID"))), ]

    # only keep PIDs that appear in 'polyProps'
    p <- p[is.element(p$PID, unique(polyProps$PID)), ];

    # add the SIDs to 'polyProps'
    polyProps <- merge(polyProps, p, by="PID");
  }

  return (polyProps);
}

#==============================================================================
.validateData <- function(data, className,
                          requiredCols = NULL, requiredAttr = NULL,
                          noFactorCols = NULL, noNACols = NULL, keyCols = NULL,
                          numericCols = NULL)
# An element of noNACols and keyCols will only be used if it exists in the
# data.  To ensure it exists in the data, make it a requiredCol.
{
  # convert matrix to data frame
  if (is.matrix(data)) {
    data <- .mat2df(data);
  }

  if (is.data.frame(data) && (nrow(data) > 0)) {
    # validate optional class name
    if (!is.null(className) && (class(data)[1] != "data.frame")) {
      if (class(data)[1] != className) {
        return(paste("Unexpected class (", class(data)[1], ").\n", sep=""));
      }
    }

    # ensure all the required columns exist in the PolySet
    if (!is.null(requiredCols) &&
        !all(is.element(requiredCols, names(data)))) {
      return(paste("One or more of the required columns is missing.\n",
      "Required columns: ", paste(requiredCols, collapse = ", "), ".\n", sep=""));
    }

    # ensure all the required attributes exists in the PolySet
    if (!is.null(requiredAttr) &&
        !all(is.element(requiredAttr, names(attributes(data))))) {
      return(paste("One or more of the required attributes is missing.\n",
      "Required attributes: ", paste(requiredAttr, collapse = ", "), ".\n", sep=""));
    }

    # check for NAs
    presentCols <- intersect(noNACols, names(data));
    if (length(presentCols) > 0) {
      # build an expression
      exprStr <- paste(paste("any(is.na(data$", presentCols, "))", sep=""), collapse=" || ");
      if (eval(parse(text=exprStr))) {
        return(paste("One or more columns (where NAs are not allowed) contains NAs.\n",
        "Columns that cannot contain NAs: ", paste(presentCols, collapse = ", "),
        ".\n", sep=""));
      }
    }

    # check for factors
    presentCols <- intersect(noFactorCols, names(data));
    if (length(presentCols) > 0) {
      # build an expression
      exprStr <- paste(paste("is.factor(data$", presentCols, ")", sep=""), collapse=" || ");
      if (eval(parse(text=exprStr))) {
        return(paste("One or more columns contains factors where they are not allowed.\n",
        "Columns that cannot contain factors: ", paste(presentCols, collapse = ", "),
        ".\n", sep=""));
      }
    }

    # check for uniqueness of the keys
    presentCols <- intersect(keyCols, names(data));
    if (length(presentCols) > 0) {
      if (length(presentCols) == 1) {
        keys <- data[[presentCols]];
      } else if ((length(presentCols) == 2)
                 && ((length(intersect(presentCols, c("PID","SID","POS","EID"))) == 2)
                     || (all(is.integer(data[[presentCols[1]]]))
                         && all(is.integer(data[[presentCols[2]]]))))) {
        # additional tests above to "ensure" the two columns contain integers
        keys <- .createIDs(data, cols=presentCols);
      } else {
        # paste the columns together
        exprStr <- paste("paste(",
                   paste(paste("data$", presentCols, sep=""), collapse=", "),");", sep="");
        keys <- eval(parse(text=exprStr));
      }

      # at this point, 'keys' is a vector
      if (any(duplicated(keys))) {
        return(paste("The 'key' for each record is not unique.\n",
        "Columns in key: ", paste(presentCols, collapse = ", "), ".\n", sep=""));
      }
    }

    # check for numeric columns
    presentCols <- intersect(numericCols, names(data));
    if (length(presentCols) > 0) {
      exprStr <- paste(paste("any(!is.numeric(data$", presentCols, "))", sep=""), collapse=" || ");
      if (eval(parse(text=exprStr))) {
        return(paste("One or more columns requires numeric values, but contains non-numerics.\n",
        "Columns that must contain numerics: ", paste(presentCols, collapse = ", "),
        ".\n", sep=""));
      }
    }

    # check for increasing/descreasing POS
    if (!is.null(className) && className == "PolySet") {
      idx <- .createIDs(data, cols=c("PID", "SID"));
      idxFirst <- which(!duplicated(idx));
      idxLast <- c((idxFirst-1)[-1], length(idx));
      # identify the holes
      holes <- (data$POS[idxFirst] > data$POS[idxLast])
      # outer/inner contour indices
      idxOuter <- rep(!holes, times=((idxLast-idxFirst)+1))
      idxInner <- !idxOuter;

      # POS[i] < POS[i+1]?
      lt <- c(data$POS[1:(nrow(data)-1)] < data$POS[2:(nrow(data))], FALSE);

      # check outer contours; change last vertex of each polygon to
      # what we expect for valid outer contours
      lt[idxLast] <- TRUE;
      # check for any that aren't in order
      j <- any(!lt[idxOuter])
      if (j) {
        j <- !lt;
        j[idxInner] <- FALSE;
        # add 1 because it's actually the next row that break it
        j <- which(j) + 1;

        return(paste("POS column must contain increasing values for outer contours.\n",
        "Offending rows: ", paste(j, collapse=", "), ".\n",  sep=""));
      }

      # check inner contours; change last vertex of each polygon to what
      # we expect for valid inner contours
      lt[idxLast] <- FALSE;
      # check for any that aren't in order
      j <- any(lt[idxInner])
      if (j) {
        j <- lt;
        j[idxOuter] <- FALSE;
        # do not add 1
        j <- which(j);

        return(paste("POS column must contain decreasing values for inner contours.\n",
        "Offending rows: ", paste(j, collapse=", "), ".\n",  sep=""));
      }
    }
  } else {
    return(paste("The object must be either a matrix or a data frame.\n"));
  }
  return(data);
}

#==============================================================================
.validateEventData <- function(EventData)
  # Perform some simple tests on the object to see if it can possibly be
  # an EventData object.
  # If the object is invalid, returns the error message.
{
  return(.validateData(EventData,
                       className = "EventData",
                       requiredCols = c("EID", "X", "Y"),
                       requiredAttr = NULL,
                       noFactorCols = c("EID", "X", "Y"),
                       noNACols = c("EID", "X", "Y"),
                       keyCols = c("EID"),
                       numericCols = c("EID", "X", "Y")));
}

#==============================================================================
.validateLocationSet <- function(LocationSet)
  # Perform some simple tests on the object to see if it can possibly be
  # an EventData object.
  # If the object is invalid, returns the error message.
{
  return(.validateData(LocationSet,
                       className = "LocationSet",
                       requiredCols = c("EID", "PID", "Bdry"),
                       requiredAttr = NULL,
                       noFactorCols = c("EID", "PID", "SID", "Bdry"),
                       noNACols = c("EID", "PID", "SID", "Bdry"),
                       keyCols = c("EID", "PID", "SID"),
                       numericCols = c("EID", "PID", "SID")));
}

#==============================================================================
.validatePolyData <- function(PolyData)
  # Perform some simple tests on the object to see if it can possibly be
  # a PolyData object.
  # If the object is invalid, returns the error message.
{
  return(.validateData(PolyData,
                       className = "PolyData",
                       requiredCols = c("PID"),
                       requiredAttr = NULL,
                       noFactorCols = c("PID", "SID"),
                       noNACols = c("PID", "SID"),
                       keyCols = c("PID", "SID"),
                       numericCols = c("PID", "SID")));
}

#==============================================================================
# 'polyProps' can be event properties, rather than polygon properties
.validatePolyProps <- function(polyProps, parCols = NULL)
{
  if (is.null(polyProps))
    return (NULL);

  return (.validateData(polyProps,
                        className = "PolyData",
                        requiredCols = NULL,            # PID/SID or EID
                        requiredAttr = NULL,
                        noFactorCols = parCols,
                        noNACols = NULL,
                        keyCols = NULL,
                        numericCols = NULL));
}

#==============================================================================
.validatePolySet <- function(polys)
  # Perform some simple tests on the object to see if it can possibly be
  # a PolySet object.
  # If the object is invalid, returns the error message.
{
  return(.validateData(polys,
                       className = "PolySet",
                       requiredCols = c("PID", "POS", "X", "Y"),
                       requiredAttr = NULL,
                       noFactorCols = c("PID", "SID", "POS", "X", "Y"),
                       noNACols = c("PID", "SID", "POS", "X", "Y"),
                       keyCols = c("PID", "SID", "POS"),
                       numericCols = c("PID", "SID", "POS", "X", "Y")));
}

#==============================================================================
.validateXYData <- function(xyData)
  # Perform some simple tests on the object to see if it can possibly be
  # an xyData object.
  # If the object is invalid, returns the error message.
{
  return(.validateData(xyData,
                       className = NULL,
                       requiredCols = c("X", "Y"),
                       requiredAttr = NULL,
                       noFactorCols = c("X", "Y"),
                       noNACols = c("X", "Y"),
                       keyCols = NULL,
                       numericCols = c("X", "Y")));
}

#==============================================================================
#addBubbles-----------------------------2010-01-07
#  addBubbles takes EventData and optional arguments to draw
#    bubbles of radius proportional to the z variable.
#
#  Modified (and for the legend, strongly inspired) from:
#    S. Tanimura, C. Kuroiwa, and T. Mizota. Proportional symbol
#    mapping in R. Journal of Statistical Software, 15(5):1-7,
#    Jan. 2006. [http://www.jstatsoft.org]
#
#  Modifications by Denis Chabot allow it to work with PBSmapping,
#  add one type of bubble (z proportional to volume of bubble)
#  and make it possible to draw several maps with bubbles that all
#  have the same radii (instead of each bubble plot having a radii
#  that depends on the max z value for that plot). Can add a
#  legend in one of 4 corners or at a specific x-y positiion.
#--------------------------------------------DC/RH
addBubbles <- function(events, type = c("perceptual", "surface", "volume"),
                       z.max = NULL, min.size = 0, max.size = 0.8, symbol.zero = "+",
                       symbol.fg = rgb(0,0,0,0.60), symbol.bg = rgb(0,0,0,0.30),
                       legend.pos = "bottomleft", legend.breaks = NULL,
                       show.actual = FALSE,
                       legend.type = c("nested", "horiz", "vert"),
                       legend.title = "Abundance", legend.cex = .8, ...) {
  # validate events
  events <- .validateEventData(events)
  if (is.character(events))
    stop(paste("Invalid EventData 'events'.\n", events, sep=""));
  if (!is.element("Z", names(events)))
    stop ("EventData is missing required column 'Z'.\n");

  # check arguments before we get too far
  type <- match.arg(type)
  if (!is.null(legend.pos))
    legend.type <- match.arg(legend.type)

  # set z.max if necessary
  if (is.null(z.max) || is.na(z.max))
    z.max <- max(events$Z, na.rm=TRUE)

  # adjust legend breaks if necessary
  if (is.null(legend.breaks) || is.na(legend.breaks))
    legend.breaks <- pretty(range(events$Z), 3)[-1]
  else if (is.vector(legend.breaks) && length(legend.breaks) == 1)
    legend.breaks <- pretty(range(events$Z), legend.breaks)[-1]

  if (show.actual)
    legend.breaks <- signif(legend.breaks / max(legend.breaks)
                            * max(events$Z, na.rm=TRUE), 3)

  # determine x/y range of plotting region
  usr.xdiff <- par("usr")[2] - par("usr")[1]
  usr.ydiff <- par("usr")[4] - par("usr")[3]

  # for sizing in inches, it's important to use the X rather than Y axis
  #
  # max.size is diameter (inches); /2 for radius; /par()$pin[1] (inches)
  # for fraction of width (inches); *usr.xdiff to convert to width to
  # user coordinates
  #
  # min.size is diameter (inches)
  stand.rad <- (max.size / 2) / par("pin")[1] * usr.xdiff
  stand.rad.min <- (min.size / 2) / par("pin")[1] * usr.xdiff

  # sorting from large to small ensures that small bubbles will not be hidden
  # behind large bubbles
  events <- events[order(events$Z, decreasing=TRUE), ]

  # determine the size of each circle/legend circle based on the selected type
  type <- match.arg(type)
  switch(type,
         volume = {
           radii     <- stand.rad.min + ((events$Z      / z.max)^(1/3)) * (stand.rad - stand.rad.min)
           radii.leg <- stand.rad.min + ((legend.breaks / z.max)^(1/3)) * (stand.rad - stand.rad.min)
         },
         surface = {
           radii     <- stand.rad.min + sqrt(events$Z      / z.max) * (stand.rad - stand.rad.min)
           radii.leg <- stand.rad.min + sqrt(legend.breaks / z.max) * (stand.rad - stand.rad.min)
         },
         perceptual = {
           # default (if type unspecified)
           radii     <- stand.rad.min + ((events$Z      / z.max)^0.57) * (stand.rad - stand.rad.min)
           radii.leg <- stand.rad.min + ((legend.breaks / z.max)^0.57) * (stand.rad - stand.rad.min)
         }
         )

  # handle multiple colours
  if (is.vector (symbol.bg) && length(symbol.bg) > 1)
    getColour <- colorRamp (symbol.bg)
  else
    # ensure that the function returns X colours when called with X values
    getColour <- function(x) { t(col2rgb(symbol.bg, alpha=TRUE))[rep(1,length(x)), ] }

  # obtain colours for the background (as a matrix)
  bgs <- getColour((events$Z - min(legend.breaks)) / (max(legend.breaks) - min(legend.breaks)))
  if (ncol(bgs) == 3)
    bgs <- cbind(bgs, 255) # add the alpha channel if necessary
  # ... now deal with Z values outside of the range
  if (is.vector (symbol.bg) && length(symbol.bg) > 1)
    {
      outside <- events$Z < min(legend.breaks) | events$Z > max(legend.breaks)
      if (sum(outside) > 0)
        {
          bgs[outside,] <- matrix(c(255,255,255,0), ncol=4)[rep(1,sum(outside)), ]
          warning(sum(outside),
                  " events were outside the legend range and were plotted with",
                  " transparent interiors.  Consider using the addBubbles",
                  " arguments 'legend.breaks', 'min.size', and 'symbol.zero'",
                  " to improve the output.")
        }
    }

  # obtain colours for the legend (as a matrix)
  bgs.leg <- getColour((legend.breaks - min(legend.breaks)) / (max(legend.breaks) - min(legend.breaks)))
  if (ncol(bgs.leg) == 3)
    bgs.leg <- cbind(bgs.leg, 255) # add the alpha channel if necessary

  # convert the matrices to hex values (#RRGGBBAA)
  bgs <- rgb(bgs[,1], bgs[,2], bgs[,3], bgs[,4], maxColorValue=255)
  bgs.leg <- rgb(bgs.leg[,1], bgs.leg[,2], bgs.leg[,3], bgs.leg[,4], maxColorValue=255)

  # compare events$Z to 0; cannot simply use "== 0" given floating-point type
  isZero <- unlist(lapply(events$Z, all.equal, current = 0)) == "TRUE"

  # plot the circles (data with non-zero radii)
  symbols(events$X[!isZero], events$Y[!isZero], circles = radii[!isZero],
          inches = FALSE, bg = bgs[!isZero], fg = symbol.fg, add = TRUE)

  # plot the zero symbol for points (where necessary)
  if (any(isZero) && (!is.logical(symbol.zero) || symbol.zero)) {
    if (is.logical(symbol.zero))
      symbol.zero <- "+" # set to default

    dots <- list(...);
    if (!is.null(dots$pch))
      stop("Specify 'pch' through 'symbol.zero'")

    points(events$X[isZero], events$Y[isZero], pch = symbol.zero, ...)
  }

  # plot the legend if there's a position specified for it
  if (!is.null(legend.pos)) {
    # only plot zero symbol if used
    if (!any(isZero))
      symbol.zero <- FALSE;
    .addBubblesLegend (radii.leg, usr.xdiff, usr.ydiff, symbol.zero, symbol.fg,
                       bgs.leg, legend.pos, legend.breaks, legend.type,
                       legend.title, legend.cex, ...)
  }

  invisible()
}

#==============================================================================
# '...' contains arguments for the 'text()' function
addLabels <- function(data, xlim = NULL, ylim = NULL, polyProps = NULL,
                      placement = "DATA", polys = NULL, rollup = 3,
                      cex = NULL, col = NULL, font = NULL, ...)
{
  # performs the validation below when we know what type of data we have

  # load the limits, so we can clip before plotting
  if (is.null(xlim))
    xlim <- options()$map.xlim;
  if (is.null(ylim))
    ylim <- options()$map.ylim;
  if (is.null(xlim) || is.null(ylim)) {
    stop(
"You must create a plot before you can add to one.\n");
  }

  projectionPlot <- options()$map.projection;
  if (placement == "DATA") {
    projectionData <- attr(data, "projection");
  } else {
    projectionData <- attr(polys, "projection");
  }
  .checkProjection(projectionPlot, projectionData);
  zoneData <- attr(data, "zone");

  # assume PolyData
  EventData <- FALSE;
  type <- "p";

  # test for PolyData/EventData
  if (!any(is.element(names(data), c("PID", "EID")))) {
    stop(
"'data' must be PolyData or EventData.\n");
  }
  # test for 'label' column
  if (!is.element("label", names(data))) {
    stop(
"'data' is missing a column named 'label'.\n");
  }

  # test for EventData
  if (is.EventData (data)) {
    # ensure placement = "DATA"
    if (placement != "DATA") {
      stop(
"When 'data' is EventData, the 'placement' attribute must equal \"DATA\".\n");
    }

    data <- .validateEventData(data);
    if (is.character(data))
      stop(paste("Invalid EventData 'data'.\n", data, sep=""));
    EventData <- TRUE;
    type <- "e";
  }
  # otherwise, PolyData...
  else {
    data <- .validatePolyData(data);
    if (is.character(data))
      stop(paste("Invalid PolyData 'data'.\n", data, sep=""));
  }

  # if placement is 'DATA', test for 'X'/'Y' columns (they might possibly
  # be missing at this point if 'data' is PolyData)
  if (placement == "DATA") {
    if (any(!is.element(c("X", "Y"), names(data)))) {
      stop(
"'data' must contain columns 'X' and 'Y' when placement = \"DATA\"'.\n");
    }

    # at this point, placement is DATA and 'data' has label, X, and Y columns
    # fall through for general handling...

  }
  # otherwise, add 'X'/'Y' columns (only possible for PolyData)
  else {
    if (is.null(polys)) {
      stop(
"Must supply 'polys' when placement is not equal to \"DATA\".\n");
    }
    polys <- .validatePolySet(polys);
    if (is.character(polys))
      stop(paste("Invalid PolySet 'polys'.\n", polys, sep=""));

    # if 'data' contains SID, 'polys' must as well
    if ((is.element("SID", names(data)) && !is.element("SID", names(polys)))) {
      stop(
"Since 'data' contains an SID field, 'polys' must as well.\n");
    }

    # at this point, placement != DATA, 'data' has label column, and 'polys'
    # is a valid PolySet

    # reduce 'polys' if possible to avoid unnecessary processing (we only
    # need to calculate X/Y for PIDs/SIDs that appear in 'data')
    polys <- polys[is.element(paste(polys$PID, polys$SID),
                              unique(paste(data$PID, data$SID))), ];

    # 'polys': PolySet >> PolyData
    if (placement == "CENTROID") {
      polys <- calcCentroid(polys, rollup = rollup);
    } else if (placement == "MEAN_RANGE") {
      polys <- calcMidRange(polys, rollup = rollup);
    } else if (placement == "MEAN_XY") {
      polys <- calcSummary(polys, rollup = rollup, FUN = mean);
    } else {
      stop(paste(
"Unknown 'placement'.  Must be one of \"DATA\", \"CENTROID\",",
"\"MEAN_RANGE\", or \"MEAN_XY\".\n",
                 sep = "\n"));
    }

    # at this point, 'polys' is PolyData with X/Y columns; merge label column
    # from 'data' with PolyData 'polys'

    data <- merge(data,
                  polys,
                  by=intersect(intersect(names(data), names(polys)),
                    c("PID", "SID")),
                  all.x = TRUE);

    # check if any entries were not found in PolyData 'polys'
    n <- nrow(data);
    data <- na.omit(data);
    if (nrow(data) != n) {
      warning(
"Omitted some polygons from 'data' because they did not exist in 'polys'.\n");
    }

    # we can continue as if placement was DATA, as 'data' now has
    # PID, (optional SID), X, Y, and label columns
    # fall through for general handling...

  }

  # clip the data
  data <- data[data$X >= xlim[1] & data$X <= xlim[2] &
               data$Y >= ylim[1] & data$Y <= ylim[2], ];

  # add the labels
  .addFeature(feature = "labels", data = data, polyProps = polyProps,
              isEventData = EventData,
              cex = cex, col = col, font = font, ...);

  # prepare the return value
  if (!is.null(data)) {
    # add projection and zone
    attr(data, "projection") <- projectionData;
    attr(data, "zone") <- zoneData;

    # add to the class attribute
    if (EventData && !is.EventData(data, fullValidation = FALSE)) {
      class(data) <- c("EventData", class(data));
    } else if (!EventData && !is.PolyData(data, fullValidation = FALSE)) {
      class(data) <- c("PolyData", class(data));
    }
  }

  invisible (data);
}

#addLines=========================================2012-03-01
#  Add polylines to an existing map plot.
#  '...' contains arguments for the 'lines()' function
#---------------------------------------------------------NB
addLines <- function(polys, xlim = NULL, ylim = NULL, polyProps = NULL,
                     lty = NULL, col = NULL, arrows = FALSE, ...)
{
  # check 'polys'
  polys <- .validatePolySet(polys);
  if (is.character(polys))
    stop(paste("Invalid PolySet 'polys'.\n", polys, sep=""));

  # load/check 'xlim'/'ylim'
  if (is.null(xlim))
    xlim <- options()$map.xlim;
  if (is.null(ylim))
    ylim <- options()$map.ylim;
  if (is.null(xlim) || is.null(ylim)) {
    stop(
"You must create a plot before you can add to one.\n");
  }

  # load/check 'projection'
  projectionPlot <- options()$map.projection;
  projectionPoly <- attr(polys, "projection");
  .checkProjection(projectionPlot, projectionPoly);

  # clip
  polys <- .clip(polys, xlim, ylim, isPolygons = FALSE, keepExtra = FALSE);
  if (is.null(polys)) {
    warning("All vertices were clipped: no polylines to plot.\n");
    return (NULL);
  }

  # valid columns for PolyProps (described in ?lines under ... and in ?arrows)
  if (arrows)
    pPropsCols <- c("lty", "lwd", "col", "pch", "lend", "ljoin", "lmitre",
                    "angle", "length", "code")
  else
    pPropsCols <- c("lty", "lwd", "col", "pch", "lend", "ljoin", "lmitre",
                    "type")
  # properties that aren't included in par
  nonParProps <- c("angle", "length", "code", "type")

  # validate the polyProps argument
  polyProps <- .validatePolyProps(polyProps, parCols = pPropsCols);
  if (is.character(polyProps))
    stop(paste("Invalid PolyData 'polyProps'.\n", polyProps, sep=""));
  if (is.element("SID", names(polyProps))
      && !is.element("SID", names(polys))) {
    stop(
"No lines to plot since 'polyProps' contains SIDs but 'polys' does not.\n");
  }

  # obtain default values for 'polyProps': !par values from ?lines/?arrows
  propDefaults <- append(par(setdiff(pPropsCols, nonParProps)),      # in par
                         list(type="l",length=0.25,angle=30,code=2)) # !in par
  # don't replace attributes already in 'polyProps' with default values
  propDefaults <- propDefaults[setdiff(names(propDefaults), c(names(polyProps)))]
  # use this function's arguments to override defaults obtained from par
  if (!is.null(col)) propDefaults[["col"]] <- col;
  if (!is.null(lty)) propDefaults[["lty"]] <- lty;
  # look into ..., too, and evaluate statements like
  #   if (!is.null(list(...)$lwd)) propDefaults[[\"lwd\"]] <- list(...)$lwd;
  pPropsColsTmp <- setdiff(pPropsCols, c("col", "lty"));
  expr <- paste("if (!is.null(list(...)$", pPropsColsTmp,
                ")) propDefaults[[\"", pPropsColsTmp,
                "\"]] <- list(...)$", pPropsColsTmp,
                sep="", collapse="; ");
  eval(parse(text=expr))

  # adds an SID column to 'polyProps' if one exists in 'polys'
  polyProps <- .preparePolyProps(polys$PID, polys$SID, polyProps);

  # flesh out 'polyProps' columns (add propDefaults to polyProps, cycling by PID)
  if (length(propDefaults) > 0)
    polyProps <- .addProps(type = "p", polyProps = polyProps, propDefaults);
  polyPropsReturn <- polyProps;

  # create an index for the properties
  # use an expression like
  #   polyProps$props <- paste(polyProps$lty, polyProps$lwd, ..., sep = "-");
  expr <- paste("polyProps$props <- paste(",
                paste("polyProps$", pPropsCols, sep="", collapse=", "),
                ", sep=\"-\")", sep="");
  eval(parse(text=expr));

  # determine if we can use "fast IDs" (i.e., no 'paste')
  fastIDdig <- .createFastIDdig(polysA=polys, polysB=polyProps,
                                cols=c("PID", "SID"));
  # if fastIDdig > 0, we can use "fast IDs"; must divide SID by 10^fastIDdig

  polyProps$IDs <- .createIDs(polyProps, cols=c("PID", "SID"), fastIDdig);
  polyPropsIdx <- split(polyProps$IDs, polyProps$props);

  # create an 'IDs' column for 'polys'
  polys$IDs <- .createIDs(polys, cols=c("PID", "SID"), fastIDdig);

  # reduce before split
  toSplit <- polyProps[!duplicated(polyProps$props), c(pPropsCols, "props")];
  res <- split(toSplit[, pPropsCols], toSplit$props);

  # create lists of X/Y vertices
  polyx <- split(polys$X, polys$IDs);
  polyy <- split(polys$Y, polys$IDs);

  for (prop in names(polyPropsIdx)) {
    toSet <- as.vector(res[prop][[1]]);

    # separate lines with NAs for plotting;
    # for arrows, a start/end of NA does not produce a warning/error
    new.x <- .insertNAs(polyx, polyPropsIdx[[prop]]);
    new.y <- .insertNAs(polyy, polyPropsIdx[[prop]]);

    # when we call lines, set par values like
    #   col = as.vector(toSet$col), lwd = as.vector(toSet$lwd), ...
    if (arrows) {
      eval(parse(text=paste("arrows(",
                   "x0 = new.x[1:(length(new.x)-1)], ",
                   "y0 = new.y[1:(length(new.y)-1)], ",
                   "x1 = new.x[2:length(new.x)], ",
                   "y1 = new.y[2:length(new.y)], ",
                   paste(pPropsCols, " = as.vector(toSet$", pPropsCols, ")",
                         sep="", collapse=", "), ")",
                   sep="")))
    } else {
      eval(parse(text=paste("lines(x = new.x, y = new.y, ",
                   paste(pPropsCols, " = as.vector(toSet$", pPropsCols, ")",
                         sep="", collapse=", "), ")",
                   sep="")))
    }
  }

  # add to the class attribute
  if (!is.null(polyPropsReturn) &&
      !is.PolyData(polyPropsReturn, fullValidation = FALSE))
    class(polyPropsReturn) <- c("PolyData", class(polyPropsReturn));

  invisible(polyPropsReturn);
}

#addPoints========================================2010-11-10
#  Add points to an existing map plot.
#  '...' contains arguments for the '.addFeature()' function,
#  which will then pass those arguments to 'points' function.
#---------------------------------------------------------NB
addPoints <- function(data, xlim = NULL, ylim = NULL, polyProps = NULL,
                      cex = NULL, col = NULL, pch = NULL, ...)
{
  # performs the validation below when we know what type of data we have

  if (is.null(xlim))
    xlim <- options()$map.xlim;
  if (is.null(ylim))
    ylim <- options()$map.ylim;
  if (is.null(xlim) || is.null(ylim)) {
    stop(
"You must create a plot before you can add to one.\n");
  }

  projectionPlot <- options()$map.projection;
  projectionData <- attr(data, "projection");
  .checkProjection(projectionPlot, projectionData);

  # perform some additional validation
  # convert matrix to data frame
  if (is.matrix(data)) {
    data <- .mat2df(data);
  }
  if (!any(is.element(names(data), c("PID", "EID")))) {
    stop("'data' must be either PolyData or EventData.\n");
  }

  # if it contains neither an EID nor PID column, add an EID column
  # and produce a warning
  if (!any(is.element(c("EID", "PID"), names(data)))) {
    data$EID <- 1:nrow(data);
    warning(paste(
"'data' contains neither an 'EID' nor 'PID' column.  Added an 'EID' column",
"for plotting only.\n",
                  sep = "\n"));
  }

  # is it event data? -- validate it here
  isEventData <- FALSE;
  type <- "p";
  if (is.element("EID", names(data)) &&
      !is.element("PID", names(data))) {
    isEventData <- TRUE;
    type <- "e";
    data <- .validateEventData(data);
    if (is.character(data))
      stop(paste("Invalid EventData 'data'.\n", data, sep=""));
    # we know the EventData has X and Y columns
  } else {
    data <- .validatePolyData(data);
    if (is.character(data))
      stop(paste("Invalid PolyData 'data'.\n", data, sep=""));
    # does the PolyData have X and Y columns?
    if (any(!is.element(c("X", "Y"), names(data))))
      stop(
"'data' must contain columns 'X' and 'Y'.\n");
  }

  # clip
  data <- data[data$X >= xlim[1] & data$X <= xlim[2]
               & data$Y >= ylim[1] & data$Y <= ylim[2], ];

  polyPropsReturn <-
    .addFeature(feature = "points", data = data, polyProps = polyProps,
                isEventData = isEventData,
                cex = cex, col = col, pch = pch, ...);

  if (!is.null(polyPropsReturn) &&
      !is.PolyData(polyPropsReturn, fullValidation = FALSE))
    class(polyPropsReturn) <- c("PolyData", class(polyPropsReturn));

  invisible (polyPropsReturn);
}

#addPolys=========================================2012-03-01
#  Add polygons to an existing map plot.
#---------------------------------------------------------NB
addPolys <- function(polys, xlim = NULL, ylim = NULL, polyProps = NULL,
                     border = NULL, lty = NULL, col = NULL, colHoles = NULL,
                     density = NA, angle = NULL, ...)
{
  # check 'polys'
  polys <- .validatePolySet(polys);
  if (is.character(polys))
    stop(paste("Invalid PolySet 'polys'.\n", polys, sep=""));

  # load/check 'xlim'/'ylim'
  if (is.null(xlim))
    xlim <- options()$map.xlim;
  if (is.null(ylim))
    ylim <- options()$map.ylim;
  if (is.null(xlim) || is.null(ylim)) {
    stop(
"You must create a plot before you can add to one.\n");
  }

  # load/check 'projection'
  projectionPlot <- options()$map.projection;
  projectionPoly <- attr(polys, "projection");
  .checkProjection(projectionPlot, projectionPoly);

  # clip
  polys <- .clip(polys, xlim, ylim, isPolygons = TRUE, keepExtra = FALSE);
  if (is.null(polys)) {
    warning("All vertices were clipped: no polygons to plot.\n");
    return (NULL);
  }

  # "polyProps": valid columns
  pPropsCols <- c("angle", "border", "col", "colHoles", "lty", "density")

  # validate the polyProps argument
  pProps <- .validatePolyProps(polyProps, parCols=pPropsCols);
  if (is.character(pProps))
    stop(paste("Invalid PolyData 'polyProps'.\n", pProps, sep=""));
  if (is.element("SID", names(pProps))
      && !is.element("SID", names(polys))) {
    stop(
"No polygons to plot since 'polyProps' contains SIDs but 'polys' does not.\n");
  }

  # default values
  pPropsDefs <- list(angle = 45, border = 1, col = 0, colHoles = NULL,
                     lty = par("lty"), density = NA);
  # ... used only when not already in "pProps"
  pPropsDefs <- pPropsDefs[!is.element(names(pPropsDefs), names(pProps))];

  # allow arguments to override defaults
  for (p in setdiff(pPropsCols, "density")) {
    # create an expression like
    #   if (!is.null(angle)) pPropsDefaults[["angle"]] <- angle;
    # for each of the properties and evaluate it
    expr <- paste("if (!is.null(", p, ")) pPropsDefs[[\"", p, "\"]] <- ", p,
                  sep="");
    eval(parse(text=expr));
  }
  # handle "density" with care: in R (1.8.1+), density == NULL is equivalent to
  # density == 0, and density == NA equivalent to density < 0; make the
  # conversion here to match S-PLUS
  if (is.null(density)) {
    pPropsDefs[["density"]] <- 0;
  } else if (is.na(density)) {
    pPropsDefs[["density"]] <- -1;
  } else {
    pPropsDefs[["density"]] <- density;
  }

  # add an SID column to 'polyProps' (if one exists in 'polys')
  pProps <- .preparePolyProps(polys$PID, polys$SID, pProps);

  # add all of the new attributes to 'polyProps'; if "colHoles" not already
  # in "pProps", this function will not add it (since def. value is NULL)
  if (length(pPropsDefs) > 0)
    pProps <- .addProps(type = "p", polyProps = pProps, pPropsDefs);

  # only play with "pProps$colHoles" if it's a character vector; otherwise, cannot possibly
  # contain the string "transparent"
  if (is.element("colHoles", names(pProps)) && is.character (pProps$colHoles)) {
    pProps$colHoles[pProps$colHoles == "transparent"] <- NA;
  }

  # save a copy to return
  polyPropsReturn <- pProps;

  # create an index for the properties
  # use an expression like
  #   pProps$props <- paste(pProps$angle, pProps$border, ..., sep = "-");
  # that includes each of the properties and evaluate it
  expr <- paste("pProps$props <- paste(",
                paste("pProps$", pPropsCols, sep="", collapse=", "),
                ", sep=\"-\")", sep="");
  eval(parse(text=expr));

  # attempt to use fast IDs (i.e., no "paste"); calculate it here so we can
  # reuse it and also use consistent indices between calls to .createIDs
  fastIDdig <- .createFastIDdig(polysA=polys, polysB=pProps,
                                cols=c("PID", "SID"));

  #----------------------------------------------------------------------------
  # Create a data structure that links each plotting style with a list
  # of IDs to plot using that style.
  polys$IDs <- .createIDs(polys, cols=c("PID", "SID"), fastIDdig);
  idxS <- which(!duplicated(polys$IDs))
  idxE <- c((idxS-1)[-1], length(polys$IDs))
  holes <- (polys$POS[idxS] > polys$POS[idxE])
  holesIDs <- (polys[idxS, "IDs"])[holes];

  pProps$IDs <- .createIDs(pProps, cols=c("PID", "SID"), fastIDdig)
  holes <- is.element(pProps$IDs, holesIDs)

  pPropsOuter <- pProps[!holes, c("IDs", "props")]
  pPropsOuter <- split(pPropsOuter$IDs, pPropsOuter$props)

  # the only case where we (completely) do not plot holes as when
  # all of the "colHole"s are trasparent
  if (is.element("colHoles", names(pProps)) &&
      all(is.na(pProps$colHoles))) {
    pPropsHoles <- NA
  } else if (any(holes)) {
    pPropsHoles <- pProps[holes, c("IDs", "props")]
    pPropsHoles <- split(pPropsHoles$IDs, pPropsHoles$props)
  } else {
    pPropsHoles <- NULL
  }

  integrateHoles <- !is.element("colHoles", names(pProps));
  #----------------------------------------------------------------------------

  #----------------------------------------------------------------------------
  # Create a data structure that links each plotting style with the style's
  # associated settings.
  pProps <- pProps[!duplicated(pProps$props), ];
  pProps <- split(pProps[, intersect(pPropsCols, names(pProps))], pProps$props);
  #----------------------------------------------------------------------------

  # create 'polylines' and 'polys'
  polylines <- polys;

  # "integrate holes" for 'polys' if necessary
  # note that "polylines" is not "rolled up"; it contains original PIDs/SIDs
  if (integrateHoles && length(unique(polys$SID)) > 1) {
    polys <- .rollupPolys(polys, rollupMode = 2, exteriorCCW = -1,
                          closedPolys = 1, addRetrace = TRUE);
  }

  # create an 'IDs' column for both 'polys' and 'polylines'
  polylines$IDs <- .createIDs(polylines, cols=c("PID", "SID"), fastIDdig);
  polys$IDs <- .createIDs(polys, cols=c("PID", "SID"), fastIDdig);

  polylinex <- split(polylines$X, polylines$IDs);
  polyliney <- split(polylines$Y, polylines$IDs);
  polyx <- split(polys$X, polys$IDs);
  polyy <- split(polys$Y, polys$IDs);

  #- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # WORK-AROUND:
  # Problem: the R polygon command behaves different from the S polygon
  #   command.  The R command draws filling lines according to the 'lty' par
  #   parameter, whereas the S polygon command ignores 'lty' for filling
  #   lines.
  #- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Problem:
  #   S-PLUS draws retrace lines when 'lty' != 0.
  #   R refused to draw shading lines (i.e., density) when 'lty' == 0.
  # To stop S-PLUS from plotting retrace lines, the polygon command needs:
  # 1) lty = 0 *AND*
  # 2) border = F
  # To stop R from plotting retrace lines, the polygon command needs:
  # 1) lty = 0 *OR*
  # 2) border = F
  #- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (!is.null(version$language) && (version$language == "R")) {
    fillBorder <- NA;   # in R, NA => omit borders
    polyLty <- 1;
  } else {
    fillBorder <- FALSE;
    polyLty <- 0;
  }

  # walk through the properties, plotting the appropriate fills and edges
  # as we go
  mode <- "o"
  pPropsIdx <- pPropsOuter

  # A) if "pPropsHoles" == NA, all holes were transparent and we have
  #    absolutely nothing to plot for holes; in this instance, range need to
  #    include the holes
  # B) if "pPropsHoles" != NA, we either plot borders (if we integrated holes)
  #    or a filled polygon; in these instances, range much include the holes

  # pPropsHoles is NULL when no holes
  # pPropsHoles is NA when "colHoles" are all transparent
  if (!is.null(pPropsHoles) && !any(is.na(pPropsHoles))) {
    range <- c(1:length(pPropsOuter),NA,1:length(pPropsHoles))
  } else {
    range <- c(1:length(pPropsOuter),NA,NA)
  }

  for (propIdx in range) {
    if (is.na(propIdx)) {
      if (mode == "o") {
        # switch modes
        mode <- "h"
        pPropsIdx <- pPropsHoles
        if (is.null(pPropsIdx))
          break
      }
      next
    }

    # set up the drawing parameters
    prop <- names(pPropsIdx)[propIdx]
    toSet <- as.vector(pProps[prop][[1]]);

    # need to as.vector function to remove levels (remove the factors!)
    border <- as.vector(toSet$border);
    lty <- as.vector(toSet$lty);
    if (mode == "o" || is.null(toSet$colHoles)) {
      col <- as.vector(toSet$col);
    } else if (mode == "h") {
      col <- as.vector(toSet$colHoles);
    }
    density <-  as.vector(toSet$density);
    angle <- as.vector(toSet$angle);

    if (!is.na(col)) {
      # plot the fill
      new.polyx <- .insertNAs(polyx, pPropsIdx[[propIdx]])
      new.polyy <- .insertNAs(polyy, pPropsIdx[[propIdx]])

      # Subtle differences exist between the S-PLUS and R "polygon" command.
      # The S-PLUS command fails if not given at least three points, whereas the
      # R command draws nothing for 0-1 point and a line for 2 points.
      # The following tests for 0 points...
      if (!all(is.na(new.polyx)))
        polygon(x = new.polyx, y = new.polyy, col = col, lty = polyLty,
                border = fillBorder, density = density, angle = angle, ...);

      # plot the edges
      if (is.na(as.logical(border)) || as.logical(border)) {
        if (is.character(border) && !is.na(as.logical(border)))
          border <- as.logical(border);

        new.polylinex <- .insertNAs(polylinex, pPropsIdx[[propIdx]])
        new.polyliney <- .insertNAs(polyliney, pPropsIdx[[propIdx]])

        polygon(x = new.polylinex, y = new.polyliney, lty = lty,
                border = border, density = 0, ...);
      }
    }
  }

  if (!is.null(polyPropsReturn) &&
      !is.PolyData(polyPropsReturn, fullValidation = FALSE))
    class(polyPropsReturn) <- c("PolyData", class(polyPropsReturn));

  invisible(polyPropsReturn);
}

#==============================================================================
# side = -1 implies outside
# side =  0 implies both sides
# side = +1 implies inside
# Note:
#   It is up to the caller to replot the polygons if the stipples cover part
#   of it.
# 'distance' measured as a percentage of the abs(diff(xlim))
addStipples <- function(polys, xlim = NULL, ylim = NULL, polyProps = NULL,
                        side = 1, density = 1, distance = 4, ...)
{
  # initialization
  polys <- .validatePolySet(polys);
  if (is.character(polys))
    stop(paste("Invalid PolySet 'polys'.\n", polys, sep=""));

  if (is.null(xlim))
    xlim <- options()$map.xlim;
  if (is.null(ylim))
    ylim <- options()$map.ylim;
  if (is.null(xlim) || is.null(ylim)) {
    stop(
"You must create a plot before you can add to one.\n");
  }

  projectionPlot <- options()$map.projection;
  projectionPoly <- attr(polys, "projection");
  .checkProjection(projectionPlot, projectionPoly);

  # always thicken polys to a distance 2% of the plot region (xlim)
  # since thickenPolys() wants a tolerance in kilometers, convert degrees
  # to kilometers if necessary
  if (projectionPlot == "LL") {
    # Equatorial radius 6,378.14 km
    # Polar radius 6,356.78 km
    # Mean radius 6,371.3 km
    # Sources:
    #   http://en.wikipedia.org/wiki/Earth
    #   http://en.wikipedia.org/wiki/Earth_radius
    R <- 6371.3;

    # use arc length along the line of latitude (mean(ylim)) (NOT a Great
    # Circle distance in this case)

    # this is the radius of the line of latitude, perpedicular to the axis
    # of rotation
    radiusAtLat <- R * sin(pi/2.0 - (abs(mean(ylim))*pi/180.0));
    arcLen <- radiusAtLat * (abs(diff(xlim))*pi/180.0);
    tol <- arcLen * 0.02;
  } else {
    tol <- abs(diff(xlim)) * 0.02;
  }

  # distance is measured as a percentage abs(diff(xlim)); adjust for
  # latitude
  distance <- distance / 100.0;
  distX <- abs(diff(xlim)) * distance;
  if (projectionPlot == "LL") {
    xyRatio <- cos((mean(ylim) * pi) / 180);
    distY <- distX * xyRatio;
  } else {
    distY <- distX;
  }

  # 1) clip polys to the window to avoid excess computations
  # 2) extend the limits to avoid plotting stipples along the border where
  #    a border only artificially exists because of the clipping; 1.1 fudge-
  #    factor added to ensure points will be off the plot region
  polys <- .clip(polys,
                 xlim=xlim + c(-distX, distX)*1.1,
                 ylim=ylim + c(-distY, distY)*1.1,
                 isPolygons = TRUE, keepExtra = FALSE);
  thick <- thickenPolys(polys, tol = tol, keepOrig = FALSE);
  x <- rep(thick$X, each=(15 * density));
  y <- rep(thick$Y, each=(15 * density));
  events <- data.frame(EID = 1:length(x),
                       X = x + runif(length(x), -distX, distX),
                       Y = y + runif(length(y), -distY, distY));

  # keep only those events that fall within the proper limits
  events <- events[events$X > xlim[1] & events$X < xlim[2]
                   & events$Y > ylim[1] & events$Y < ylim[2], ];

  link <- findPolys(events, polys)
  if (side == -1) {
    # if outside
    stip <- events[!is.element(events$EID, unique(link$EID)), ]
  } else if (side == 1) {
    # if inside
    link <- link[link$Bdry == FALSE, ]
    stip <- events[is.element(events$EID, unique(link$EID)), ]
  } else {
    # else both sides
    stip <- events;
  }

  stip <- merge(stip, link[, intersect(names(link), c("EID", "PID", "SID"))],
                all.x = TRUE, by = "EID");
  stip <- stip[, !is.element(names(stip), "EID")];

  .addFeature(feature = "points", data = stip, polyProps = polyProps,
              isEventData = FALSE, ...);

  if (!is.null(polyProps) &&
      !is.PolyData(polyProps, fullValidation = FALSE))
    class(polyProps) <- c("PolyData", class(polyProps));
  invisible(polyProps);
}

#==============================================================================
appendPolys <- function(polys, mat, PID = NULL, SID = NULL, isHole = FALSE)
{
  # validate the arugments
  if (missing(polys) || missing(mat)) {
    stop(paste(
"Requires arguments 'polys' and 'mat'.  To create a new PolySet, pass",
"NULL as the 'polys' argument.\n",
               sep = "\n"));
  }
  # validate 'polys'
  if (!is.null(polys)) {
    polys <- .validatePolySet(polys);
    if (is.character(polys))
      stop(paste("Invalid PolySet 'polys'.\n", polys, sep=""));
  }
  # validate 'mat'
  if (!is.matrix(mat) || (ncol(mat) != 2) || !is.numeric(mat)) {
    stop(paste(
"'mat' must be a matrix with two numeric columns.  The first column contains",
"the X values and the second contains the Y values.\n",
               sep = "\n"));
  }
  if ((!is.null(PID) && (!is.numeric(PID) || length(PID) != 1)) ||
      (!is.null(SID) && (!is.numeric(SID) || length(SID) != 1)) ||
      (!is.logical(isHole) || length(isHole) != 1)) {
    stop(paste(
"If supplied as arguments, the length of 'PID', 'SID', and 'isHole' must",
"equal one.  'PID' and 'SID' must be numerics, and 'isHole' must be a",
"logical.\n",
               sep = "\n"));
  }

  # save the attributes of the data frame (.validatePolySet returns a data
  # frame); in this case, SAVE ALL ATTRIBUTES since we are "modifying" an
  # existing data set
  attrNames <- setdiff(names(attributes(polys)),
                       c("names", "row.names", "class"));
  attrValues <- attributes(polys)[attrNames];

  # get PID
  if (is.null(PID)) {
    if (is.null(polys)) {
      PID <- 1;
    } else {
      PID <- max(polys$PID) + 1;
    }
  }

  # get SID; only set it if there exists an SID column
  if (is.null(SID) && !is.null(polys) && is.element("SID", names(polys))) {
    curSIDs <- polys[polys$PID == PID, "SID"];
    if (!is.null(curSIDs) && (length(curSIDs) > 0)) {
      SID <- max(curSIDs) + 1;
    } else {
      SID <- 1;
    }
  }

  # build the POS column
  POS <- 1:nrow(mat);
  if (isHole)
    POS <- rev(POS);

  # build the new rows
  if (!is.null(SID)) {
    newRows <- data.frame(PID = rep(PID, nrow(mat)),
                          SID = rep(SID, nrow(mat)),
                          POS = POS, X = mat[, 1], Y = mat[, 2]);
  } else {
    newRows <- data.frame(PID = rep(PID, nrow(mat)),
                          POS = POS, X = mat[, 1], Y = mat[, 2]);
  }

  # if SID not in 'polys', but it is specified here, add it to 'polys'
  if (!is.null(polys) && !is.element("SID", names(polys)) && !is.null(SID))
    polys$SID <- rep(1, nrow(polys));

  # create ordering for 'poly' columns
  extraCols <- vector();
  if (!is.null(polys)) {
    extraCols <- setdiff(names(polys), c("PID", "SID", "POS", "X", "Y"));
    cols <- c("PID", "POS", "X", "Y", extraCols);
  } else {
    cols <- c("PID", "POS", "X", "Y");
  }
  if (is.element("SID", names(polys)) || !is.null(SID))
    cols <- c(cols[1], "SID", cols[-1]);

  # append NA for extraCols in our 'newRows'
  if (length(extraCols) >  0)
    newRows <- cbind(newRows, matrix(NA, nrow=nrow(newRows),
                                     ncol=length(extraCols)));

  # bind 'em
  if (!is.null(polys)) {
    ret <- (rbind(polys[, cols], newRows));
  } else {
    ret <- (newRows);
  }

  if (!is.null(ret)) {
    # restore the attributes (including projection and zone)
    attributes(ret) <- c(attributes(ret), attrValues);

    if (!is.PolySet(ret, fullValidation = FALSE))
      class(ret) <- c("PolySet", class(ret));
  }

  return (ret);
}

#==============================================================================
as.EventData <- function(x, projection = NULL, zone = NULL)
{
  # if a data frame, remove all other classes
  if (is.data.frame(x))
    class(x) <- "data.frame"

  x <- .validateEventData(x);
  if (is.character(x))
    stop(paste("Cannot coerce into EventData.\n", x, sep=""));;

  if (!is.EventData(x, fullValidation = FALSE))
    class(x) <- c("EventData", class(x));
  if (!is.null(projection))
    attr(x, "projection") <- projection;
  if (!is.null(zone))
    attr(x, "zone") <- zone;

  return (x);
}

#==============================================================================
as.LocationSet <- function(x)
{
  # if a data frame, remove all other classes
  if (is.data.frame(x))
    class(x) <- "data.frame"

  x <- .validateLocationSet(x);
  if (is.character(x))
    stop(paste("Cannot coerce into a LocationSet.\n", x, sep=""));;

  if (!is.LocationSet(x, fullValidation = FALSE))
    class(x) <- c("LocationSet", class(x));

  return (x);
}

#==============================================================================
as.PolyData <- function(x, projection = NULL, zone = NULL)
{
  # if a data frame, remove all other classes
  if (is.data.frame(x))
    class(x) <- "data.frame"

  x <- .validatePolyData(x);
  if (is.character(x))
    stop(paste("Cannot coerce into PolyData.\n", x, sep=""));;

  if (!is.PolyData(x, fullValidation = FALSE))
    class(x) <- c("PolyData", class(x));
  if (!is.null(projection))
    attr(x, "projection") <- projection;
  if (!is.null(zone))
    attr(x, "zone") <- zone;

  return (x);
}

#==============================================================================
as.PolySet <- function(x, projection = NULL, zone = NULL)
{
  # if a data frame, remove all other classes
  if (is.data.frame(x))
    class(x) <- "data.frame"

  x <- .validatePolySet(x);
  if (is.character(x))
    stop(paste("Cannot coerce into a PolySet.\n", x, sep=""));;

  if (!is.PolySet(x, fullValidation = FALSE))
    class(x) <- c("PolySet", class(x));
  if (!is.null(projection))
    attr(x, "projection") <- projection;
  if (!is.null(zone))
    attr(x, "zone") <- zone;

  return (x);
}

#==============================================================================
calcArea <- function(polys, rollup = 3)
{
  polys <- .validatePolySet(polys);
  if (is.character(polys))
    stop(paste("Invalid PolySet 'polys'.\n", polys, sep=""));

  # copied the following code to calcCentroid()
  msg <- ""
  # check if projection is LL
  if (!is.null(attr(polys, "projection"))
      && !is.na(attr(polys, "projection"))
      && attr(polys, "projection") == "LL") {
    msg <- paste(
"To calculate the areas of polygons defined by longitude-latitude",
"coordinates, this routine first projects them in UTM.\n",
                  sep = "\n");
    # check if zone needs to be set
    # SIMILAR code to convUL()
    if (is.null(attr(polys, "zone")) || is.na(attr(polys, "zone"))) {
      m <- mean(polys$X)
      if (m <= -180 || m > 180) {
        msg <- paste(msg, paste(
"\nAttempted to automatically calculate the missing 'zone' attribute, but",
"that failed because the mean longitude falls outside the range",
"-180 < x <= 180.\n",
                                sep = "\n"),
                     sep = "");
        stop(msg);
      }
      attr(polys, "zone") <- ceiling((m + 180) / 6);
      msg <- paste(msg, paste(
"\nFor the UTM conversion, automatically detected zone ", attr(polys, "zone"),
".\n",
                              sep = ""),
                   sep = "");
    } else {
      msg <- paste(msg, paste(
"\nUsing 'zone' attribute (", attr(polys, "zone"), ").\n", sep = ""),
                   sep = "");
    }
    warning(msg);
    # now the zone is set
    polys <- convUL(polys);
  }

  # rollup 'polys' to fix CW/CCW and rollup polys
  # since the calcArea function returns signed areas, this is important
  polys <- .rollupPolys(polys, rollupMode = rollup, exteriorCCW = 1,
                        closedPolys = 1, addRetrace = TRUE);

  inRows <- nrow(polys);
  outCapacity <- inRows;

  # create the data structure that the C functions expect
  if (!is.element("SID", names(polys))) {
    inID <- c(polys$PID, integer(length = inRows), polys$POS);
  } else {
    inID <- c(polys$PID, polys$SID, polys$POS);
  }
  inXY <- c(polys$X, polys$Y);

  # call the C function
  results <- .C("calcArea",
                inID = as.integer(inID),
                inXY = as.double(inXY),
                inVerts = as.integer(inRows),
                outID = integer(2 * outCapacity),
                outArea = double(outCapacity),
                outRows = as.integer(outCapacity),
                outStatus = integer(1),
                PACKAGE = "PBSmapping");
  # note: outRows is set to how much space is allocated -- the C function
  #       should consider this

  if (results$outStatus == 1) {
    stop(
"Insufficient physical memory for processing.\n");
  }
  if (results$outStatus == 2) {
    stop(paste(
"Insufficient memory allocated for output.  Please upgrade to the latest",
"version of the software, and if that does not fix this problem, please",
"file a bug report.\n",
               sep = "\n"));
  }

  # determine the number of rows in the result
  outRows <- as.vector(results$outRows);

  # extract the data from the C function results
  if (outRows > 0) {
    d <- data.frame(PID = results$outID[1:outRows],
                    SID = results$outID[(outCapacity+1):(outCapacity+outRows)],
                    area = results$outArea[1:outRows]);

    if (!is.element("SID", names(polys)) || rollup == 1)
      d$SID <- NULL;

    if (!is.PolyData(d, fullValidation = FALSE))
      class(d) <- c("PolyData", class(d));

    return(d);
  } else {
    return(NULL);
  }
}

#------------------------------------------------------------------------------
# calcCentroid:
#
# Capabilities:
#   [Y] Holes
#   [-] Projection
# Legend:
#   [Y] yes   [N] no   [-] N/A   [?] good question
#
# Returns:
#   PolyData with centroids of all polygons.
#   Values for 'rollup' include:
#     1 - rollup to the PID level
#     2 - rollup to the hole level; i.e., integrate holes
#     3 - do not rollup
#
#   Setting onlyPID = TRUE results in only one centroid per PID.
#------------------------------------------------------------------------------
calcCentroid <- function(polys, rollup = 3)
{
  polys <- .validatePolySet(polys);
  if (is.character(polys))
    stop(paste("Invalid PolySet 'polys'.\n", polys, sep=""));

  # to ensure correct orientation of the exterior/interior vertices, run a
  # fixPOS() on the PolySet
  if ((rollup == 1) || (rollup == 2)) {
    polys <- .rollupPolys(polys, rollupMode = rollup, exteriorCCW = 1,
                          closedPolys = 1, addRetrace = TRUE);
    # As for the overall clockwise/counter-clockwise orientation of the
    # merged PolySet, it shouldn't matter.  The retrace lines don't change
    # the area, and hence shouldn't change the centroid.
  }

  inRows <- nrow(polys);
  outCapacity <- inRows;

  # create the data structure that the C functions expect
  if (!is.element("SID", names(polys))) {
    inID <- c(polys$PID, integer(length = inRows), polys$POS);
  } else {
    inID <- c(polys$PID, polys$SID, polys$POS);
  }
  inXY <- c(polys$X, polys$Y);

  # call the C function
  results <- .C("calcCentroid",
                inID = as.integer(inID),
                inXY = as.double(inXY),
                inVerts = as.integer(inRows),
                outID = integer(2 * outCapacity),
                outXY = double(2 * outCapacity),
                outRows = as.integer(outCapacity),
                outStatus = integer(1),
                PACKAGE = "PBSmapping");
  # note: outRows is set to how much space is allocated -- the C function
  #       should consider this

  if (results$outStatus == 1) {
    stop(
"Insufficient physical memory for processing.\n");
  }
  if (results$outStatus == 2) {
    stop(paste(
"Insufficient memory allocated for output.  Please upgrade to the latest",
"version of the software, and if that does not fix this problem, please",
"file a bug report.\n",
               sep = "\n"));
  }

  # determine the number of rows in the result
  outRows <- as.vector(results$outRows);

  # extract the data from the C function results
  if (outRows > 0) {
    d <- data.frame(PID = results$outID[1:outRows],
                    SID = results$outID[(outCapacity+1):(outCapacity+outRows)],
                    X = results$outXY[1:outRows],
                    Y = results$outXY[(outCapacity+1):(outCapacity+outRows)]);

    if (!is.element("SID", names(polys)))
      d$SID <- NULL;

    # add projection and zone
    attr(d, "projection") <- attr(polys, "projection");
    attr(d, "zone") <- attr(polys, "zone");

    if (!is.PolyData(d, fullValidation = FALSE))
      class(d) <- c("PolyData", class(d));

    return(d);
  } else {
    return(NULL);
  }
}

#==============================================================================
# New calcConvexHull using chull() in package grDevices
calcConvexHull <- function(xydata,keepExtra=FALSE) {
	xydata <- .validateXYData(xydata);
	if (is.character(xydata))
		stop(paste("Invalid X/Y data 'xydata'.\n", xydata, sep=""));

	inRows <- nrow(xydata);
	ii     <- grDevices::chull(xydata$X,xydata$Y); n <- length(ii)
	hull0  <- xydata[ii,]
	zXY    <- is.element(names(hull0),c("X","Y"))
	zExtra <- !zXY
	hull   <- cbind(PID=rep(1,n),POS=1:n,hull0[,zXY])
	if (keepExtra && any(zExtra==TRUE)) {
			hull <- cbind(hull,hull0[,zExtra])
			names(hull)[5:ncol(hull)] <- names(hull0)[zExtra] }

	# add projection and zone attributes
	attr(hull, "inRows") <- inRows;
	attr(hull, "projection") <- attr(xydata, "projection");
	attr(hull, "zone") <- attr(xydata, "zone");

	if (!is.PolySet(hull, fullValidation = FALSE))
		class(hull) <- c("PolySet", class(hull));
	hull <- rbind(hull,hull[1,]); hull[n+1,2] <- n+1; # close polygon
	return(hull)
}

#==============================================================================
calcLength <- function (polys, rollup = 3, close = FALSE)
{
  polys <- .validatePolySet(polys);
  if (is.character(polys))
    stop(paste("Invalid PolySet 'polys'.\n", polys, sep=""));

  # according to documentation
  if (rollup == 2) {
    rollup <- 3;
    warning(
"'rollup = 2' has no meaning in 'calcLength'; reset it to 'rollup = 3'.\n");
  }

  # STEP 1: compute a distance column
  if (!is.null(attr(polys, "projection"))
      && !is.na(attr(polys, "projection"))
      && ((attr(polys, "projection") == "UTM")
          || (attr(polys, "projection") == "LL")
          || (attr(polys, "projection") == 1))) {
    # if necessary, close polygons; don't roll them up until after
    # the distance calculation
    if (close) {
      polys <- .rollupPolys(polys, rollupMode = 3, exteriorCCW = -1,
                            closedPolys = 1, addRetrace = FALSE);
    }

    polys$D <- .calcDist(polys);
  }
  else {
    stop(paste(
"Invalid projection attribute.  Supported projections include \"LL\",",
"\"UTM\", and 1.\n"));
  }

  # STEP 2: prepare the result
  if (is.element("D", names(polys))) {
    # simplify polys
    polys <- polys[, which(is.element(names(polys), c("PID", "SID", "D")))]

    # build index and remove "unnecessary" distances (i.e., distances
    # between polygons)
    polys$IDX <- .createIDs(polys, cols=c("PID", "SID"));
    polys$D[which(!duplicated(polys$IDX)[c(2:nrow(polys), 1)])] <- 0;

    # split and sum
    # if rollup equals 1, split by PID rather than IDX
    if (rollup == 1) {
      idx <- "PID";
    } else {
      idx <- "IDX";
    }
    result <- split(polys$D, polys[[idx]])
    result <- lapply(result, sum);
    polys <- polys[!duplicated(polys[[idx]]), ]
    toMerge <- data.frame(names(result), length=unlist(result))
    names(toMerge)[1] <- idx;
    polys <- merge(polys, toMerge, by=idx)
  }

  if (is.element("SID", names(polys)) && (rollup != 1)) {
    d <- data.frame(PID = polys$PID,
                    SID = polys$SID,
                    length = polys$length);
  } else {
    d <- data.frame(PID = polys$PID,
                    length = polys$length);
  }

  if (!is.null(d) &&
      !is.PolyData(d, fullValidation = FALSE))
    class(d) <- c("PolyData", class(d));

  return (d);
}

#==============================================================================
calcMidRange <- function(polys, rollup = 3)
{
  polys <- .validatePolySet(polys);
  if (is.character(polys))
    stop(paste("Invalid PolySet 'polys'.\n", polys, sep=""));

  # calculate range
  polys <- calcSummary (polys, FUN="range", rollup = rollup);

  # calculate mean range
  polys$X <- (polys$X1 + polys$X2) / 2;
  polys$Y <- (polys$Y1 + polys$Y2) / 2;

  # drop obsolete columns
  polys$X1 <- NULL;
  polys$X2 <- NULL;
  polys$Y1 <- NULL;
  polys$Y2 <- NULL;

  # no need to assign 'projection' and 'zone' attributes; 'calcSummary'
  # maintains them

  # assign class if necessary
  if (!is.null(polys) &&
      !is.PolyData(polys, fullValidation = FALSE))
    class(polys) <- c("PolyData", class(polys));

  return (polys);
}

#calcSummary======================================2012-03-01
#  Apply functions to polygons in a PolySet.
#---------------------------------------------------------NB
calcSummary <- function(polys, rollup = 3, FUN, ...)
  # Maintains extra attributes.
{
  polys <- .validatePolySet(polys);
  if (is.character(polys))
    stop(paste("Invalid PolySet 'polys'.\n", polys, sep=""));

  # roll up the PolySet if necessary; don't change the number of points
  if (rollup == 1 || rollup == 2) {
    polys <- .rollupPolys(polys, rollupMode = rollup, exteriorCCW = -1,
                          closedPolys = -1, addRetrace = FALSE);
  } else if (rollup != 3) {
    stop(
"Invalid value for argument 'rollup'.  Must equal 1, 2 or 3.\n");
  }

  # create an index
  idx <- .createIDs(polys, cols=c("PID", "SID"), fastIDdig = NULL);
  uniqueIdx <- unique(idx);

  # get data value (run FUN on each vector of X and vector of Y)
  dataX <- lapply(split(polys$X, idx), FUN, ...);
  dataY <- lapply(split(polys$Y, idx), FUN, ...);

  nCol <- length(dataX[[1]]);

  # as soon as you bind in a string with numbers in a matrix, the whole
  # matrix becomes characters (all items in matrix have same type)

  # cbind() can result in calling the data frame method if the arguments
  # are all either data frames or vectors, and this will result in the
  # coercion of character vectors to factors. (source: R cbind() help file)

  dataX <- cbind(data.frame(IDX=names(dataX)),
                 as.data.frame(matrix(unlist(dataX, use.names = FALSE),
                                      ncol = nCol, byrow = TRUE)));
  dataX$IDX <- as.character(dataX$IDX); # remove factor
  dataY <- cbind(data.frame(IDX=names(dataY)),
                 as.data.frame(matrix(unlist(dataY, use.names = FALSE),
                                      ncol = nCol, byrow = TRUE)));
  dataY$IDX <- as.character(dataY$IDX); # remove factor

  if (nCol > 1) {
    names(dataX) <- c("IDX", paste("X", seq(1, nCol), sep=""));
    names(dataY) <- c("IDX", paste("Y", seq(1, nCol), sep=""));
  } else {
    names(dataX) <- c("IDX", "X");
    names(dataY) <- c("IDX", "Y");
  }

  polys <- polys[!duplicated(idx), ];

  # build PolyData result
  result <- data.frame(IDX=idx[!duplicated(idx)]);
  result$PID <- polys$PID;
  result$SID <- polys$SID;
  result <- merge(result, dataX, by="IDX");
  result <- merge(result, dataY, by="IDX");
  result <- result[, -1];

  if (!is.null(result)) {
    # copy the projection and zone attribute from input
    attr(result, "projection") <- attr(polys, "projection");
    attr(result, "zone") <- attr(polys, "zone");

    if (!is.PolyData(result, fullValidation = FALSE))
      class(result) <- c("PolyData", class(result));
  }

  return (result);
}

#==============================================================================
calcVoronoi <- function(xydata, xlim = NULL, ylim = NULL, eps = 1e-09,
                        frac = 0.0001) {
  .checkRDeps ("calcVoronoi", c("deldir"))

  xydata <- .validateXYData(xydata)
  if (is.character(xydata))
    stop(paste("Invalid X/Y data 'xydata'.\n", xydata, sep=""))
  if (nrow(xydata) < 2)
    stop("This function requires two or more input points.\n")

  if (!is.null(xlim) && (length(xlim) != 2 || diff(xlim) < 0))
    stop("Invalid 'xlim' argument.\n")
  if (!is.null(ylim) && (length(ylim) != 2 || diff(ylim) < 0))
    stop("Invalid 'xlim' argument.\n")

  # if limits null, calculate them by increasing each extent by 10% of its
  # range
  if (is.null(xlim))
    xlim <- range(xydata$X) + (rep(diff(range(xydata$X)) * 0.10, 2) * c(-1, 1))
  if (is.null(ylim))
    ylim <- range(xydata$Y) + (rep(diff(range(xydata$Y)) * 0.10, 2) * c(-1, 1))

  if (!is.null(attr(xydata, "projection"))
      && attr(xydata, "projection") == "LL")
    warning(paste(
"When the projection is \"LL\", this function naively treats the data as",
"a one-to-one projection. Consider converting it to UTM before running",
"this function.", sep="\n"));

  dd <- deldir::deldir(xydata$X, xydata$Y, rw=c(xlim, ylim), eps=eps, frac=frac,
               digits=7)
  if (is.null(dd))
    stop(paste(
"The call to \"deldir\", this function's dependency, failed.  The function",
"does not support linear (horizontal/vertical) data and that could be the",
"cause.", sep="\n"));

  # dd$dirsgs:
  #  (x1, y1) -> (x2, y2): a line segment in a Dirichlet tile
  #  ind1 and ind2: indices of the two points separated by this line segment
  dirsgs <- dd$dirsgs

  # dd$summary:
  #  (x, y): points used in the tesselation; may be less than input if some
  #          points judged to be the same (see "frac" argument)
  #  n.tside: number of side of the Dirichlet tile surrounding the point
  summary <- dd$summary

  # we will use "summary" below AND in ".addCorners" to determine closest
  # polygon for adding corners

  # each line segment belongs to two tiles (ind1 and ind2); create a data
  # structure with two copies of each point (i.e., a copy for each tile)
  pt1tile1 <- dirsgs[, c("x1", "y1", "ind1", "bp1")]
  pt1tile2 <- dirsgs[, c("x1", "y1", "ind2", "bp1")]
  pt2tile1 <- dirsgs[, c("x2", "y2", "ind1", "bp2")]
  pt2tile2 <- dirsgs[, c("x2", "y2", "ind2", "bp2")]
  names(pt1tile1) <- names(pt1tile2) <- names(pt2tile1) <-
    names(pt2tile2) <- c("X", "Y", "tile", "bp")
  res <- rbind(pt1tile1, pt1tile2, pt2tile1, pt2tile2)
  res <- res[order(res$tile), ]

  # the points in "res" are ordered by "tile" -- the same order of the
  # points in "summary"

  # add the "origin" point (twice) for each tile to the data structure
  res$Xo <- rep(summary[, "x"], times=2*summary[, "n.tside"])
  res$Yo <- rep(summary[, "y"], times=2*summary[, "n.tside"])
  # compute offset from each "origin"
  res$Xdiff <- res$X - res$Xo
  res$Ydiff <- res$Y - res$Yo
  # compute the arctan for each offset, so that we know the order
  # when building polygons
  res$atan2 <- atan2(res$Ydiff, res$Xdiff)

  # NOTE: the order set below is paramount to the next section of code.
  # DO NOT change the ordering unless you revise the next section as well.

  # order them and then eliminate duplicates
  res <- res[order(res$tile, res$atan2), ]
  res <- res[!duplicated(paste(res$tile, res$atan2)), ]

  # number of vertices in each polygon
  nVerts <- diff(c(which(!duplicated(res$tile)), length(res$tile)+1))

  #----------------------------------------------------------------------------
  # At this point, the order of points in each tile is sometimes inappropriate
  # for plotting. To obtain the best results, tiles with boundary points should
  # have one boundary point appear first in the tile and one last.
  #
  # Here is the idea:
  # We have two types of tiles:
  #   (1) "correct" tiles (completely interior OR boundary points in the first
  #       and last position) and
  #   (2) "incorrect" tiles (boundary points adjacent to each other in a given
  #       tile).
  # We only need to reorder vertices of the "incorrect" tiles. We can think of
  # the reordering as "rotating" the vertices until the first and last are the
  # boundary points.
  #
  # If F is a non-boundary point and T is a boundary point, a tile with seven
  # vertices may initially look like FFFTTFF.  If we rotate these vertices right
  # three times
  #   FFFFTTF (1st rotation)
  #   FFFFFTT (2nd rotation)
  #   TFFFFFT (3rd rotation)
  # we achieve the desired ordering.
  #
  # To calculate the rotation, we look at the index of the first boundary point
  # in the tile.  If we subtract its index from the number of vertices and add
  # one, we obtain the number of right rotations (e.g., 7 - 5 + 1 = 3 in the
  # above example).
  #
  # Given the number of right rotations, we can build a new index vector that
  # will "rotate" points within the desired tiles.

  # divide 'res' into two data frames, one where tiles touch two or fewer
  # boundary points and one where tiles touch more than two boundary points
  bpByTile <- split(res$bp, res$tile)
  noRotate <- rep((lapply(bpByTile, "sum") > 2),
                  times=lapply(bpByTile, "length"))
  resNoRotate <- res[noRotate, ]        # touch more than two
  res <- res[!noRotate, ]               # touch two or fewer

  # vertices in each tile
  tileLen <- diff(c(which(!duplicated(res$tile)), nrow(res)+1))
  # find the tiles to rotate
  bpByTile <- split(res$bp, res$tile)
  bpByTileWhich <- lapply(bpByTile, "which")
  toRotate <- (lapply(bpByTileWhich, "diff") == 1)
  toRotate[is.na(toRotate)] <- FALSE

  # focus on the tiles to fix
  toFix <- bpByTileWhich[toRotate]
  toFixLen <- tileLen[toRotate]
  toFixMin <- lapply(toFix, "min")

  # setup a "rotation" vector
  rotR <- unlist(toFixLen) - unlist(toFixMin)
  rotR <- -1 - rotR

  # a "-1" rotation means no change, "-2" means one to the right, ...
  rot <- rep(-1, length(tileLen))
  rot[toRotate] <- rotR

  # start with an index of 1 .. len for each tile
  newIDX <- unlist(lapply(tileLen, "seq"))

  # rotate those indices
  rot <- rep(rot, times=unlist(tileLen))
  newIDX <- newIDX + rot
  mod <- rep(unlist(tileLen), times=unlist(tileLen))
  newIDX <- (newIDX %% mod) + 1

  # adjust indices from 1 .. len to reflect their position in the
  # data frame
  base <- rep(which(!duplicated(res$tile)) - 1, times=unlist(tileLen))
  newIDX <- newIDX + base

  # reorder
  res <- res[newIDX, ]

  # combine the two 'res' data frames (the no rotate and rotate one)
  if (any(noRotate)) {
    res <- rbind(resNoRotate, res)
    # it is safe to 'order' by tile: according to the R help on 'order',
    # "any unresolved ties will be left in their original ordering."
    res <- res[order(res$tile), ]
  }
  #----------------------------------------------------------------------------

  # build the POS column
  maxPos <- max(nVerts);
  POS <- rep(seq(from=1, to=maxPos), times=length(nVerts));
  n <- vector();
  n[seq(1, by=2, length.out=length(nVerts))] <- nVerts;
  n[seq(2, by=2, length.out=length(nVerts))] <- maxPos - nVerts;
  b <- rep(c(TRUE, FALSE), length.out=length(n));
  POS <- POS[rep(b, n)];

  names(res)[is.element(names(res), "tile")] <- "PID"
  res$POS <- POS

  res <- res[, c("PID", "POS", "X", "Y")]
  class(res) <- c("PolySet", class(res))

  # polygons that straddle corners will be missing the corner points; add
  # the missing corner points now
  res <- .addCorners(res, summary)

  # ensure edges reach the proper extents
  res <- .expandEdges(res, data.frame(X=summary[, "x"], Y=summary[, "y"]),
                      xlim, ylim)

  # set attributes appropriate to the result
  attr(res, "projection") <- 1;
  attr(res, "zone") <- NULL;

  return (res);
}

#==============================================================================
clipLines <- function(polys, xlim, ylim, keepExtra = FALSE)
{
  polys <- .validatePolySet(polys);
  if (is.character(polys))
    stop(paste("Invalid PolySet 'polys'.\n", polys, sep=""));

  ret <- .clip(polys, xlim, ylim, isPolygons = FALSE, keepExtra);

  # .clip retains extra attributes

  if (!is.null(ret) &&
      !is.PolySet(ret, fullValidation = FALSE))
    class(ret) <- c("PolySet", class(ret));

  return (ret);
}

#==============================================================================
clipPolys <- function(polys, xlim, ylim, keepExtra = FALSE)
{
  polys <- .validatePolySet(polys);
  if (is.character(polys))
    stop(paste("Invalid PolySet 'polys'.\n", polys, sep=""));

  ret <- .clip(polys, xlim, ylim, isPolygons = TRUE, keepExtra);

  # .clip retains extra attributes

  if (!is.null(ret) &&
      !is.PolySet(ret, fullValidation = FALSE))
    class(ret) <- c("PolySet", class(ret));

  return (ret);
}

#==============================================================================
closePolys <- function(polys)
{
  polys <- .validatePolySet(polys);
  if (is.character(polys))
    stop(paste("Invalid PolySet 'polys'.\n", polys, sep=""));

  # save attributes, so they can be added to the resulting data frame
  attrNames <- setdiff(names(attributes(polys)),
                       c("names", "row.names", "class"));
  attrValues <- attributes(polys)[attrNames];

  xlim <- range(polys$X);
  ylim <- range(polys$Y);

  inRows <- nrow(polys);
  outCapacity <- as.integer(2 * inRows);

  # create the data structure that the C function expects
  if (!is.element("SID", names(polys))) {
    inID <- c(polys$PID, integer(length = inRows), polys$POS);
  } else {
    inID <- c(polys$PID, polys$SID, polys$POS);
  }
  inXY <- c(polys$X, polys$Y);
  limits <- c(xlim, ylim);

  # call the C function
  results <- .C("closePolys",
                inID = as.integer(inID),
                inXY = as.double(inXY),
                inVerts = as.integer(inRows),
                limits = as.double(limits),
                outID = integer(3 * outCapacity),
                outXY = double(2 * outCapacity),
                outRows = as.integer(outCapacity),
                outStatus = integer(1),
                PACKAGE = "PBSmapping");
  # note: outRows is set to how much space is allocated -- the C function
  #       should consider this

  if (results$outStatus == 1) {
    stop(
"Insufficient physical memory for processing.\n");
  }
  if (results$outStatus == 2) {
    stop(paste(
"Insufficient memory allocated for output.  Please upgrade to the latest",
"version of the software, and if that does not fix this problem, please",
"file a bug report.\n",
               sep = "\n"));
  }

  # determine the number of rows in the result
  outRows <- as.vector(results$outRows);

  # extract the data from the C function results
  if (outRows > 0) {
    d <- data.frame(PID = results$outID[1:outRows],
                    SID = results$outID[(outCapacity+1):(outCapacity+outRows)],
                    POS = results$outID[(2*outCapacity+1):(2*outCapacity+outRows)],
                    X = results$outXY[1:outRows],
                    Y = results$outXY[(outCapacity+1):(outCapacity+outRows)]);

    if (!is.element("SID", names(polys)))
      d$SID <- NULL;

    # restore the attributes, including 'projection' and 'zone'
    attributes(d) <- c(attributes(d), attrValues);

    if (!is.PolySet(d, fullValidation = FALSE))
      class(d) <- c("PolySet", class(d));

    return(d);
  } else {
    return(NULL);
  }
}

#combineEvents====================================2012-03-01
#  Combine measurements of events.
#---------------------------------------------------------NB
combineEvents <- function(events, locs, FUN, ..., bdryOK = TRUE)
{
  events <- .validateEventData(events);
  if (is.character(events))
    stop(paste("Invalid EventData 'events'.\n", events, sep=""));
  if (!is.element("Z", names(events))) {
    stop (
"EventData is missing required column 'Z'.\n");
  }

  # filter the boundary points
  if (!bdryOK)
    locs <- locs[locs$Bdry == 0, ];

  # make the list...
  if (is.element("SID", names(locs))) {
    colIDs <- 2;
    colNames <- c("PID", "SID", "Z");
    locs <- split(locs$EID, paste(locs$PID, locs$SID, sep = "-"));
  } else {
    colIDs <- 1;
    colNames <- c("PID", "Z");
    locs <- split(locs$EID, locs$PID);
  }

  summary <- lapply(locs,
                    function(x, FUN, events, ...)
                    {
                      FUN(events[is.element(events$EID, x), "Z"], ...);
                    },
                    FUN, events, ...);

  # The as.numeric function calls (below) are very important.  They prevent
  # elements in the PID and SID columns from becoming factors.
  output <- as.numeric(unlist(strsplit(names(summary), "-")));
  output <- as.data.frame(matrix(output, byrow = TRUE, ncol = colIDs));

  output <- cbind(output, unlist(summary));
  names(output) <- colNames;

  if (!is.null(output) &&
      !is.PolyData(output, fullValidation = FALSE))
    class(output) <- c("PolyData", class(output));

  return(output);
}

#==============================================================================
combinePolys <- function(polys)
{
  polys <- .validatePolySet(polys)
  if (is.character(polys))
    stop(paste("Invalid PolySet 'polys'.\n", polys, sep=""))

  # save the attributes for the data frame (.validatePolySet returns a data
  # frame)
  attrNames <- setdiff(names(attributes(polys)),
                       c("names", "row.names", "class"));
  attrValues <- attributes(polys)[attrNames];

  IDs <- .createIDs(polys, cols = c("PID", "SID"))

  # create/insert new PID/SID columns
  nVerts <- diff(c(which(!duplicated(IDs)), length(IDs) + 1))
  nPolys <- length(nVerts)
  polys$SID <- rep(1:nPolys, times=nVerts)
  polys$PID <- 1

  # reorder columns
  columns <- c("PID", "SID", "POS", "X", "Y")
  polys <- polys[, c(columns, setdiff(names(polys), columns))]

  # restore the attributes
  attributes(polys) <- c(attributes(polys), attrValues);

  return (polys)
}

#==============================================================================
convCP <- function(data, projection = NULL, zone = NULL)
{
  # do a preliminary test on the 'data' to see if it resulted from a call to
  # contourLines()
  if ((length(data) < 1) ||
      (length(data[[1]]) != 3) ||
      (any(!is.element(names(data[[1]]), c("level", "x", "y"))))) {
    stop(paste(
"convCP() accepts data from R's contourLines() function.  The format of",
"'data' does not match the expected format.\n",
               sep = "\n"));
  }

  # 'vData': 'data' converted to a vector
  vData <- unlist(data);

  levelsIdx <- which(names(vData) == "level");
  levelsIdx <- c(levelsIdx, length(vData) + 1);
  # 'nVerts': # of vertices for each polyline
  nVerts <- (diff(levelsIdx) - 1) / 2;
  # 'levelsIdx': indices to each 'level' element in vector 'vData'
  levelsIdx <- levelsIdx[-length(levelsIdx)];

  # create the X column
  # 'n': number of times to repeat each Boolean
  # 'b': vector of Booleans
  n <- vector();
  n[seq(1, by=2, length=length(nVerts))] <- nVerts;
  n[seq(2, by=2, length=length(nVerts))] <- nVerts + 1;
  n <- c(1, n[-length(n)], n[length(n)-1]);
  b <- rep(c(FALSE, TRUE), length.out=length(n));
  x <- vData [rep(b, n)];

  # create the Y column
  n <- vector();
  n[seq(1, by=2, length=length(nVerts))] <- nVerts + 1;
  n[seq(2, by=2, length=length(nVerts))] <- nVerts;
  b <- rep(c(FALSE, TRUE), length.out=length(n));
  y <- vData [rep(b, n)];

  # 'dupLevel': Boolean vector; TRUE when a level is a duplicate
  dupLevel <- duplicated(vData[levelsIdx]);

  # create the PID column
  pid <- vector()
  pid[which(!dupLevel)] <- 1:length(unique(vData[levelsIdx]));
  pid[which(dupLevel)] <-
    rep(pid[!dupLevel], diff(c(which(!dupLevel), length(levelsIdx)+1))-1);
  PolyData <- data.frame(PID=pid);
  pid <- rep(PolyData$PID, nVerts);

  # create the SID column
  d <- diff(c(which(!dupLevel), length(levelsIdx)+1));
  maxSid <- max(d);
  sid <- rep(seq(from=1, to=maxSid), times=length(d));
  n <- vector();
  n[seq(1, by=2, length.out=length(d))] <- d;
  n[seq(2, by=2, length.out=length(d))] <- maxSid - d;
  b <- rep(c(TRUE, FALSE), length.out=length(n));
  PolyData$SID <- sid[rep(b, n)];
  sid <- rep(PolyData$SID, nVerts);

  # create the POS column
  maxPos <- max(nVerts);
  pos <- rep(seq(from=1, to=maxPos), times=length(nVerts));
  n <- vector();
  n[seq(1, by=2, length.out=length(nVerts))] <- nVerts;
  n[seq(2, by=2, length.out=length(nVerts))] <- maxPos - nVerts;
  b <- rep(c(TRUE, FALSE), length.out=length(n));
  pos <- pos[rep(b, n)];

  # add level column to 'PolyData'
  PolyData$level <- vData[levelsIdx];

  # put the data frames together in a list
  PolySet <- data.frame(PID=pid, SID=sid, POS=pos, X=x, Y=y);

  if (!is.null(PolySet)) {
    # add attributes from arguments; we cannot add them from 'data' since
    # it isn't a PBS Mapping type
    if (!is.null(projection))
      attr(PolySet, "projection") <- projection;
    if (!is.null(zone))
      attr(PolySet, "zone") <- zone;

    # add class attributes as necessary
    if (!is.PolySet(PolySet, fullValidation = FALSE))
      class(PolySet) <- c("PolySet", class(PolySet));
  }
  if (!is.null(PolyData) &&
      !is.PolyData(PolyData, fullValidation = FALSE))
    class(PolyData) <- c("PolyData", class(PolyData));

  return (list(PolySet=PolySet, PolyData=PolyData));
}

#==============================================================================
convDP <- function(data, xColumns, yColumns)
{
  # validate all the arguments
  if (!any(is.element(names(data), c("PID", "EID")))) {
    stop(
"'data' must be either PolyData or EventData.\n");
  }
  # if 'data' contains both 'PID' and 'EID' columns, assume it is PolyData
  # and ignore 'EID'
  if (is.element("EID", names(data)) &&
      !is.element("PID", names(data))) {
    # replace EID with PID so we can .validatePolyData(), which is
    # less restrictive than .validateEventData()
    names(data)[which(is.element(names(data), "EID"))] <- "PID";
  }
  data <- .validatePolyData(data);
  if (is.character(data))
    stop(paste("Invalid PolyData 'data'.\n", data, sep=""));

  if (missing(xColumns) || missing(yColumns))
    stop(
"Must specify 'xColumns' and 'yColumns' vectors.\n");
  if ((length(xColumns) != length(yColumns))
      || (length(xColumns) == 0))
    stop(
"Length of 'xColumns' and 'yColumns' must be the same and greater than 0.\n");
  if (!all(is.element(xColumns, names(data)))
      || !all(is.element(yColumns, names(data))))
    stop(paste(
"One or more column specified in 'xColumns' or 'yColumns' does not",
"exist in 'data'.\n",
               sep = "\n"));

  result <- matrix(nrow=(nrow(data) * length(xColumns)), ncol=2);
  for (i in 1:length(xColumns)) {
    result[seq(i, by=length(xColumns), length=nrow(data)), 1] <-
      data[, xColumns[i]];
    result[seq(i, by=length(yColumns), length=nrow(data)), 2] <-
      data[, yColumns[i]];
  }
  result <- cbind(rep(1:length(xColumns), nrow(data)), result);
  if (is.element("SID", names(data)))
    result <- cbind(rep(data$SID, each=length(xColumns)), result);
  result <- cbind(rep(data$PID, each=length(xColumns)), result);

  result <- as.data.frame(result);
  if (is.element("SID", names(data))) {
    names(result) <- c("PID", "SID", "POS", "X", "Y");
  } else {
    names(result) <- c("PID", "POS", "X", "Y");
  }

  if (!is.null(result)) {
    # copy the 'projection' and 'zone' attribute from input
    attr(result, "projection") <- attr(data, "projection");
    attr(result, "zone") <- attr(data, "zone");

    if (!is.PolySet(result, fullValidation = FALSE))
      class(result) <- c("PolySet", class(result));
  }

  invisible(result);
}

#==============================================================================
convLP <- function(polyA, polyB, reverse = TRUE)
{
  polyA <- .validatePolySet(polyA);
  if (is.character(polyA))
    stop(paste("Invalid PolySet 'polyA'.\n", polyA, sep=""));
  polyB <- .validatePolySet(polyB);
  if (is.character(polyB))
    stop(paste("Invalid PolySet 'polyB'.\n", polyB, sep=""));

  # test for obvious user errors
  if ((length(unique(polyA$PID)) > 1 ||
       length(unique(polyB$PID)) > 1) ||
      (is.element("SID", names(polyA)) &&
       length(unique(polyA$SID)) > 1) ||
      (is.element("SID", names(polyB)) &&
       length(unique(polyB$SID)) > 1)) {
    warning(
"Each PolySet should contain only one polyline.\n");
  }

  result <- polyA[, c("X", "Y")];
  if (reverse) {
    result <- rbind(result, polyB[nrow(polyB):1, c("X", "Y")]);
  } else {
    result <- rbind(result, polyB[, c("X", "Y")]);
  }
  result <- cbind(1:nrow(result), result);
  result <- cbind(rep(polyA[1, "PID"], nrow(result)), result);
  names(result) <- c("PID", "POS", "X", "Y");

  if (!is.null(result)) {
    # add 'projection' and 'zone' attributes
    if (is.null(attr(polyA, "projection"))
        && !is.null(attr(polyB, "projection"))) {
      attr(result, "projection") <- attr(polyB, "projection");
    } else if ((!is.null(attr(polyA, "projection"))
                && is.null(attr(polyB, "projection")))
               || ((!is.null(attr(polyA, "projection")) &&
                    !is.null(attr(polyB, "projection")))
                   && (attr(polyA, "projection") ==
                       attr(polyB, "projection")))) {
      attr(result, "projection") <- attr(polyA, "projection");
    }

    if (is.null(attr(polyA, "zone"))
        && !is.null(attr(polyB, "zone"))) {
      attr(result, "zone") <- attr(polyB, "zone");
    } else if ((!is.null(attr(polyA, "zone"))
                && is.null(attr(polyB, "zone")))
               || ((!is.null(attr(polyA, "zone")) &&
                    !is.null(attr(polyB, "zone")))
                   && (attr(polyA, "zone") ==
                       attr(polyB, "zone")))) {
      attr(result, "zone") <- attr(polyA, "zone");
    }

    if (!is.PolySet(result, fullValidation = FALSE))
      class(result) <- c("PolySet", class(result));
  }

  invisible (result);
}

#==============================================================================
convUL <- function(xydata, km=TRUE, southern=NULL)
{
  xydata <- .validateXYData(xydata);
  if (is.character(xydata))
    stop(paste("Invalid X/Y data 'xydata'.\n", xydata, sep=""));
  # check for valid projection/zone attributes
  if (!is.element("projection", names(attributes(xydata)))
      || (!is.element(attr(xydata, "projection"), c("LL", "UTM")))) {
    stop(
"Missing or invalid projection attribute.\n");
  }

  # automagically set the zone attribute, if possible
  # SIMILAR code to calcArea()
  if ((attr(xydata, "projection") == "LL")
      && (is.null(attr(xydata, "zone")) || is.na(attr(xydata, "zone")))) {

    m <- mean(xydata$X);
    if ((m <= -180) || (m > 180)) {
      stop(paste(
"Attempted to automatically calculate the missing 'zone' attribute, but",
"that failed because the mean longitude falls outside the range",
"-180 < x <= 180.  Please manually set the attribute and then try again.\n",
                 sep = "\n"));
    }
    attr(xydata, "zone") <- ceiling((m + 180) / 6);
    message("convUL: For the UTM conversion, automatically detected zone ",
            attr(xydata, "zone"), ".");
  }
  else if (!is.element("zone", names(attributes(xydata)))
           || (attr(xydata, "zone") < 1) || (attr(xydata, "zone") > 60)) {
    stop(
"Invalid or missing zone attribute; possibly out of valid range.\n");
  }

  # determine whether in the northern or southern hemisphere
  if (is.null(southern)) {
    if (attr(xydata, "projection") == "LL")
      southern <- (mean(range(xydata$Y)) < 0)
    else if (attr(xydata, "projection") == "UTM")
      southern <- FALSE
    else
      stop("Projection attribute must be UTM or LL.")

    message("convUL: Converting coordinates within the ",
            ifelse(southern, "southern", "northern"), " hemisphere.")
  }

  inXY <- c(xydata$X, xydata$Y);
  outCapacity <- inVerts <- nrow(xydata);
  inProj <- attr(xydata, "projection");
  inZone <- attr(xydata, "zone");
  if (inProj == "UTM") {
    # C function expects data to be in metres...
    if (km) {
    	inXY <- inXY * 1000;
    }
    toUTM <- FALSE;
  } else {
    toUTM <- TRUE;
  }

  # call the C function
  results <- .C("convUL",
                inXY = as.double(inXY),
                inVerts = as.integer(inVerts),
                toUTM = as.integer(toUTM),
                zone = as.integer(inZone),
                southern = as.integer(southern),
                outXY = double(2 * outCapacity),
                outVerts = as.integer(outCapacity),
                outStatus = integer(1),
                PACKAGE = "PBSmapping");
  # note: outVerts is set to how much space is allocated -- the C function
  #       should consider this

  if (results$outStatus == 1) {
    stop(
"Insufficient physical memory for processing.\n");
  }
  if (results$outStatus == 2) {
    stop(paste(
"Insufficient memory allocated for output.  Please upgrade to the latest",
"version of the software, and if that does not fix this problem, please",
"file a bug report.\n",
               sep = "\n"));
  }

  # determine the number of rows in the result
  outRows <- as.vector(results$outVerts);

  # extract the data from the C function results
  if (outRows > 0) {
    # C function returns data in metres...
    if (inProj == "LL" && km)
      results$outXY <- results$outXY / 1000;

    xydata$X <- results$outXY[1:outRows];
    xydata$Y <- results$outXY[(outCapacity+1):(outCapacity+outRows)];

    if (inProj == "UTM") {
      attr(xydata, "projection") <- "LL";
    } else {
      attr(xydata, "projection") <- "UTM";
    }

    # no need to set zone attribute, same as input

    # no need to set class, same as input

    return(xydata);
  } else {
    return(NULL);
  }
}

#==============================================================================
dividePolys <- function(polys)
{
  polys <- .validatePolySet(polys)
  if (is.character(polys))
    stop(paste("Invalid PolySet 'polys'.\n", polys, sep=""))

  # save the attributes for the data frame (.validatePolySet returns a data
  # frame)
  attrNames <- setdiff(names(attributes(polys)),
                       c("names", "row.names", "class"));
  attrValues <- attributes(polys)[attrNames];

  IDs <- .createIDs(polys, cols = c("PID", "SID"))
  startsIdx <- which(!duplicated(IDs))

  POSstart <- polys[startsIdx, "POS"]
  # the calculation for POSend works even with one polygon
  POSend <- polys[c(startsIdx[-1] - 1, nrow(polys)), "POS"]
  isHole <- POSstart > POSend
  nVerts <- diff(c(startsIdx, length(IDs) + 1))

  # create new PID/SID columns
  newPID <- rep(1:sum(!isHole),
                times=diff(c(which(!(isHole)), length(isHole)+1)))
  newSID <- (1:length(newPID))
  subSID <- rep(newSID[!duplicated(newPID)] - 1,
                times=diff(c(which(!duplicated(newPID)), length(newPID)+1)))
  newSID <- newSID - subSID

  # insert new PID/SID columns
  polys$PID <- rep(newPID, times=nVerts)
  polys$SID <- rep(newSID, times=nVerts)

  # reorder columns
  columns <- c("PID", "SID", "POS", "X", "Y")
  polys <- polys[, c(columns, setdiff(names(polys), columns))]

  # restore the attributes
  attributes(polys) <- c(attributes(polys), attrValues);

  return (polys)
}

#==============================================================================
extractPolyData <- function(polys)
{
  polys <- .validatePolySet(polys);
  if (is.character(polys))
    stop(paste("Invalid PolySet 'polys'.\n", polys, sep=""));

  # filter out POS, oldPOS, X, Y columns
  polys <- polys[, !is.element(names(polys), c("POS", "oldPOS", "X", "Y"))];
  # check to see if it was reduced to a vector... (or if only PID/SID remained)
  if (length(setdiff(names(polys), c("PID", "SID"))) == 0) {
    stop(
"No extra fields to include in a PolyData object.\n");
  }

  # before building an index, ensure one doesn't already exist
  tempIDX <- date();
  if (is.element(tempIDX, names(polys))) {
    stop(paste(
"Failed attempt to create temporary index column; column already",
"exists.\n",
               sep = "\n"));
  }

  polys[[tempIDX]] <- .createIDs(polys, cols = c("PID", "SID"));
  pdata <- polys[!duplicated(polys[[tempIDX]]),
                 intersect(names(polys), c("PID", "SID", tempIDX))];

  # a function used to process each list element below
  processElement <- function(a) {
    a <- unique(as.vector(a));

    if (length(a) > 1) {
      warning("Added non-unique values to PolyData as 'NA'.\n");
      return (NA);
    } else {
      return (a);
    }
  }

  # process each column individually
  for (column in setdiff(names(polys), c(tempIDX, "PID", "SID"))) {
    temp <- split(polys[[column]], polys[[tempIDX]]);
    temp <- lapply(temp, processElement);
    temp <- data.frame(names(temp), unlist(temp));
    names(temp) <- c(tempIDX, column);
    pdata <- merge(pdata, temp, by = tempIDX, all.x = TRUE);
  }

  # remove index column
  pdata[[tempIDX]] <- NULL;

  # order the results
  if (is.element("SID", names(pdata))) {
    pdata <- pdata[order(pdata$PID, pdata$SID), ];
  } else {
    pdata <- pdata[order(pdata$PID), ];
  }

  # add a class, if necessary
  if (!is.null(pdata) &&
      !is.PolyData(pdata, fullValidation = FALSE))
    class(pdata) <- c("PolyData", class(pdata));

  return (pdata);
}

#findCells------------------------------2014-12-15
# Find events in grid cells
#--------------------------------------------NB/RH
findCells <- function (events, polys, includeBdry=NULL)
{
	events <- .validateEventData(events);
	if (is.character(events))
		stop(paste("Invalid EventData 'events'.\n", events, sep=""));
	polys <- .validatePolySet(polys);
	if (is.character(polys))
		stop(paste("Invalid PolySet 'polys'.\n", polys, sep=""));

	# get parameters to makeGrid: we'll use them to adjust the PID/SID
	# values in the output
	gridPars <- .getGridPars(polys, fullValidation = TRUE);
	if (is.null(gridPars))
		stop("Invalid 'polys'; appears altered since call to 'makeGrid'.\n");

	# prepare data for the C function
	nEvents <- nrow(events);
	brksX <- sort(unique(polys$X));
	brksY <- sort(unique(polys$Y));

	# call the C function for the X's
	resultsX <- .C("findCells",
		inPt = as.double(events$X),
		inPts = as.integer(nEvents),
		inBrk = as.double(brksX),
		inBrks = as.integer(length(brksX)),
		outCell = integer(nEvents),
		outBdry = integer(nEvents),
		outStatus = integer(1),
		PACKAGE = "PBSmapping");
	# call the C function for the Y's
	resultsY <- .C("findCells",
		inPt = as.double(events$Y),
		inPts = as.integer(nEvents),
		inBrk = as.double(brksY),
		inBrks = as.integer(length(brksY)),
		outCell = integer(nEvents),
		outBdry = integer(nEvents),
		outStatus = integer(1),
		PACKAGE = "PBSmapping");

	# build the output data structures
	d <- data.frame(EID = events$EID,
		PID = resultsX$outCell,
		PIDbdry = resultsX$outBdry,
		SID = resultsY$outCell,
		SIDbdry = resultsY$outBdry);

	bdryX <- d[d$PIDbdry == 1 & d$PID > 1 & d$PID < length(brksX), ];
	bdryX$PID <- bdryX$PID - 1;
	d <- rbind(d, bdryX);
	bdryY <- d[d$SIDbdry == 1 & d$SID > 1 & d$SID < length(brksY), ];
	bdryY$SID <- bdryY$SID - 1;
	d <- rbind(d, bdryY);

	d$Bdry <- as.integer(d$PIDbdry | d$SIDbdry);

	# filter out the bad stuff
	d <- d[d$PID != -1 & d$SID != -1, c("EID", "PID", "SID", "Bdry")];
	if (nrow(d)==0) return(NULL)
	d[d$PID == length(brksX), "PID"] <- length(brksX) - 1;
	d[d$SID == length(brksY), "SID"] <- length(brksY) - 1;

	# at this point:
	# - d uses PID/SID as if called with byrow = T, addSID = T
	# - transform to correct settings...
	if (gridPars$addSID && !gridPars$byrow) {
		# swap
		tmp <- d$PID;
		d$PID <- d$SID;
		d$SID <- tmp;
	} else if (!gridPars$addSID && gridPars$byrow) {
		d$PID <- (d$SID - 1) * (length(brksX) - 1) + d$PID;
		d$SID <- NULL;
	} else if (!gridPars$addSID && !gridPars$byrow) {
		d$PID <- (d$PID - 1) * (length(brksY) - 1) + d$SID;
		d$SID <- NULL;
	}

	# Check new argument `includeBdry' (2014-12-15)
	if (!is.null(includeBdry) && any(d$Bdry==1)) {
		if (includeBdry==0)
			d = d[!is.element(d$Bdry,1),]
		else {
			if (includeBdry==1)
				dd = d[order(d$PID,d$SID),]
			else
				dd = d[rev(order(d$PID,d$SID)),]
			dd = dd[!duplicated(dd$EID),]
			d  = dd[order(dd$EID),]
		}
		if (nrow(d)==0) return(NULL)
	}

	if (!is.LocationSet(d, fullValidation = FALSE))
		class(d) <- c("LocationSet", class(d));
	return (d);
}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~findCells


#findPolys------------------------------2014-12-15
# Find events in polygons
#--------------------------------------------NB/RH
findPolys <- function(events, polys, maxRows=1e+05, includeBdry=NULL)
{
	events <- .validateEventData(events);
	if (is.character(events))
		stop(paste("Invalid EventData 'events'.\n", events, sep=""));
	polys <- .validatePolySet(polys);
	if (is.character(polys))
		stop(paste("Invalid PolySet 'polys'.\n", polys, sep=""));

	# create the data structure that the C function expects
	inEvents <- nrow(events);
	inEventsID <- events$EID;
	inEventsXY <- c(events$X, events$Y);

	inPolys <- nrow(polys);
	if (!is.element("SID", names(polys))) {
		inPolysID <- c(polys$PID, integer(length = inPolys), polys$POS);
	} else {
		inPolysID <- c(polys$PID, polys$SID, polys$POS);
	}
	inPolysXY <- c(polys$X, polys$Y);

	# the maximum number of rows in the results occurs when each event
	# lies in all of the polygons; we can easily calculate this case
	# but it often results in too much allocated memory

	# instead, simply use an argument that the user can tune (this makes
	# the behavior more consistent with joinPolys, too)
	outCapacity <- maxRows;

	# call the C function
	results <- .C("findPolys",
		inEventsID = as.integer(inEventsID),
		inEventsXY = as.double(inEventsXY),
		inEvents = as.integer(inEvents),
		inPolysID = as.integer(inPolysID),
		inPolysXY = as.double(inPolysXY),
		inPolys = as.integer(inPolys),
		outID = integer(4 * outCapacity),
		outRows = as.integer(outCapacity),
		outStatus = integer(1),
		PACKAGE = "PBSmapping");
	# note: outRows is set to how much space is allocated -- the C function should consider this

	if (results$outStatus == 1) {
		stop(paste(
"Insufficient physical memory for processing.",
"Try reducing this function's 'maxRows' argument to reduce its memory",
"requirements.\n", sep = "\n"));
	}
	if (results$outStatus == 2) {
		stop(paste(
"Insufficient memory allocated for output.",
"Try increasing this function's 'maxRows' argument to increase the memory",
"that it allocates for the output.\n", sep = "\n"));
	}

	# determine the number of rows in the result
	outRows <- as.vector(results$outRows);
	if (outRows == 0) return(NULL)

	# extract the data from the C function results
	d <- data.frame(EID = results$outID[1:outRows],
		PID = results$outID[(outCapacity+1):(outCapacity+outRows)],
		SID = results$outID[(2*outCapacity+1):(2*outCapacity+outRows)],
		Bdry = results$outID[(3*outCapacity+1):(3*outCapacity+outRows)]);

	# Check new argument `includeBdry' (2014-12-15)
	if (!is.null(includeBdry) && any(d$Bdry==1)) {
		if (includeBdry==0)
			d = d[!is.element(d$Bdry,1),]
		else {
			d0 = d[is.element(d$Bdry,0),]
			d1 = d[is.element(d$Bdry,1),]
			if (includeBdry==1)
				dd = d1[order(d1$PID,d1$SID),]
			else
				dd = d1[rev(order(d1$PID,d1$SID)),]
			dd = dd[!duplicated(dd$EID),]
			dd = rbind(d0,dd)
			d  = dd[order(dd$EID,dd$PID,dd$SID),]
		}
		if (nrow(d)==0) return(NULL)
	}

	if (!is.element("SID", names(polys)))
		d$SID <- NULL;

	if (!is.LocationSet(d, fullValidation = FALSE))
		class(d) <- c("LocationSet", class(d));

	return(d);
}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~findPolys


fixBound <- function(polys, tol)
{
  polys <- .validatePolySet(polys);
  if (is.character(polys))
    stop(paste("Invalid PolySet 'polys'.\n", polys, sep=""));

  if (length(tol) == 1) {
    tolX <- tolY <- tol;
  } else if (length(tol) == 2) {
    tolX <- tol[1]; tolY <- tol[2];
  } else {
    stop(
"Invalid value for 'tol'.  It must be a vector of length 1 or 2.\n");
  }

  xlim <- range(polys$X);
  rangeX <- abs(diff(xlim));
  xDiff <- c(-(tolX * rangeX), (tolX * rangeX));
  leftBnd  <- c(xDiff + xlim[1]);
  rightBnd <- c(xDiff + xlim[2]);

  ylim <- range(polys$Y);
  rangeY <- abs(diff(ylim));
  yDiff <- c(-(tolY * rangeY), (tolY * rangeY));
  lowerBnd <- c(yDiff + ylim[1]);
  upperBnd <- c(yDiff + ylim[2]);

  polys$X[leftBnd[1]  < polys$X & polys$X < leftBnd[2]] <- xlim[1];
  polys$X[rightBnd[1] < polys$X & polys$X < rightBnd[2]] <- xlim[2];
  polys$Y[lowerBnd[1] < polys$Y & polys$Y < lowerBnd[2]] <- ylim[1];
  polys$Y[upperBnd[1] < polys$Y & polys$Y < upperBnd[2]] <- ylim[2];

  # retains attributes including  'projection' and 'zone'

  if (!is.null(polys) &&
      !is.PolySet(polys, fullValidation = FALSE))
    class(polys) <- c("PolySet", class(polys));

  return(polys);
}

#==============================================================================
fixPOS <- function(polys, exteriorCCW = NA)
{
  polys <- .validatePolySet(polys);
  if (is.character(polys))
    stop(paste("Invalid PolySet 'polys'.\n", polys, sep=""));

  if (is.na(exteriorCCW)) {
    exteriorCCW <- -1;
  } else {
    exteriorCCW <- as.logical(exteriorCCW);
    if (is.na(exteriorCCW)) {
      stop(
"Invalid value for 'exteriorCCW'.  Assign it a Boolean value or NA.\n");
    }
  }

  polys <- .rollupPolys(polys, rollupMode = 3, exteriorCCW = exteriorCCW,
                        closedPolys = -1, addRetrace = FALSE);

  # .rollupPolys() result retains attributes of 'polys'

  if (!is.null(polys) &&
      !is.PolySet(polys, fullValidation = FALSE))
    class(polys) <- c("PolySet", class(polys));

  return (polys);
}

#=============================================================================
importPolys <- function(PolySet, PolyData=NULL, projection=NULL, zone=NULL)
{
  x <- as.PolySet(read.table(PolySet, header=TRUE), projection=projection,
                  zone=zone)
  if (!is.null(PolyData)) {
    xDat <- as.PolyData(read.table(PolyData, header=TRUE))
    attr(x, "PolyData") <- xDat
  }
  return(x)
}

#=============================================================================
importEvents <- function(EventData, projection=NULL, zone=NULL)
{
  return(as.EventData(read.table(EventData, header=TRUE),
                      projection = projection, zone = zone))
}

#=============================================================================
importLocs <- function(LocationSet)
{
  return(as.LocationSet(read.table(LocationSet, header=TRUE)))
}
#=============================================================================

#importGSHHS----------------------------2015-02-12
# Import a GSHHG binary file provided by Paul Wessel.
# The C call on gshhs using x-limits (0,360) is the only
# limit that ties the polygons near the Greenwhich meridian together.
# The meshing at the International date line is not yet ideal.
#--------------------------------------------NB/RH
importGSHHS <- function(gshhsDB, xlim, ylim, maxLevel=4, n=0, useWest=FALSE)
{
	# Transform all X-coordinates to lie between 0 and 360
	normAngle = function(x){(x + 360.)%%360.}

	# Determine if the specified limits lie in the western hemisphere
	isitWest  = function(lim){
		world = 1:360
		if (diff(lim)>0) span = lim[1]:lim[2]
		else span = setdiff(world,lim[2]:lim[1])
		sum(span>180) > sum(span<=180) }

	# Get gshhsDB filename (with and without path information)
	if (missing(gshhsDB))
		gshhsDB <- paste(system.file(package="PBSmapping"), "gshhs_f.b", sep="/")
	else
		gshhsDB <- path.expand(gshhsDB)
	if (!file.exists(gshhsDB))
		stop("Unable to find gshhsDB \"", gshhsDB, "\"")
	gshhsDB.name = basename(substitute(gshhsDB))
	isBorder = grepl("borders",gshhsDB.name)

	#-----FORMAT-BORDERS-CHANGES-------------------
	# Sometimes border longitudes lie between (-180,180) or (0,360).
	# Easiest solution is to import all and use refocusWorld.
	#----------------------------------------------
	if (isBorder) {
		.checkClipLimits(c(xlim, ylim))
		#DBdir  = dirname(substitute(gshhsDB))
		#DBinfo = readLines(paste0(DBdir,"/README.TXT"))
		#DBver  = strsplit(sub("^[[:space:]]+","",DBinfo[grep("Version",DBinfo)[1]]),split=" ")[[1]][2]
		xres <- .Call("importGSHHS", as.character(gshhsDB), as.numeric(c(0,360,-90,90)),
			as.integer(maxLevel), as.integer(n), PACKAGE = "PBSmapping")
		xPS  = refocusWorld(as.PolySet(as.data.frame(xres),projection="LL"),xlim=xlim,ylim=ylim)
		if (is.null(xPS) || !length(xPS$PID)) return(NULL)
#browser()
		attr(xPS, "PolyData") <- attr(xres, "PolyData")
		if (useWest) xlim <- xlim - 360.
		if (!all(xlim==attributes(xPS)$rf.xlim)) {  # refocusWorld adjusts xlims
			xadj = (xlim-attributes(xPS)$rf.xlim)[1]
			xPS$X = xPS$X + xadj
		}
		return(xPS)
	}
	#----------------------------------------------
	# Headers in GSHHS database range from 0 to 360 while the underlying data
	# range from -180 to 180; if xlim spans the globe use c(0,360) and
	# refocusWorld the world
	#----------------------------------------------
	xlim.orig <- NULL
	if (abs(diff(xlim)) >= 360) {
		xlim.orig <- xlim
		xres <- .Call("importGSHHS", as.character(gshhsDB), as.numeric(c(0,360,ylim)),
			as.integer(maxLevel), as.integer(n), PACKAGE = "PBSmapping")
		xPS  = refocusWorld(as.PolySet(as.data.frame(xres),projection="LL"),xlim=xlim,ylim=ylim)
		if (is.null(xPS) || !length(xPS$PID)) return(NULL)
#browser();return()
		PolyData    <- attr(xres, "PolyData")
		clipAsPolys <- attr(xres, "clipAsPolys")
		if (useWest) xlim <- xlim - 360.
		if (!all(xlim==attributes(xPS)$rf.xlim)) {  # refocusWorld adjusts xlims
			xadj = (xlim-attributes(xPS)$rf.xlim)[1]
			xPS$X = xPS$X + xadj
		}
	}
	else {
		#----------------------------------------------
		# Try to subset the workd using the C code (tricky)
		# Set up limits: for lon, (-20, 360)
		#----------------------------------------------
		limits <- list(c(xlim, ylim))
		.checkClipLimits(limits[[1]])
		overlap = ifelse(grepl("gshhs",gshhsDB.name),-20.,0.)
		nxlim   = normAngle(xlim)

		# Check if xlim spans the Greenwich meridian
		#----------------------------------------------
		green0 = (nxlim[1] >= nxlim[2]) | nxlim[2] > 340.
		green1 = nxlim[1]==0  && grepl("gshhs",gshhsDB.name)
		green2 = nxlim[2]==0  && grepl("gshhs",gshhsDB.name)
		green  = green0 | green1 | green2
		if (green) {
			if (green1) {
				limits[[1]] = c(0.,normAngle(xlim[2]),ylim)
				#limits[[2]] = c(340.,360.,ylim) # maybe not necessary
				bighalf = "R"
			} else if (green2) {
				limits[[1]] = c(normAngle(xlim[1]),360.,ylim)
				limits[[2]] = c(-20.,0.,ylim)
				bighalf = "L"
			} else {
				if ( (360-nxlim[1]) >= nxlim[2] ) {
					limits[[1]] = c(nxlim[1],360.,ylim)
					limits[[2]] = c(overlap,nxlim[2],ylim)
					bighalf = "L"
				}
				else {
					limits[[1]] = c(overlap,xlim[2],ylim)
					limits[[2]] = c(nxlim[1],360.,ylim)
					bighalf = "R"
				}
			}
		}
		isWest = isitWest(limits[[1]][1:2]) # defined solely by user's xlim
#browser();return()

		# Initilialize
		x = xbit = list()
		pidlen = 0
		clipAsPolys = logical()
		PolyData = data.frame()

		#Go through 1 or 2 limits, depending on whether `xlim' includes 0 degrees
		#----------------------------------------------
		for (i in 1:length(limits)) {
			.checkClipLimits(limits[[i]])
			# return an R object -- C call unfortunately converts western hemisphere to negative
			#xres <- .Call("importGSHHS", as.character(gshhsDB), as.numeric(c(0,360,0,90)), as.integer(maxLevel), as.integer(n), PACKAGE = "PBSmapping") # debug only
			#browser();return()
			xres <- .Call("importGSHHS", as.character(gshhsDB), as.numeric(limits[[i]]),
				as.integer(maxLevel), as.integer(n), PACKAGE = "PBSmapping")
			if (is.null(xres) || !length(xres$PID)) next
#if(i==2){browser();return()}

			# Attempt to reverse the mess created by the C-code.
			# If xlim spans the Internation Date Line but not Greenwich
			#----------------------------------------------
			dateline = (normAngle(xlim[1]) >= 45. & normAngle(xlim[2]) <= 315.) & !green
			#if ( any(is.element(c(180,-180),floor(limits[[i]][1]):ceiling(limits[[i]][2]))) & !green ) {
			if (dateline) {
				xneg  = sapply(split(xres[["X"]],xres[["PID"]]),function(x){all(x<0)})
				pneg  = names(xneg[xneg])  # PIDs with all-negative coordinates
				xmove = is.element(xres[["PID"]],pneg)
				xres[["X"]][xmove] = normAngle(xres[["X"]][xmove])
			}
#browser();return()
			#isWestMess = median(xres[["X"]]) < 0  # `isWestMess' can be different than `isWest'

			#collect attributes also from each C call
			pdata <- as.data.frame(attr(xres, "PolyData"))
			clipAsPolys <- c(clipAsPolys, attr(xres, "clipAsPolys"))

			# Renumber the PIDs to accommodate a possible second set
			pidnam = unique(xres$PID)
			pid    = (pidlen+1):(pidlen+length(pidnam))
			names(pid) = pidnam
			xres[["oldPID"]] = xres[["PID"]]
			xres[["PID"]] = as.vector(pid[as.character(xres[["PID"]])])
			pdata$oldPID = pdata$PID
			pdata$PID = as.vector(pid[as.character(pdata$PID)])
			pidlen = rev(pid)[1]
			for (j in names(xres)){
				if (length(x)==0) x[[j]] = xres[[j]]
				else              x[[j]] = c(x[[j]],xres[[j]]) # add on results from second limit
			}
			PolyData = rbind(PolyData,pdata)
			xbit[[i]] = xres  # collect the separate results for debugging only
		}
#browser();return()
		#if ( "PBSmodelling" %in% rownames(installed.packages()) ) { # primarily for debugging
		#	if ( exists("XBIT",envir=.PBSmapEnv) )
		#		PBSmodelling::tget(XBIT, tenv=.PBSmapEnv) else XBIT = list()
		#	XBIT[[gshhsDB.name]] = xbit
		#	PBSmodelling::tput(XBIT, tenv=.PBSmapEnv)
		#	PBSmodelling::tput(isWest, tenv=.PBSmapEnv)
		#}
		if (length(x)==0 || !length(x$PID))
			return(NULL)

		# headers in the GSHHS database range from 0 to 360 while the underlying data
		# ranges from -180 to 180; if our xlim > 180, shift it
		xoff = 0.
		#else if (max(xlim) > 180) { # too harsh for places like NZ
		if (median(xlim) > 180) {
			xlim.orig <- xlim
			if (useWest) {
				xlim <- xlim - 360.
				xoff = 360.
				if (median(x[["X"]])>0) x[["X"]] = x[["X"]] - 360.
			} else {
				if (median(x[["X"]])<0) x[["X"]] = x[["X"]] + 360.  # adjust for .Call("importGSHHS")
			}
		}
		XLIM = range(x[["X"]]) # imported data from C-function

		# Convert the list to a PolySet
		xPS <- as.PolySet(as.data.frame(x), projection = "LL")
	} #----- end else
#browser();return()

	# Clip the PolySet
	#----------------------------------------------
	if (any(clipAsPolys))
		xPS <- clipPolys(xPS, xlim=xlim, ylim=ylim)
	else
		xPS <- clipLines(xPS, xlim=xlim, ylim=ylim)
	if (is.null(xPS)) return(NULL)
	xPS$oldPOS <- NULL

	# Sort so that holes immediately follow their parents
	xPS <- xPS[order(xPS$PID, xPS$SID), ]

	# Display a message if the output longitudes differ from xlim
	if (!is.null(xlim.orig))
		message("importGSHHS: input xlim was (", paste(xlim.orig, collapse=", "),
					") and the longitude range of the extracted data is (",
					paste(range(xPS$X), collapse=", "), ").")
	attr(xPS, "PolyData") <- PolyData
	return (xPS)
}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~importGSHHS

#==============================================================================
is.EventData <- function(x, fullValidation = TRUE)
{
  if (fullValidation) {
    msg <- .validateEventData(x)
    if (is.character(msg))
      return (FALSE)
  }

  return (inherits(x, "EventData", which = TRUE) == 1);
}

#==============================================================================
is.LocationSet <- function(x, fullValidation = TRUE)
{
  if (fullValidation) {
    msg <- .validateLocationSet(x)
    if (is.character(msg))
      return (FALSE)
  }

  return (inherits(x, "LocationSet", which = TRUE) == 1);
}

#==============================================================================
is.PolyData <- function(x, fullValidation = TRUE)
{
  if (fullValidation) {
    msg <- .validatePolyData(x)
    if (is.character(msg))
      return (FALSE)
  }

  return (inherits(x, "PolyData", which = TRUE) == 1);
}

#==============================================================================
is.PolySet <- function(x, fullValidation = TRUE)
{
  if (fullValidation) {
    msg <- .validatePolySet(x)
    if (is.character(msg))
      return (FALSE)
  }

  return (inherits(x, "PolySet", which = TRUE) == 1);
}

#==============================================================================
isConvex <- function(polys)
{
  polys <- .validatePolySet(polys);
  if (is.character(polys))
    stop(paste("Invalid PolySet 'polys'.\n", polys, sep=""));

  # create the data structure that the C function expects
  inPolys <- nrow(polys);
  if (!is.element("SID", names(polys))) {
    inPolysID <- c(polys$PID, vector(length = inPolys));
  } else {
    inPolysID <- c(polys$PID, polys$SID);
  }
  inPolysXY <- c(polys$X, polys$Y);

  # if each point appears in more than two polygons, this number will be
  # too low
  outCapacity <- inPolys;

  # call the C function
  results <- .C("isConvex",
                inID = as.integer(inPolysID),
                inXY = as.double(inPolysXY),
                inVerts = as.integer(inPolys),
                outID = integer(2 * outCapacity),
                outResult = integer(outCapacity),
                outRows = as.integer(outCapacity),
                outStatus = integer(1),
                PACKAGE = "PBSmapping");
  # note: outRows is set to how much space is allocated -- the C function
  #       should consider this

  if (results$outStatus == 1) {
    stop(
"Insufficient physical memory for processing.\n");
  }
  if (results$outStatus == 2) {
    stop(paste(
"Insufficient memory allocated for output.  Please upgrade to the latest",
"version of the software, and if that does not fix this problem, please",
"file a bug report.\n",
               sep = "\n"));
  }

  # determine the number of rows in the result
  outRows <- as.vector(results$outRows);

  # extract the data from the C function results
  if (outRows > 0) {
    d <- data.frame(PID = results$outID[1:outRows],
                    SID = results$outID[(outCapacity+1):(outCapacity+outRows)],
                    convex = results$outResult[1:outRows]);

    # when I had this conversion in the about data.frame() call, the 'convex'
    # column became factors...
    d$convex <- as.logical(d$convex);

    if (!is.element("SID", names(polys)))
      d$SID <- NULL;

    if (!is.PolyData(d, fullValidation = FALSE))
      class(d) <- c("PolyData", class(d));

    return(d);
  } else {
    return(NULL);
  }
}

#==============================================================================
isIntersecting <- function(polys, numericResult = FALSE)
{
  polys <- .validatePolySet(polys);
  if (is.character(polys))
    stop(paste("Invalid PolySet 'polys'.\n", polys, sep=""));

  # create the data structure that the C function expects
  inPolys <- nrow(polys);
  if (!is.element("SID", names(polys))) {
    inPolysID <- c(polys$PID, integer(length = inPolys), polys$POS);
  } else {
    inPolysID <- c(polys$PID, polys$SID, polys$POS);
  }
  inPolysXY <- c(polys$X, polys$Y);

  # if ecah point appears in more than two polygons, this number will be
  # too low
  outCapacity <- inPolys;

  # call the C function
  results <- .C("isIntersecting",
                inPolysID = as.integer(inPolysID),
                inPolysXY = as.double(inPolysXY),
                inPolys = as.integer(inPolys),
                numericResult = as.integer(numericResult),
                outID = integer(2 * outCapacity),
                outResult = integer(outCapacity),
                outRows = as.integer(outCapacity),
                outStatus = integer(1),
                PACKAGE = "PBSmapping");
  # note: outRows is set to how much space is allocated -- the C function
  #       should consider this

  if (results$outStatus == 1) {
    stop(
"Insufficient physical memory for processing.\n");
  }
  if (results$outStatus == 2) {
    stop(paste(
"Insufficient memory allocated for output.  Please upgrade to the latest",
"version of the software, and if that does not fix this problem, please",
"file a bug report.\n",
               sep = "\n"));
  }

  # determine the number of rows in the result
  outRows <- as.vector(results$outRows);

  # extract the data from the C function results
  if (outRows > 0) {
    d <- data.frame(PID = results$outID[1:outRows],
                    SID = results$outID[(outCapacity+1):(outCapacity+outRows)],
                    intersecting = results$outResult[1:outRows]);

    if (!numericResult)
      d$intersecting <- as.logical(d$intersecting);

    if (!is.element("SID", names(polys)))
      d$SID <- NULL;

    if (!is.PolyData(d, fullValidation = FALSE))
      class(d) <- c("PolyData", class(d));

    return(d);
  } else {
    return(NULL);
  }
}


#joinPolys------------------------------2014-03-10
# Join one or two PolySets using a logic operation.
#-----------------------------------------------NB
joinPolys <- function(polysA, polysB = NULL, operation = "INT")
{
  polysA <- .validatePolySet(polysA)
  if (is.character(polysA))
    stop(paste("Invalid PolySet 'polysA'.\n", polysA, sep=""))

  if (!is.null(polysB)) {
    polysB <- .validatePolySet(polysB)
    if (is.character(polysB))
      stop(paste("Invalid PolySet 'polysB'.\n", polysB, sep=""))
  }

  validOps <- c("INT", "UNION", "DIFF", "XOR")
  if (!is.element(operation, validOps)) {
    stop(paste(
"Invalid \"operation.\"  Must be one of: ", paste(validOps, collapse=", "),
"\n"))
  }
  op <- which(is.element(validOps, operation)) - 1

  inputHasSID <- is.element("SID", names(polysA))
  if (!is.element("SID", names(polysA))) {
    polysA$SID <- 1
  }
  if (!is.element("SID", names(polysB))) {
    polysB$SID <- 1
  }

  # call the C function
  results <- .Call("joinPolys",
                   operation = as.integer(op),
                   sPID = as.integer(polysA$PID),
                   sSID = as.integer(polysA$SID),
                   sPOS = as.integer(polysA$POS),
                   sX   = as.numeric(polysA$X),
                   sY   = as.numeric(polysA$Y),
                   cPID = as.integer(polysB$PID),
                   cSID = as.integer(polysB$SID),
                   cPOS = as.integer(polysB$POS),
                   cX   = as.numeric(polysB$X),
                   cY   = as.numeric(polysB$Y),
                   PACKAGE = "PBSmapping")

  if (is.null(results)){
    return(NULL)
  } else if (nrow(results) > 0) {
    # the C routine might add an SID column when one previously didn't exist...
    if (!inputHasSID && all(results$SID == 1)) {
      results$SID <- NULL
    }

    # copy the 'projection' and 'zone' attribute from the input
    attr(results, "projection") <- attr(polysA, "projection")
    attr(results, "zone") <- attr(polysA, "zone")

    if (!is.PolySet(results, fullValidation = FALSE))
      class(results) <- c("PolySet", class(results))

    return(results)
  } else {
    return(NULL)
  }
}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~joinPolys

locateEvents <- function(EID, n = 512, type = "p", ...)
{
  if (!missing(EID)) {
    if (any(!is.numeric(EID))) {
      stop(paste(
"When the 'EID' argument exists, it must contain a vector of numeric",
"values.\n",
                 sep = "\n"));
    }
    n <- length(EID);
  }

  events <- locator(n, type, ...);
  # when 0 pts., S-PLUS returns an object of class "missing", whereas
  # R returns NULL; when > 0 pts., both return an object of class "list"

  if (!missing(EID) || (data.class(events) == "list")) {
    if (missing(EID)) {
      EventData <- data.frame(EID = 1:length(events$x),
                              X = events$x,
                              Y = events$y);
    } else {
      EventData <- data.frame(EID = EID,
                              X = c(events$x, rep(NA, n - length(events$x))),
                              Y = c(events$y, rep(NA, n - length(events$y))));
    }

    if (!is.null(EventData)) {
      # set the projection attribute
      attr(EventData, "projection") <- options()$map.projection;
      # does not set "zone" attribute as per documentation

      if (!is.EventData(EventData, fullValidation = FALSE))
        class(EventData) <- c("EventData", class(EventData));
    }

    return (EventData);
  } else {
    return (NULL);
  }
}

#==============================================================================
locatePolys <- function(pdata, n = 512, type = "o", ...)
{
  if (!missing(pdata)) {
    pdata <- .validatePolyData(pdata);
    if (is.character(pdata))
      stop(paste("Invalid PolyData 'pdata'.\n", pdata, sep=""));
    polys <- nrow(pdata);
  } else {
    polys <- 1;
  }

  PID <- 1;
  SID <- NULL;
  output <- NULL;
  for (i in 1:polys) {
    # replace as necessary from 'pdata'
    if (!missing(pdata)) {
      PID <- pdata[i, "PID"];
      if (is.element("SID", names(pdata)))
        SID <- pdata[i, "SID"];
      if (is.element("n", names(pdata)))
        n <- pdata[i, "n"];
      if (is.element("type", names(pdata)))
        type <- pdata[i, "type"]
    }

    pts <- locator(n = abs(n), type = type, ...);
    if (length(pts$x) == 0) {
      stop(
"You must locate at least one point.\n");
    }
    if (type == "l" || type == "o")
      lines(pts$x[c(1, length(pts$x))], pts$y[c(1, length(pts$y))], ...);

    if (missing(pdata) || (!missing(pdata) && !is.element("n", names(pdata))))
      # keep the sign in case we are generating a hole with no parent
      n <- sign(n) * length(pts$x);

    # create vectors for data frame
    PID <- rep(PID, abs(n));
    if (!is.null(SID))
      SID <- rep(SID, abs(n));
    if (n > 0) POS <- 1:n else POS <- abs(n):1;
    if (!missing(pdata) && is.element("n", names(pdata))) {
      pts$x <- c(pts$x, rep(NA, abs(n) - length(pts$x)));
      pts$y <- c(pts$y, rep(NA, abs(n) - length(pts$y)));
    }

    if (!is.null(SID)) {
      output <- rbind(output, data.frame(PID = PID, SID = SID, POS = POS,
                                         X = pts$x, Y = pts$y));
    } else {
      output <- rbind(output, data.frame(PID = PID, POS = POS,
                                         X = pts$x, Y = pts$y));
    }
  }

  if (!is.null(output)) {
    # set the projection attribute
    attr(output, "projection") <- options()$map.projection;
    # does not set "zone" attribute as per documentation

    if (!is.PolySet(output, fullValidation = FALSE))
      class(output) <- c("PolySet", class(output));
  }

  return(output)
}

#==============================================================================
# The code to check the Xs and Ys for errors, and then generate the Xs and Ys,
# came from Rowan Haigh.
makeGrid <- function(x, y, byrow = TRUE, addSID = TRUE, projection = NULL,
                     zone = NULL)
{
  # check Xs and Ys for errors
  dX <- diff(x);
  dY <- diff(y);
  if (any(c(dX, dY) < 0)) {
    stop(
"Either the 'x' vector or the 'y' vector was not sequential.\n");
  }
  if (length(x) < 2 || length(y) < 2) {
    stop(
"Both x and y must have a length > 1.\n");
  }

  # generate the Xs and Ys
  numX <- length(x);
  numY <- length(y);
  x1 <- rep(x[1:(numX - 1)], times = numY - 1);
  y1 <- rep(y[1:(numY - 1)], each = numX - 1);
  x2 <- x1 + rep(dX, times = numY - 1);
  y2 <- y1 + rep(dY, each = numX - 1);

  xP <- as.vector(rbind(x1, x2, x2, x1))
  yP <- as.vector(rbind(y1, y1, y2, y2))

  # set up for (byrow && addSID)
  pid <- rep(rep(1:(numX - 1), each = 4), times = numY - 1);
  sid <- rep(1:(numY - 1), each = 4 * (numX - 1));
  pos <- rep(1:4, times = ((numX - 1) * (numY - 1)));

  # make data frame
  pSet <- data.frame(PID = pid, SID = sid, POS = pos, X = xP, Y = yP);
  if (!addSID || !byrow) {
    # adjust for alternative indexing
    pSet <- .createGridIDs(pSet, addSID, byrow);
  }

  if (!is.null(pSet)) {
    if (!is.null(projection))
      attr(pSet, "projection") <- projection;
    if (!is.null(zone))
      attr(pSet, "zone") <- zone;

    if (!is.PolySet(pSet, fullValidation = FALSE))
      class(pSet) <- c("PolySet", class(pSet));
  }

  return(pSet);
}

#==============================================================================
# the default value for 'propVals' is less than optimal if 'breaks' is a vector
# of length 1...
makeProps <- function(pdata, breaks, propName = "col",
                      propVals = 1:(length(breaks) - 1))
{
  pdata <- .validatePolyData(pdata);
  if (is.character(pdata))
    stop(paste("Invalid PolyData 'pdata'.\n", pdata, sep=""));

  if (!is.element("Z", names(pdata)))
    stop(
"'Z' column is missing in 'pdata'.\n");

  if ((length(breaks) == 1 && breaks < 2) || (length(breaks) == 0))
    stop(
"You must have at least two breaks.\n");

  if (is.null(propVals))
    stop(
"'propVals' cannot be NULL.\n");

  if (length(breaks) == 1) {
    numBreaks <- breaks;
  } else {
    numBreaks <- length(breaks) - 1;
  }

  propVals <- rep(propVals, length.out = numBreaks);
  propIdx <- cut(pdata$Z, breaks, include.lowest = TRUE);
  pdata[[propName]] <- propVals[unclass(propIdx)];

  if (!is.null(pdata) &&
      !is.PolyData(pdata, fullValidation = FALSE))
    class(pdata) <- c("PolyData", class(pdata));

  return(pdata);
}

#==============================================================================
plotLines <- function(polys, xlim = NULL, ylim = NULL, projection = FALSE,
                      plt = c(0.11, 0.98, 0.12, 0.88), polyProps = NULL,
                      lty = NULL, col = NULL, bg = 0, axes = TRUE,
                      tckLab = TRUE, tck = 0.014, tckMinor = 0.5 * tck, ...)
  # The only layout graphics parameter accepted by 'plotLines()' is 'plt'.
  # Since the function changes the plot region to satisfy the aspect ratio,
  # it is assigned a default value to improve consistency among plots.
{
  r <- .plotMaps(polys=polys, xlim=xlim, ylim=ylim, projection=projection,
                 plt=plt, polyProps=polyProps, border=NULL, lty=lty, col=col,
                 density=NULL, angle=NULL, bg=bg, axes=axes, tckLab=tckLab,
                 tck=tck, tckMinor=tckMinor, isType="lines", ...);

  if (!is.null(r) &&
      !is.PolyData(r, fullValidation = FALSE))
    class(r) <- c("PolyData", class(r));

  invisible(r);
}

#==============================================================================
plotMap <- function(polys, xlim = NULL, ylim = NULL, projection = TRUE,
                    plt = c(0.11, 0.98, 0.12, 0.88), polyProps = NULL,
                    border = NULL, lty = NULL, col = NULL, colHoles = NULL,
                    density = NA, angle = NULL, bg = 0, axes = TRUE,
                    tckLab = TRUE, tck = 0.014, tckMinor = 0.5 * tck, ...)
  # The only layout graphics parameter accepted by 'plotMap()' is 'plt'.
  # Since the function changes the plot region to satisfy the aspect ratio,
  # it is assigned a default value to improve consistency among plots.
{
  r <- .plotMaps(polys=polys, xlim=xlim, ylim=ylim, projection=projection,
                 plt=plt, polyProps=polyProps, border=border, lty=lty, col=col,
                 colHoles=colHoles, density=density, angle=angle, bg=bg,
                 axes=axes, tckLab=tckLab, tck=tck, tckMinor=tckMinor,
                 isType="polygons", ...);

  if (!is.null(r) &&
      !is.PolyData(r, fullValidation = FALSE))
    class(r) <- c("PolyData", class(r));

  invisible(r);
}

#==============================================================================
plotPoints <- function(data, xlim = NULL, ylim = NULL, projection = FALSE,
                       plt = c(0.11, 0.98, 0.12, 0.88), polyProps = NULL,
                       cex = NULL, col = NULL, pch = NULL, axes = TRUE,
                       tckLab = TRUE, tck = 0.014, tckMinor = 0.5 * tck, ...)
  # The only layout graphics parameter accepted by 'plotPoints()' is 'plt'.
  # Since the function changes it to satisfy the aspect ratio, it's assigned
  # a default value.
{
  r <- .plotMaps(polys=data, xlim=xlim, ylim=ylim, projection=projection,
                 plt=plt, polyProps=polyProps, border=NULL, lty=NULL, col=col,
                 density=NULL, angle=NULL, bg=NULL,
                 cex=cex, pch=pch, # both become part of '...'
                 axes=axes, tckLab=tckLab, tck=tck, tckMinor=tckMinor,
                 isType="points", ...);

  if (!is.null(r) &&
      !is.PolyData(r, fullValidation = FALSE))
    class(r) <- c("PolyData", class(r));

  invisible(r)
}

#==============================================================================
plotPolys <- function(polys, xlim = NULL, ylim = NULL, projection = FALSE,
                      plt = c(0.11, 0.98, 0.12, 0.88), polyProps = NULL,
                      border = NULL, lty = NULL, col = NULL, colHoles = NULL,
                      density = NA, angle = NULL, bg = 0, axes = TRUE,
                      tckLab = TRUE, tck = 0.014, tckMinor = 0.5 * tck, ...)
  # The only layout graphics parameter accepted by 'plotPolys()' is 'plt'.
  # Since the function changes it to satisfy the aspect ratio, it's assigned
  # a default value.
{
  r <- .plotMaps(polys=polys, xlim=xlim, ylim=ylim, projection=projection,
                 plt=plt, polyProps=polyProps, border=border, lty=lty, col=col,
                 colHoles=colHoles, density=density, angle=angle, bg=bg,
                 axes=axes, tckLab=tckLab, tck=tck, tckMinor=tckMinor,
                 isType="polygons", ...);

  if (!is.null(r) &&
      !is.PolyData(r, fullValidation = FALSE))
    class(r) <- c("PolyData", class(r));

  invisible(r)
}

#==============================================================================
print.EventData <- function(x, ...)
{
  if (exists("PBSprint") && (PBSprint || .PBSmapEnv$PBSprint)) {
    print(summary.EventData(x), ...);
  } else {
    print(data.frame(unclass(x)));
  }
  invisible(x);
}

#==============================================================================
print.LocationSet <- function(x, ...)
{
  if (exists("PBSprint") && (PBSprint || .PBSmapEnv$PBSprint)) {
    print(summary.LocationSet(x), ...);
  } else {
    print(data.frame(unclass(x)));
  }
  invisible(x);
}

#==============================================================================
print.PolyData <- function(x, ...)
{
  if (exists("PBSprint") && (PBSprint || .PBSmapEnv$PBSprint)) {
    print(summary.PolyData(x), ...);
  } else {
    print(data.frame(unclass(x)));
  }
  invisible(x);
}

#==============================================================================
print.PolySet <- function(x, ...)
{
  if (exists("PBSprint") && (PBSprint || .PBSmapEnv$PBSprint)) {
    print(summary.PolySet(x), ...);
  } else {
    print(data.frame(unclass(x)));
  }
  invisible(x);
}

#==============================================================================
print.summary.PBS <- function(x, ...)
{
  str <-
    c(x$type, "",
      "", "",
      "Records         : ", ifelse(is.null(x$records),"NULL", x$records),
      "  Contours      : ", ifelse(x$type == "EventData",
                                   "NA",
                                   ifelse(is.null(x$contours), "NULL",
                                          x$contours)),
      "    Holes       : ", ifelse(x$type != "PolySet",
                                   "NA",
                                   ifelse(is.null(x$holes), "NULL", x$holes)),
      "  Events        : ", ifelse(x$type == "EventData",
                                   ifelse(is.null(x$events), "NULL", x$events),
                                   "NA"),
      "    On boundary : ", ifelse(x$type != "LocationSet",
                                   "NA",
                                   ifelse(is.null(x$boundary), "NULL",
                                          x$boundary)),
      "", "",
      "Ranges", "",
      "  X             : ", ifelse(is.null(x$range.x), "NULL",
                                   paste("[",
                                         paste(range(x$range.x),collapse=", "),
                                         "]", sep="")),
      "  Y             : ", ifelse(is.null(x$range.y), "NULL",
                                   paste("[",
                                         paste(range(x$range.y),collapse=", "),
                                         "]", sep="")),
      "", "",
      "Attributes", "",
      "  Projection    : ", ifelse(is.null(x$attr.projection), "NULL",
                                   x$attr.projection),
      "  Zone          : ", ifelse(is.null(x$attr.zone), "NULL",
                                   x$attr.zone),
      "", "",
      "Extra columns   : ", ifelse(is.null(x$columns), "NULL", x$columns));

  for (i in seq(1, by=2, length.out=(length(str)/2))) {
    cat (str[i], str[i+1], "\n", sep="");
  }

  invisible(x);
}

#refocusWorld---------------------------2015-02-11
#  Refocus the 'worldLL'/'worldLLhigh' data sets.
#--------------------------------------------NB/RH
refocusWorld <- function(polys, xlim = NULL, ylim = NULL) {
  # define a local function to assist with shifting
  .shiftRegion <- function(polys, shift = 0) {
    npolys <- polys;

    npolys$xmin  <- polys$xmin + (360 * shift);
    npolys$xmax  <- polys$xmax + (360 * shift);
    npolys$shift <- polys$shift + shift;

    return (npolys);
  }

  # validate input
  polys <- .validatePolySet(polys);
  if (is.character(polys))
    stop(paste("Invalid PolySet 'polys'.\n", polys, sep=""));

  # default values
  if (is.null(xlim))
    xlim <- range(polys$X)
  if (is.null(ylim))
    ylim <- range(polys$Y)

  # validate/adjust xlim
  if (abs(diff(xlim)) > 360)
    stop("Invalid 'xlim' range: cannot exceed 360 degrees.\n");
  warn <- FALSE;
  while (min(xlim) < -360) {
    xlim <- xlim + 360;
    warn <- TRUE;
  }
  while (max(xlim) > 360) {
    xlim <- xlim - 360;
    warn <- TRUE;
  }
  if (warn)
    warning(paste("Shifted 'xlim' to ", paste(xlim, collapse=","),
                  ": use new limits in future plots.\n", sep=""));

  # save the attributes of the data frame
  attrNames <- setdiff(names(attributes(polys)),
                       c("names", "row.names", "class"));
  attrValues <- attributes(polys)[attrNames];

  # build indices in case both PID *and* SID exist
  polys$IDs <- .createIDs(polys, cols=c("PID", "SID"));

  # build a data frame with the ranges
  xrng <- lapply(split(polys[, c("X")], polys$IDs), range)
  yrng <- lapply(split(polys[, c("Y")], polys$IDs), range)
  ppid <- lapply(split(polys[, c("PID")], polys$IDs), unique)
  rng <- data.frame(IDs=names(xrng),
                    PID=unlist(ppid),
                    xmin=unlist(xrng)[c(TRUE,FALSE)],
                    xmax=unlist(xrng)[c(FALSE,TRUE)],
                    ymin=unlist(yrng)[c(TRUE,FALSE)],
                    ymax=unlist(yrng)[c(FALSE,TRUE)],
                    shift=0);

  # 'rng' currently contains one row per polygon; to wrap the region,
  # let's make the range [-360,360] contain *every* polygon that exists
  # within that range -- doing so will require replicating polygons;
  # for "shifted" polygons, we'll add a "shift" field to the data frame
  nrng <- rng;  # 'nrng' will be the new range data frame

  # go right until every polygon's X is greater than the range
  srng <- .shiftRegion (rng, 1);
  while (any(srng$xmin < 360)) {
    # typically one iteration
    nrng <- rbind (nrng, srng);
    srng <- .shiftRegion (srng, 1);
  }

  # go left until every polygon's X is less than the range
  srng <- .shiftRegion (rng, -1);
  while (any(srng$xmax > -360)) {
    # typically one iteration
    nrng <- rbind(nrng, srng);
    srng <- .shiftRegion (srng, -1);
  }

  # we've added polygons that are completely outside of the range;
  # remove those polygons
  nrng <- nrng[nrng$xmax > -360 & nrng$xmin < 360, ]

  # identify all of the polygons within the region of interest (x/ylim)
  nrng <- nrng[!(nrng$xmax < xlim[1] | nrng$xmin > xlim[2] |
                 nrng$ymax < ylim[1] | nrng$ymin > ylim[2]), ];

  # eliminate duplicated polys for now, i.e., Antarctica;
  # in the future (TO DO), we should probably merge them into a single
  # polygon; some code that might help...
  #   nrng <- nrng[order(nrng$IDs), ]
  #   dups <- duplicated(nrng$IDs)
  #   dups <- dups | c(dups[-1], FALSE)
  #   drng <<- nrng[dups, ]
  #   nrng <- nrng[!dups, ]
  # but for now, let's simply...
  dups <- duplicated(nrng$IDs);
  if (any(dups)) {
    warning(paste("Removed duplicates of following polygons (Antarctica?): ",
                  paste(unique(nrng$PID[dups]), collapse=", "),
                  "\n", sep=""));
    nrng <- nrng[!dups, ];
  }

  # given relevant PIDs and their shift factors, let's create a new
  # data frame containing only those polygons, where we've shifted the X
  # values
  idx <- split(1:nrow(polys), polys$IDs) # to extract relevant polygons
  len <- lapply(idx, length)           # num vertices in each polygon

  # indices for relevant polys
  idx <- unlist(idx[as.character(nrng$IDs)])
  # num vertices in relevant polys
  len <- unlist(len[as.character(nrng$IDs)])

  # build the return value
  npolys <- polys[idx, setdiff (names (polys), "IDs")];
  npolys$X <- npolys$X + rep(nrng$shift * 360, len);
  row.names(npolys) <- NULL;

  # restore the attributes
  attributes(npolys) <- c(attributes(npolys), attrValues);
  attr(npolys,"rf.xlim") = xlim # RH added for use in importGSHHS

  invisible(npolys);
}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~refocusWorld

#==============================================================================
summary.EventData <- function(object, ...)
{
  z <- list();
  z$type <- "EventData";
  z$records <- nrow(object);
  z$events <- length(unique(object$EID));
  z$columns <- paste(setdiff(names(object), c("EID", "X", "Y")), collapse=", ");
  z$range.x <- range(object$X);
  z$range.y <- range(object$Y);
  z$attr.projection <- attr(object, "projection");
  z$attr.zone <- attr(object, "zone");
  class(z) <- "summary.PBS";

  z;
}

#==============================================================================
summary.LocationSet <- function(object, ...)
{
  # create 'IDX' vector used later for 'unique'/'duplicated'
  IDX <- .createIDs(object, cols = c("PID", "SID"));

  z <- list();
  z$type <- "LocationSet";
  z$events <- length(unique(object$EID));
  z$contours <- length(unique(IDX));
  z$records <- nrow(object);
  z$boundary <- sum(object$Bdry);
  z$attr.projection <- attr(object, "projection");
  z$attr.zone <- attr(object, "zone");
  z$columns <- paste(setdiff(names(object), c("EID", "PID", "SID", "Bdry")),
                     collapse=", ");
  class(z) <- "summary.PBS";

  z;
}

#==============================================================================
summary.PolyData <- function(object, ...)
{
  # create 'IDX' vector used later for 'unique'/'duplicated'
  IDX <- .createIDs(object, cols = c("PID", "SID"));

  z <- list();
  z$type <- "PolyData";
  z$records <- nrow(object);
  z$contours <- length(unique(IDX));
  z$columns <- paste(setdiff(names(object), c("PID", "SID")), collapse=", ");
  z$attr.projection <- attr(object, "projection");
  z$attr.zone <- attr(object, "zone");
  class(z) <- "summary.PBS";

  z;
}

#==============================================================================
summary.PolySet <- function(object, ...)
{
  # create 'IDX' vector used later for 'unique'/'duplicated'
  # create 'holes' (Boolean) vector to count how many holes exist
  IDX <- .createIDs(object, cols=c("PID", "SID"));
  if (is.element("SID", names(object))) {
    idxFirst <- which(!duplicated(IDX));
    idxLast <- c((idxFirst-1)[-1], length(IDX));
    holes <- object$POS[idxFirst] > object$POS[idxLast];
  } else {
    holes <- NULL;
  }

  z <- list();
  z$type <- "PolySet";
  z$records <- nrow(object);
  z$contours <- length(unique(IDX));
  z$holes <- sum(holes);
  z$range.x <- range(object$X);
  z$range.y <- range(object$Y);
  z$attr.projection <- attr(object, "projection");
  z$attr.zone <- attr(object, "zone");
  z$columns <- paste(setdiff(names(object), c("PID", "SID", "POS", "X", "Y")),
                     collapse=", ");
  class(z) <- "summary.PBS";

  z;
}

#==============================================================================
thickenPolys <- function(polys, tol = 1, filter = 3, keepOrig = TRUE,
                         close = TRUE)
{
  polys <- .validatePolySet(polys);
  if (is.character(polys))
    stop(paste("Invalid PolySet 'polys'.\n", polys, sep=""));

  # PART 1 --------------------------------------------------------------------
  # - estimate 'outCapacity'
  # - this should use the same distance formula as the underlying C code;
  #   otherwise, the estimation may be wrong

  # create an 'idx' vector for working with 'polys'
  idx <- .createIDs(polys, cols = c("PID", "SID"));

  # if closing polygons, ensure they all close before we calculate
  # 'outCapacity'; close them all, even if already closed, as a
  # precaution
  if (close) {
    # find start and end of each polygon
    starts <- which(!duplicated(idx));
    # special case of a single polygon
    if (length(starts) == 1) {
      ends <- length(idx);
    } else {
      ends <- c((starts[-1]) - 1, length(idx));
    }

    # create a vector of [start #1, end #1, start #2, end #2, ...]
    toDupe <- vector();
    toDupe[seq(1, by=2, length=length(starts))] <- starts;
    toDupe[seq(2, by=2, length=length(ends))] <- ends;

    # determine next avail. PID
    newPID <- max(polys$PID) + 1;

    # create a new PolySet to use in these calculations
    newRows <- polys[toDupe, ];
    newRows$PID <- rep((newPID:(newPID + length(starts) - 1)), each=2)
    tempPolys <- rbind(polys, newRows);

    # update 'idx'; it won't matter that the new 'idx' values
    # are a different format since they're using new PIDs (i.e., they
    # don't need to match anything); they just need to be relative to
    # earlier 'idx' value, and they will be since we're using a 'paste'
    if (is.element("SID", names(newRows))) {
      idx <- c(idx, paste(newRows$PID, newRows$SID));
    } else {
      idx <- c(idx, newRows$PID);
    }
  } else {
    tempPolys <- polys;
  }

  # compute a distance vector
  projection <- attr(polys, "projection");
  if (!is.null(projection) && !is.na(projection)
      && ((projection == "UTM") || (projection == "LL")
          || (projection == 1))) {
    dist <- .calcDist(tempPolys);
  }
  else {
    stop(paste(
"Invalid projection attribute.  Supported projections include \"LL\",",
"\"UTM\", and 1.\n"));
  }
  # filter the distances between polygons (make them not count...)
  dist[((which(!duplicated(idx)) - 1)[-1])] <- 0;

  if (keepOrig == TRUE) {
    # since original points are kept, calculate the capacity based on how
    # many points will be added between the original points
    outCapacity <- sum(ceiling(dist / tol));
  } else {
    # since original points are not kept, calculate the length of each
    # polygon, and then determine how many points will be required to
    # represent that polygon
    dist <- unlist(lapply(split(dist, idx), sum), use.names = FALSE);
    outCapacity <- sum(ceiling(dist / tol));
  }

  # add a fudge factor
  outCapacity <- outCapacity + 1000;

  # PART 2 --------------------------------------------------------------------
  # - thicken the PolySet

  pLL <- 0;
  pUTM <- 1;

  projection <- attributes(polys)$projection;
  if (!is.null(projection) && !is.na(projection) && projection == "LL") {
    proj <- pLL;
  } else {
    proj <- pUTM;
  }

  # save the attributes of the data frame (.validatePolySet returns a data
  # frame); in this case, SAVE ALL ATTRIBUTES since we are "modifying" an
  # existing data set
  attrNames <- setdiff(names(attributes(polys)),
                       c("names", "row.names", "class"));
  attrValues <- attributes(polys)[attrNames];

  inRows <- nrow(polys);

  # create the data structures that the C function expects
  if (!is.element("SID", names(polys))) {
    inID <- c(polys$PID, integer(length = inRows), polys$POS);
  } else {
    inID <- c(polys$PID, polys$SID, polys$POS);
  }
  inXY <- c(polys$X, polys$Y);

  # call the C function
  results <- .C("thickenPolys",
                inID = as.integer(inID),
                inXY = as.double(inXY),
                inVerts = as.integer(inRows),
                tolerance = as.double(tol),
                filter = as.integer(filter),
                units = as.integer(proj),
                keepOrig = as.integer(keepOrig),
                close = as.integer(close),
                outID = integer(3 * outCapacity),
                outXY = double(2 * outCapacity),
                outRows = as.integer(outCapacity),
                outStatus = integer(1),
                PACKAGE = "PBSmapping");
  # note: outRows is set to how much space is allocated -- the C function
  #       should take this into consideration

  if (results$outStatus == 1) {
    stop(
"Insufficient physical memory for processing.\n");
  }
  if (results$outStatus == 2) {
    stop(paste(
"Insufficient memory allocated for output.  Please upgrade to the latest",
"version of the software, and if that does not fix this problem, please",
"file a bug report.\n",
               sep = "\n"));
  }

  # determine the number of rows in the result
  outRows <- as.vector(results$outRows);

  # extract the data from the C function results
  if (outRows > 0) {
    d <- data.frame(PID = results$outID[1:outRows],
                    SID = results$outID[(outCapacity+1):(outCapacity+outRows)],
                    POS = results$outID[(2*outCapacity+1):(2*outCapacity+outRows)],
                    X = results$outXY[1:outRows],
                    Y = results$outXY[(outCapacity+1):(outCapacity+outRows)]);

    if (!is.element("SID", names(polys)))
      d$SID <- NULL;

    # restore the attributes (including projection and zone)
    attributes(d) <- c(attributes(d), attrValues);

    if (!is.PolySet(d, fullValidation = FALSE))
      class(d) <- c("PolySet", class(d));

    return(d);
  } else {
    return(NULL);
  }
}

#==============================================================================
thinPolys <- function(polys, tol = 1, filter = 3)
{
  polys <- .validatePolySet(polys);
  if (is.character(polys))
    stop(paste("Invalid PolySet 'polys'.\n", polys, sep=""));

  pLL <- 0;
  pUTM <- 1;

  projection <- attributes(polys)$projection;
  if (!is.null(projection) && !is.na(projection) && projection == "LL") {
    proj <- pLL;
  } else {
    proj <- pUTM;
  }

  # save the attributes of the data frame (.validatePolySet returns a data
  # frame); in this case, SAVE ALL ATTRIBUTES since we are "modifying" an
  # existing data set
  attrNames <- setdiff(names(attributes(polys)),
                       c("names", "row.names", "class"));
  attrValues <- attributes(polys)[attrNames];

  inRows <- nrow(polys);
  outCapacity <- as.integer(inRows);

  # create the data structures that the C function expects
  if (!is.element("SID", names(polys))) {
    inID <- c(polys$PID, integer(length = inRows), polys$POS);
  } else {
    inID <- c(polys$PID, polys$SID, polys$POS);
  }

  # adjust units to what the C function expects
  if (proj == pLL) {
    # degrees >> micro-degrees
    inXY <- c(polys$X, polys$Y) * 10^6;
  } else {
    # kilometers >> meters
    inXY <- c(polys$X, polys$Y) * 10^3;
  }

  # call the C function
  results <- .C("thinPolys",
                inID = as.integer(inID),
                inXY = as.integer(inXY),
                inVerts = as.integer(inRows),
                tolerance = as.double(tol),
                filter = as.integer(filter),
                units = as.integer(proj),
                outID = integer(3 * outCapacity),
                outXY = integer(2 * outCapacity),
                outRows = as.integer(outCapacity),
                outStatus = integer(1),
                PACKAGE = "PBSmapping");
  # note: outRows is set to how much space is allocated -- the C function
  #       should take this into consideration

  if (results$outStatus == 1) {
    stop(
"Insufficient physical memory for processing.\n");
  }
  if (results$outStatus == 2) {
    stop(paste(
"Insufficient memory allocated for output.  Please upgrade to the latest",
"version of the software, and if that does not fix this problem, please",
"file a bug report.\n",
               sep = "\n"));
  }

  # determine the number of rows in the result
  outRows <- as.vector(results$outRows);

  # extract the data from the C function results
  if (outRows > 0) {
    if (proj == pLL) {
      results$outXY <- results$outXY / 10^6;
    } else {
      results$outXY <- results$outXY / 10^3;
    }

    d <- data.frame(PID = results$outID[1:outRows],
                    SID = results$outID[(outCapacity+1):(outCapacity+outRows)],
                    POS = results$outID[(2*outCapacity+1):(2*outCapacity+outRows)],
                    X = results$outXY[1:outRows],
                    Y = results$outXY[(outCapacity+1):(outCapacity+outRows)]);

    if (!is.element("SID", names(polys)))
      d$SID <- NULL;

    # restore the attributes (including projection and zone)
    attributes(d) <- c(attributes(d), attrValues);

    if (!is.PolySet(d, fullValidation = FALSE))
      class(d) <- c("PolySet", class(d));

    return(d);
  } else {
    return(NULL);
  }
}

#importShapefile------------------------2013-12-19
# importShapefile (Nick Boers)
# This function has several slow parts:
# 1) conversion of the matrix (verts) to X and Y columns in a data frame
# 2) 'lapply's to create POS columns
# Changes 2008-07-15:
#   Nick's loop to extract data from 'shapeList' has been replaced
#     by Rowan's series of 'sapply' calls.
#   Rowan added check for polygons with 0 vertices.
# 2012-04-04: Rowan created function 'placeHoles'
#   to place holes under correct solids.
#--------------------------------------------NB/RH
importShapefile <- function (fn, readDBF=TRUE, projection=NULL, zone=NULL,
     placeholes=FALSE, minverts=3)
{
	# initialization
	.checkRDeps("importShapefile", c("maptools", "foreign"))
	# call to normalizePath added to perform ~ expansion; otherwise,
	# pathnames beginning with a ~ fail in the later call to
	# Rshapeget
	fn <- normalizePath(fn, mustWork=FALSE)
	fn <- .getBasename(fn, "shp")

	# test for the required '.shx' file
	shxFile <- paste(fn, ".shx", sep="")
	if (!file.exists(shxFile))
		stop(paste(
		"Cannot find the index file (\"", shxFile, "\") required to import\n",
		"the shapefile.\n", sep=""))

	# read shapefile
	eval(parse(text="shapeList <- .Call(\"Rshapeget\",as.character(fn),as.logical(FALSE),PACKAGE = \"maptools\")"))
	if (length(shapeList) < 1)
		stop("The shapefile is empty or an error occurred while importing.\n")
	shpType=unique(sapply(shapeList,function(x){x$shp.type}))
	if (length(shpType) != 1)
		stop ("Supports only a single shape type per shapefile.\n")
	nVerts=sapply(shapeList,function(x){x$nVerts})
	v0=is.element(nVerts,0) # any shapefiles with 0 vertices?
	if (any(v0==TRUE)) {
		nVerts=nVerts[!v0]; shapeList=shapeList[!v0] }
	shpID=sapply(shapeList,function(x){x$shpID})
	nParts=sapply(shapeList,function(x){x$nParts})
	pStarts=sapply(shapeList,function(x){x$Pstart},simplify=FALSE)
	if (length(pStarts)!=length(nParts) && all((nParts==sapply(pStarts,length))!=TRUE))
		stop ("Mismatch in 'nParts' and 'pStarts'.\n")
	pStarts=unlist(pStarts)
	v1=unlist(sapply(shapeList,function(x){x$verts[,1]},simplify=FALSE))
	v2=unlist(sapply(shapeList,function(x){x$verts[,2]},simplify=FALSE))
	verts=cbind(v1,v2)

	# Keep track of parents and children
	PC=pStarts
	zP=is.element(PC,0); PC[zP]=1; PC[!zP]=0

	# reformat results
	#if (shpType == 3 || shpType == 5) {  # PolySet
	if (shpType %in% c(3,13,23, 5,15,25)) {  # PolyLine, PolyLineZ, PolyLineM, Polygon, PolygonZ, PolygonM
		# create preliminary PID/SID columns
		PID <- rep(1:(length(unique(shpID))), times=nParts)
		SID <- unlist(lapply(split(nParts, 1:(length(nParts))), "seq"))

		# to determine the number of vertices in each part, we divide the problem
		# into two cases:
		# 1) last component/hole of each polygon: the total vertices in the polygon
		#    less the starting POS of that last component/hole
		# 2) otherwise: use a "diff" on the starting POS's of each part
		lastComp <- rev(!duplicated(rev(PID)))
		nv <- vector()
		nv[lastComp] <- rep(nVerts, times=nParts)[lastComp] - pStarts[lastComp]
		nv[!lastComp] <- diff(pStarts)[diff(pStarts) > 0]

		# create PID/SID columns
		PID <- rep(PID, times=nv)
		SID <- rep(SID, times=nv)
		# create POS column; we'll fix the ordering for holes later
		POS <- unlist(lapply(split(nv, 1:(length(nv))), "seq"))
		# build the data frame
		df <- data.frame(PID=PID, SID=SID, POS=POS, X=verts[, 1], Y=verts[, 2])

		#if (shpType == 5) {
		if (shpType %in% c(5,15,25)) {
			# PolySet: polygons: reorder vertices for holes
			or <- .calcOrientation (df)
			# where "orientation" == -1, we need to reverse the POS
			or$solid <- is.element(or$orientation,1); or$hole <- !or$solid
			if (any(or$hole)) {
				or$nv <- nv
				toFix <- rep(or$hole, times=or$nv)
				o <- or$nv[or$hole]
				newPOS <- unlist(lapply(lapply(split(o, 1:length(o)), "seq"), "rev"))
				df[toFix, "POS"] <- newPOS  }

			if (placeholes) {
				# Fix to the problem where ArcPew does not put solid shapes before holes
				df = placeHoles(df,minVerts=minverts)
			}
		}
		class(df) <- c("PolySet", class(df))
	#} else if (shpType == 1) {  # EventData
	} else if (shpType %in% c(1,11,21)) {  # Point, PointZ, PointM
		EID <- 1:(length(unique(shpID)))
		df <- data.frame(EID=EID, X=verts[, 1], Y=verts[, 2])
		class(df) <- c("EventData", class(df))
	} else {
		stop ("Shape type not supported.\n");
	}

	# "cbind" the DBF for EventData or attach as an attribute for PolySets:
	# According to the "ESRI Shapefile Technical Description", any set of fields
	# may be present in the DBF file.
	# The (relevant) requirements are:
	#   (1) one record per shape feature (i.e., per PID or EID), and
	#   (2) same order as in shape (*.shp) file.
	dbfFile <- paste(fn, ".dbf", sep="")
	if (readDBF && !file.exists(dbfFile)) {
		warning(paste(
		"The argument 'readDBF' is true but the attribute database\n",
		"(\"", dbfFile, "\") does not exist.\n", sep=""))
	} else if (readDBF) {
		dbf <- foreign::read.dbf (dbfFile)
		if (shpType == 1) {  # EventData
			if (nrow(df) != nrow(dbf)) {
				warning(paste(
				"The shapefile and its associated DBF do not contain the",
				"same number of records. DBF ignored.\n", sep = "\n"))
				return (df)
			}
			df.class = class(df)
			df <- cbind(df, dbf)
			class(df) <- df.class
		} else if (shpType == 3 || shpType == 5) {
			# add index to result
			dbf <- cbind(1:nrow(dbf), dbf)
			names(dbf)[1] <- "PID"
			class(dbf) <- c("PolyData", class(dbf))
			attr(df, "PolyData") <- dbf
		}
	# At this point, shpTypes != 1, 3, 5 caused the "stop" above; we do not
	# need an "else" to check here
	}
	attr(df,"parent.child") = PC
	attr(df,"shpType") = shpType
	prjFile <- paste(fn, ".prj", sep="")
	if (file.exists(prjFile)) {
		prj=scan(prjFile,what="character",quiet=TRUE)
		prj=prj[!is.element(prj,"")][1]
		if (length(prj)==0 || is.na(prj) || is.null(prj)) prj="Unknown" }
	else prj="Unknown"
	attr(df,"prj") = prj
	xmlFile <- paste(fn, ".shp.xml", sep="")
	if (file.exists(xmlFile)) {
		xml=readLines(xmlFile); attr(df,"xml") = xml }

	if (regexpr("GEO",prj)>0 | regexpr("Degree",prj)>0) proj="LL"
	else if (regexpr("PROJ",prj)>0 && regexpr("UTM",prj)>0) proj="UTM"
	else proj=1
	attr(df,"projection")=proj

	if (proj=="UTM" && any(df$X>500))
		{df$X=df$X/1000; df$Y=df$Y/1000}
	if (!is.null(zone))
		attr(df, "zone") <- zone
	if (!is.null(projection))
		attr(df,"projection")=projection
	return (df) }
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~importShapefile


#2008-08-25----------------------------------NB/RH
# 'projection' should equal "LL" (for longitude/latitude),
# "UTM" (for Universal Transverse Mercator), or anything else (for X/Y).
# 'main', 'sub', 'xlab', 'ylab' may all be part of '...'
#-------------------------------------------------
.addLabels <- function(projection = NULL, ...) {
  dots <- list(...);
  main=dots$main;  sub=dots$sub;  xlab=dots$xlab;  ylab=dots$ylab;

  # set label defaults (if necessary)
  if (!is.null(projection) && !is.na(projection) && projection == "UTM") {
    if (is.null(xlab)) xlab <- "UTM Easting (km)";
    if (is.null(ylab)) ylab <- "UTM Northing (km)";
  } else if (!is.null(projection) && !is.na(projection) && projection == "LL") {
    if (is.null(xlab)) xlab <- expression("Longitude" ~ group("(", degree, ")"))
    if (is.null(ylab)) ylab <- expression("Latitude" ~ group("(", degree, ")"))
  } else {
    if (is.null(xlab)) xlab <- "X";
    if (is.null(ylab)) ylab <- "Y";
  }
  if (is.null(main)) main <- "";
  if (is.null(sub))  sub  <- "";
  # 'xlab'/'ylab' cannot equal NULL at this point
  title (main=main, sub=sub, xlab=xlab, ylab=ylab) # cannot add dots - conflicts with formal arguments
  invisible(NULL);
}

#makeTopography-------------------------2009-07-30
# Make topography data suitable for contourLines().
#-----------------------------------------------RH
makeTopography <- function (dat, digits=2, func=NULL) {
	if (!is.data.frame(dat) || length(dat) < 3) {
		stop("'dat' must be a data frame with at least three columns.\n") }
	N   = nrow(dat)
	X   = round(dat[[1]],digits); Y  = round(dat[[2]],digits)
	ID  = complex(length.out=N,real=X,imaginary=Y)
	tmp = split(dat[[3]], ID)
	if (is.null(func)) func = function(x){mean(x,na.rm=TRUE)} # mean value
	z1  = unlist(sapply(tmp,func,simplify=FALSE))
	id  = as.complex(names(z1))
	x1  = Re(id); y1 = Im(id)
	x   = sort(unique(x1)); y = sort(unique(y1))
	nr  = length(x); nc = length(y); nz=nr*nc
	z   = matrix(NA,nrow=nr,ncol=nc,dimnames=list(x,y))
	for (i in x) { # populate matrix row by row
		xx=as.character(i); zx=is.element(x1,i)
		yy=as.character(y1[zx])
		zz=z1[zx]; z[xx,yy]=zz  }
	result = list()
	result$x = x; result$y = y; result$z = z
	return(result) }

#placeHoles-----------------------------2012-04-05
# Place holes under correct solids.
#-----------------------------------------------RH
placeHoles = function(polyset, minVerts=3) {

	is.in =function(hole,solid) {
		eval(parse(text=".Call(\"R_point_in_polygon_sp\", as.numeric(hole$X), as.numeric(hole$Y), as.numeric(solid$X), as.numeric(solid$Y), PACKAGE = \"sp\")"))
	}
	# initialization
	.checkRDeps("placeHoles", c("sp"))
	newpoly = list(); newpoly = NULL
	pid = sort(unique(polyset$PID))
	for (i in pid) {
		ipoly = polyset[is.element(polyset$PID,i),]
		orien = .calcOrientation (ipoly)
		zsols = is.element(orien$orientation,1)
		ssids = orien$SID[zsols]
		hsids = orien$SID[!zsols]
		for (s in ssids) {
			orphans = adopted = NULL
			spoly = ipoly[is.element(ipoly$SID,s),]
			newpoly = rbind(newpoly, spoly)
			if (length(hsids)==0) next
			for (h in hsids) {
				hpoly = ipoly[is.element(ipoly$SID,h),]
				if (nrow(hpoly) < minVerts) next
				io = is.in(hpoly,spoly)
				if (any(io==0)) orphans = c(orphans,h)
				else adopted = c(adopted, h)
			}
			newpoly = rbind(newpoly, ipoly[is.element(ipoly$SID,adopted),])
			hsids = orphans
		}
	}
	Natts = names(attributes(polyset))
	natts = names(attributes(newpoly))
	Datts = setdiff(Natts,natts)
	Aatts = attributes(polyset)[Datts]
	attributes(newpoly) = c(attributes(newpoly),Aatts)
	return(newpoly) }
#---------------------------------------placeHoles


