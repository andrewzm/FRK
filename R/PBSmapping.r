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
.validateData <- function(data, className,
                          requiredCols = NULL, requiredAttr = NULL,
                          noFactorCols = NULL, noNACols = NULL, keyCols = NULL,
                          numericCols = NULL)
    # An element of noNACols and keyCols will only be used if it exists in the
    # data.  To ensure it exists in the data, make it a requiredCol.
{
    # convert matrix to data frame
    ## if (is.matrix(data)) {
    ##     data <- .mat2df(data);
    ## }

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


#==============================================================================
is.PolyData <- function(x, fullValidation = TRUE)
{
    # if (fullValidation) {
    #     msg <- .validatePolyData(x)
    #     if (is.character(msg))
    #         return (FALSE)
    # }

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

