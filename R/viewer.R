.State <- new.env()

.vector3Dim <- function (image)
{
    if (image$intent_code == 1007L && ndim(image) == 5L && dim(image)[5] == 3L)
        return (5L)
    else if (getOption("RNifti.d4vectors",FALSE) && ndim(image) == 4L && dim(image)[4] == 3L)
        return (4L)
    else
        return (NA_integer_)
}

.plotLayer <- function (layer, loc, asp = NULL, add = FALSE)
{
    vector3Dim <- .vector3Dim(layer$image)
    
    # Extract the data for the appropriate plane
    axis <- which(!is.na(loc))
    indices <- alist(i=, j=, k=, t=1, u=1, v=1, w=1)
    if (!is.na(vector3Dim))
        indices[[vector3Dim]] <- 1:3
    indices[[axis]] <- loc[axis]
    data <- do.call("[", c(layer["image"], indices[seq_len(ndim(layer$image))], list(drop=FALSE)))
    dims <- dim(data)[-c(axis, setdiff(4:7,vector3Dim))]
    dim(data) <- dims
    
    if (is.null(asp))
        asp <- dims[2] / dims[1]
    
    if (!is.na(vector3Dim))
    {
        # Show 3D vector data as coloured line segments
        if (!add)
        {
            oldPars <- par(mai=c(0,0,0,0))
            on.exit(par(oldPars))
            image(array(NA,dim=dims[1:2]), axes=FALSE, asp=asp, zlim=c(0,1))
        }
        valid <- which(apply(data, 1:2, function(x) !any(is.na(x)) && !all(x==0)), arr.ind=TRUE)
        colours <- as.character(rgbArray(abs(data[cbind(valid,1)]), abs(data[cbind(valid,2)]), abs(data[cbind(valid,3)]), max=layer$window[2]))
        data <- data[,,-axis]
        segments((valid[,1]-1)/(dims[1]-1) - data[cbind(valid,1)]/(2*dims[1]*layer$window[2]),
                 (valid[,2]-1)/(dims[2]-1) - data[cbind(valid,2)]/(2*dims[2]*layer$window[2]),
                 (valid[,1]-1)/(dims[1]-1) + data[cbind(valid,1)]/(2*dims[1]*layer$window[2]),
                 (valid[,2]-1)/(dims[2]-1) + data[cbind(valid,2)]/(2*dims[2]*layer$window[2]), col=colours, lwd=2)
    }
    else if (inherits(layer$image, "rgbArray"))
    {
        # RGB display is achieved by converting the data to an array of indices into a custom palette
        data <- as.character(structure(data, class="rgbArray"))
        palette <- unique(data)
        indices <- match(data, palette)
        dim(indices) <- dims
        oldPars <- par(mai=c(0,0,0,0))
        on.exit(par(oldPars))
        image(indices, col=palette, axes=FALSE, asp=asp, add=add, zlim=c(1,length(palette)))
    }
    else
    {
        # Other data is shown using standard image(), but zeroes are set to NA to make them transparent
        if (add)
            data <- replace(data, which(data==0), NA)
        oldPars <- par(mai=c(0,0,0,0), bg=layer$colours[1])
        on.exit(par(oldPars))
        image(data, col=layer$colours, axes=FALSE, asp=asp, add=add, zlim=layer$window)
    }
}

#' A basic 3D image viewer
#' 
#' This function displays one or more 2D or 3D images, with optional
#' click-to-navigate interactivity.
#' 
#' @param ... One or more images, or \code{"viewLayer"} objects, to display.
#'   These will be plotted in the order specified, so that the first image
#'   forms the lowest layer.
#' @param point A numeric vector giving the location to initially centre the
#'   view on. If crosshairs are in use, they will be placed at this point. For
#'   3D images, this parameter also determines the planes shown in each
#'   subview.
#' @param radiological Logical value. If \code{TRUE}, images will be displayed
#'   in the radiological convention whereby the left of the image is shown on
#'   the right; otherwise left is on the left.
#' @param interactive Logical value. If \code{TRUE}, the user can navigate
#'   around the image by repeatedly clicking on a new centre point; otherwise
#'   the view is fixed.
#' @param crosshairs Logical value, indicating whether crosshairs should be
#'   shown or not.
#' @param labels Logical value, indicating whether orientation labels should be
#'   shown or not. Ignored (defaulting to \code{FALSE}) if the image is 2D or
#'   orientation information is not available.
#' @param infoPanel A function of three arguments, which must produce a plot
#'   for the information panel of the view. \code{\link{defaultInfoPanel}} is
#'   the default, which shows the labels and values of each image at the
#'   current point. A \code{NULL} value will suppress the panel.
#' @param image The image being shown in this layer.
#' @param scale A character vector of colour values for the scale, or a single
#'   string naming a predefined scale: \code{"grey"} or \code{"gray"} for
#'   greyscale, \code{"heat"} for a heatmap, \code{"rainbow"} for a rainbow
#'   scale, or any of the scales defined in the \code{shades} package (see
#'   \code{?shades::gradient}, if that package is installed). A fixed colour
#'   can be used by wrapping a string in a call to \code{I}. Ignored for RGB
#'   images.
#' @param min,max The window minimum and maximum for the layer, i.e., the black
#'   and white points. These are ignored for RGB images. Otherwise, if
#'   \code{NULL}, the default, they are taken from the \code{cal_min} or
#'   \code{cal_max} NIfTI header fields. If either is \code{NA}, the image has
#'   no window stored in its header, or the two values are equal, then the 1st
#'   and 99th percentiles of the data are used, with values close to zero
#'   rounded to that extreme. If these values are still equal, the untrimmed
#'   range of the image data is used.
#' @param mask A optional mask array, which may be of lower dimensionality than
#'   the main image. If specified, this is converted to logical mode and pixels
#'   that evaluate \code{FALSE} will be set to \code{NA} for that layer,
#'   meaning they will not be plotted. This operation is performed last, and so
#'   will not affect auto-windowing.
#' @return \code{lyr} returns a list of class \code{"viewLayer"}, to be used
#'   in a view. \code{view} is called for its side-effect of showing a view.
#' 
#' @note Because of the way R's main run-loop interacts with graphics, it will
#'   not be possible to issue further commands to the terminal while
#'   interactive mode is enabled. Instructions for leaving this mode are shown
#'   by the default info panel; see also \code{\link{locator}}, which is the
#'   underlying core function.
#' 
#' @examples
#' im <- readNifti(system.file("extdata", "example.nii.gz", package="RNifti"))
#' view(im, interactive=FALSE)
#' view(lyr(im, max=800), interactive=FALSE)
#' view(lyr(im, mask=im<800), interactive=FALSE)
#' 
#' @author Jon Clayden <code@@clayden.org>
#' @seealso \code{\link{defaultInfoPanel}}, \code{\link{orientation}},
#'   \code{\link{locator}}
#' @export
view <- function (..., point = NULL, radiological = getOption("radiologicalView",FALSE), interactive = base::interactive(), crosshairs = TRUE, labels = TRUE, infoPanel = defaultInfoPanel)
{
    # Get the layers to display, and the expressions used to generate them
    layers <- list(...)
    layerExpressions <- sub("^[\\s\"]*(.+?)[\\s\"]*$", "\\1", sapply(substitute(list(...)), deparse, control=NULL, nlines=1)[-1], perl=TRUE)
    nLayers <- length(layers)
    if (nLayers == 0)
        stop("At least one image must be specified")
    
    # In general, data will need reordering for visual consistency
    originalXform <- NULL
    orientation <- ifelse(radiological, "LAS", "RAS")
    for (i in seq_len(nLayers))
    {
        # Images need to be converted to layer objects
        if (!inherits(layers[[i]], "viewLayer"))
        {
            layers[[i]] <- lyr(layers[[i]])
            layers[[i]]$label <- layerExpressions[i]
        }
        # The xform of the first image is used for indexing
        if (i == 1)
            originalXform <- xform(layers[[i]]$image)
        # If the image has orientation information and isn't 2D, reorient it as necessary
        header <- niftiHeader(layers[[i]]$image)
        if (ndim(layers[[i]]$image) > 2 && (header$qform_code > 0 || header$sform_code > 0))
            orientation(layers[[i]]$image) <- orientation
    }
    
    baseImage <- layers[[1]]$image
    reorientedXform <- xform(baseImage)
    ndim <- ndim(baseImage)
    dims <- c(dim(baseImage), rep(1,max(0,3-ndim)))[1:3]
    fov <- dims * c(pixdim(baseImage), rep(1,max(0,3-ndim)))[1:3]
    
    # Don't show labels if the base image is 2D or has no meaningful xform information
    if (ndim < 3L || attr(reorientedXform,"code") == 0L)
        labels <- FALSE
    
    # If no point is specified, use the origin if it's nontrivial, otherwise the centre of the image
    if (is.null(point) && any(origin(originalXform) > 1))
        point <- round(origin(originalXform))
    else if (is.null(point))
        point <- round(dims / 2)
    
    # Work out the point location in viewer space
    reorientedPoint <- round(worldToVoxel(voxelToWorld(point, originalXform), reorientedXform))
    
    positiveLabels <- unlist(strsplit(orientation, ""))
    negativeLabels <- c(R="L", A="P", S="I", L="R", P="A", I="S")[positiveLabels]
    
    # Set some graphics parameters, and make sure they get reset
    oldPars <- par(bg="black", col="white", fg="white", col.axis="white", col.lab="white", col.main="white")
    oldOptions <- options(locatorBell=FALSE, preferRaster=TRUE)
    .State$interactive <- interactive
    on.exit({
        par(oldPars)
        options(oldOptions)
        .State$interactive <- FALSE
    })
    
    repeat
    {
        reorientedPoint[reorientedPoint < 1] <- 1
        reorientedPoint[reorientedPoint > dims] <- dims[reorientedPoint > dims]
        voxelCentre <- (reorientedPoint - 1) / (dims - 1)
        
        # Work out the point location in source space
        point <- round(worldToVoxel(voxelToWorld(reorientedPoint, reorientedXform), originalXform))
        
        # The boundaries of each subview
        starts <- ends <- numeric(0)
        
        # Plot the info panel first so that we have some handle on the coordinate system when we use locator()
        if (ndim == 2)
        {
            starts <- ends <- rep(0:1, 2)
            if (is.null(infoPanel))
                layout(matrix(1))
            else
                layout(matrix(c(2,1), nrow=1))
        }
        else
            layout(matrix(c(2,3,4,1) - ifelse(is.null(infoPanel),1,0), nrow=2, byrow=TRUE))
        
        if (!is.null(infoPanel))
        {
            # For each layer, extract the data values corresponding to the current point, and pass them to the info panel function
            data <- lapply(layers, function(layer) {
                indices <- alist(i=, j=, k=, t=, u=, v=, w=)
                indices[seq_along(point)] <- reorientedPoint
                result <- do.call("[", c(list(layer$image), indices[seq_len(ndim(layer$image))]))
                if (inherits(layer$image, "rgbArray"))
                    return (as.character(structure(result, dim=c(1,length(result)), class="rgbArray")))
                else
                    return (result)
            })
            infoPanel(point, data, sapply(layers,"[[","label"))
        }
        
        for (i in 1:3)
        {
            # 2D images don't require three views
            if (ndim == 2 && i < 3)
                next
            
            inPlaneAxes <- setdiff(1:3, i)
            loc <- replace(rep(NA,3), i, reorientedPoint[i])
            
            # Plot each layer for the current plane
            for (j in seq_along(layers))
                .plotLayer(layers[[j]], loc, asp=fov[inPlaneAxes[2]]/fov[inPlaneAxes[1]], add=(j>1))
            
            # "Measure" the subview canvas
            region <- par("usr")
            starts <- c(starts, region[c(1,3)])
            ends <- c(ends, region[c(2,4)])
            width <- c(region[2]-region[1], region[4]-region[3])
            
            # Plot the crosshairs, if required
            if (crosshairs)
            {
                halfVoxelWidth <- 0.5 / (dims[inPlaneAxes] - 1)
                lines(rep(voxelCentre[inPlaneAxes[1]],2), c(-halfVoxelWidth[2],1+halfVoxelWidth[2]), col="red")
                lines(c(-halfVoxelWidth[1],1+halfVoxelWidth[1]), rep(voxelCentre[inPlaneAxes[2]],2), col="red")
            }
            
            # Plot the labels, if required
            if (labels)
            {
                # Order is left, right, bottom, top
                currentLabels <- c(negativeLabels[inPlaneAxes[1]], positiveLabels[inPlaneAxes[1]], negativeLabels[inPlaneAxes[2]], positiveLabels[inPlaneAxes[2]])
                text(c(0.1*width[1]+region[1],0.9*width[1]+region[1],0.5*width[2]+region[3],0.5*width[2]+region[3]), c(0.5*width[1]+region[1],0.5*width[1]+region[1],0.1*width[2]+region[3],0.9*width[2]+region[3]), labels=currentLabels)
            }
        }
        
        if (!interactive)
            break
        
        # Find the next point
        nextPoint <- locator(1)
        if (is.null(nextPoint))
            break
        
        # Coordinates are relative to the axial plot at this point
        nextPoint <- unlist(nextPoint)
        if (nextPoint[1] > ends[5] && nextPoint[2] <= ends[6])
            next
        else if (nextPoint[1] <= ends[5] && nextPoint[2] > ends[6])
        {
            adjustedPoint <- (nextPoint-c(starts[5],ends[6])) / (ends[5:6]-starts[5:6]) * (ends[1:2]-starts[1:2]) + starts[1:2]
            reorientedPoint[2:3] <- round(adjustedPoint * (dims[2:3] - 1)) + 1
        }
        else if (nextPoint[1] > ends[5] && nextPoint[2] > ends[6])
        {
            adjustedPoint <- (nextPoint-ends[5:6]) / (ends[5:6]-starts[5:6]) * (ends[3:4]-starts[3:4]) + starts[3:4]
            reorientedPoint[c(1,3)] <- round(adjustedPoint * (dims[c(1,3)] - 1)) + 1
        }
        else
            reorientedPoint[1:2] <- round(nextPoint * (dims[1:2] - 1)) + 1
    }
    
    invisible(NULL)
}

#' @rdname view
#' @export
lyr <- function (image, scale = "grey", min = NULL, max = NULL, mask = NULL)
{
    label <- deparse(substitute(image))
    image <- asNifti(image, internal=FALSE)
    
    if (inherits(image, "rgbArray"))
        colours <- window <- NULL
    else
    {
        colours <- NULL
        if (is.character(scale) && length(scale) == 1 && !inherits(scale,"AsIs"))
            colours <- try(switch(scale, grey=gray(0:99/99), gray=gray(0:99/99), greyscale=gray(0:99/99), grayscale=gray(0:99/99), heat=heat.colors(100), rainbow=rainbow(100,start=0.7,end=0.1), unclass(shades::gradient(scale,100))), silent=TRUE)
        
        # If the scale isn't named as per the switch() call above or recognised by shades::gradient(), the latter will have thrown an error, so treat the name as a literal colour instead
        # Ditto if it was surrounded by I() so as to make it of class "AsIs"
        if (is.null(colours) || inherits(colours, "try-error"))
            colours <- unclass(scale)
        
        if (is.null(min))
            min <- image$cal_min
        if (is.null(max))
            max <- image$cal_max
    
        window <- c(min, max)
        if (any(is.na(window)) || (min == max))
        {
            window <- quantile(image[is.finite(image)], c(0.01,0.99), na.rm=TRUE)
            if (diff(window) > abs(mean(window)))
                window[which.min(abs(window))] <- 0
            if (diff(window) == 0)
                window <- range(image, na.rm=TRUE)
            message("Setting window to (", signif(window[1],4), ", ", signif(window[2],4), ")")
        }
        
        if (is.na(.vector3Dim(image)))
        {
            image[image < window[1]] <- window[1]
            image[image > window[2]] <- window[2]
        }
    }
    
    if (!is.null(mask) && ndim(mask) > 1)
    {
        repDims <- setdiff(seq_len(ndim(image)), seq_len(ndim(mask)))
        if (length(repDims) > 0)
        {
            data <- apply(image, repDims, "[<-", !as.logical(mask), NA)
            dim(data) <- dim(image)
            image <- updateNifti(data, image)
        }
        else
            image[!as.logical(mask)] <- NA
    }
    
    return (structure(list(image=image, label=label, colours=colours, window=window), class="viewLayer"))
}

.quitInstructions <- function ()
{
    if (!isTRUE(.State$interactive))
        return ("")
    else
    {
        escapeToQuit <- isTRUE(names(dev.cur()) %in% c("quartz","RStudioGD"))
        return (paste(ifelse(escapeToQuit,"Press Esc","Right click"), "to leave interactive mode", sep=" "))
    }
}

#' Info panels for the built-in viewer
#' 
#' A default info panel for \code{\link{view}}, which shows the labels and
#' values of each image at the current point, and a panel suitable for
#' plotting four-dimensional time series images.
#' 
#' @param point A numeric vector giving the current point location.
#' @param data A list of data values for each image at the current point.
#'   Note that, for images of more than three dimensions, there will be more
#'   than one value per image.
#' @param labels A character vector of image labels.
#' 
#' @author Jon Clayden <code@@clayden.org>
#' @seealso \code{\link{view}}
#' @export
defaultInfoPanel <- function (point, data, labels)
{
    plot(NA, xlim=c(0,1), ylim=c(0,1), xlab="", ylab="", xaxt="n", yaxt="n", bty="n", main=paste("Location: [", paste(point,collapse=","), "]", sep=""))
    nImages <- min(4, length(labels))
    yLocs <- 0.95 - cumsum(c(0,rep(c(0.1,0.13),nImages)))
    yLocs[length(yLocs)] <- -0.05
    text <- .quitInstructions()
    for (i in seq_len(nImages))
    {
        text <- c(text, {
            if (is.numeric(data[[i]]) && length(data[[i]]) == 1)
                as.character(signif(data[[i]], 6))
            else if (is.numeric(data[[i]]))
                paste0(signif(data[[i]][1],6), ", ... (", length(data[[i]]), " values)")
            else
                data[[i]]
        }, labels[i])
    }
    text(0.5, yLocs, rev(text), col=c(rep(c("white","red"),nImages),"grey70"), cex=pmin(1,1/strwidth(rev(text))), xpd=TRUE)
}

#' @rdname defaultInfoPanel
#' @export
timeSeriesPanel <- function (point, data, labels)
{
    lengths <- sapply(data, length)
    suppressWarnings(range <- c(min(sapply(data,min,na.rm=TRUE)), max(sapply(data,max,na.rm=TRUE))))
    range[is.infinite(range)] <- 0
    plot(NA, xlim=c(1,max(lengths)), ylim=range, xlab=.quitInstructions(), ylab="intensity", col.lab="grey70", bty="n", main=paste("Location: [", paste(point,collapse=","), "]", sep=""))
    
    for (i in seq_along(data))
    {
        if (lengths[i] > 1)
            lines(1:lengths[i], data[[i]], col="red", lwd=2)
    }
}
