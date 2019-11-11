.plotLayer <- function (layer, loc, asp = NULL, add = FALSE)
{
    axis <- which(!is.na(loc))
    indices <- alist(i=, j=, k=, t=1, u=1, v=1, w=1)
    indices[axis] <- loc[axis]
    data <- do.call("[", c(layer["image"], indices[seq_len(ndim(layer$image))], list(drop=FALSE)))
    dims <- dim(data)[-axis]
    dim(data) <- dims
    
    if (is.null(asp))
        asp <- dims[2] / dims[1]
    
    if (inherits(layer$image, "rgbArray"))
    {
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
        if (add)
            data <- replace(data, which(data==0), NA)
        oldPars <- par(mai=c(0,0,0,0), bg=layer$colours[1])
        on.exit(par(oldPars))
        image(data, col=layer$colours, axes=FALSE, asp=asp, add=add, zlim=layer$window)
    }
}

#' @export
view <- function (..., point = NULL, radiological = FALSE, interactive = base::interactive(), crosshairs = TRUE, labels = TRUE, infoPanel = defaultInfoPanel)
{
    layers <- list(...)
    layerExpressions <- trimws(sapply(substitute(list(...)), deparse, control=NULL, nlines=1)[-1], whitespace="[ \"]")
    nLayers <- length(layers)
    if (nLayers == 0)
        stop("At least one image must be specified")
    
    originalXform <- NULL
    orientation <- ifelse(radiological, "LAS", "RAS")
    for (i in seq_len(nLayers))
    {
        if (!inherits(layers[[i]], "viewLayer"))
        {
            layers[[i]] <- layer(layers[[i]])
            layers[[i]]$label <- layerExpressions[i]
        }
        if (i == 1)
            originalXform <- xform(layers[[i]]$image)
        header <- niftiHeader(layers[[i]]$image)
        if (ndim(layers[[i]]$image) > 2 && (header$qform_code > 0 || header$sform_code > 0))
            orientation(layers[[i]]$image) <- orientation
    }
    
    baseImage <- layers[[1]]$image
    reorientedXform <- xform(baseImage)
    ndim <- ndim(baseImage)
    dims <- c(dim(baseImage), rep(1,max(0,3-ndim)))[1:3]
    fov <- dims * c(pixdim(baseImage), rep(1,max(0,3-ndim)))[1:3]
    
    if (is.null(point) && any(origin(baseImage) > 1))
        point <- round(origin(baseImage))
    else if (is.null(point))
        point <- round(dims / 2)
    reorientedPoint <- round(worldToVoxel(voxelToWorld(point, originalXform), reorientedXform))
    
    positiveLabels <- unlist(strsplit(orientation, ""))
    negativeLabels <- c(R="L", A="P", S="I", L="R", P="A", I="S")[positiveLabels]
    
    oldPars <- par(bg="black", col="white", fg="white", col.axis="white", col.lab="white", col.main="white")
    oldOptions <- options(locatorBell=FALSE, preferRaster=TRUE)
    on.exit({
        par(oldPars)
        options(oldOptions)
        if (interactive)
            dev.off()
    })
    
    repeat
    {
        reorientedPoint[reorientedPoint < 1] <- 1
        reorientedPoint[reorientedPoint > dims] <- dims[reorientedPoint > dims]
        voxelCentre <- (reorientedPoint - 1) / (dims - 1)
        
        point <- round(worldToVoxel(voxelToWorld(reorientedPoint, reorientedXform), originalXform))
        
        starts <- ends <- numeric(0)
        
        # Plot the info panel first so that we have some handle on the coordinate system when we use locator()
        if (ndim == 2)
        {
            starts <- ends <- rep(0:1, 2)
            layout(matrix(c(2,1), nrow=1))
        }
        else
            layout(matrix(c(2,3,4,1), nrow=2, byrow=TRUE))
        
        data <- lapply(layers, function(layer) {
            indices <- alist(i=, j=, k=, t=, u=, v=, w=)
            indices[seq_along(point)] <- reorientedPoint
            do.call("[", c(list(layer$image), indices[seq_len(ndim(layer$image))]))
        })
        infoPanel(point, data, sapply(layers,"[[","label"))
        
        for (i in 1:3)
        {
            if (ndim == 2 && i < 3)
                next
            
            inPlaneAxes <- setdiff(1:3, i)
            loc <- replace(rep(NA,3), i, reorientedPoint[i])
            
            for (j in seq_along(layers))
                .plotLayer(layers[[j]], loc, asp=fov[inPlaneAxes[2]]/fov[inPlaneAxes[1]], add=(j>1))
            
            region <- par("usr")
            starts <- c(starts, region[c(1,3)])
            ends <- c(ends, region[c(2,4)])
            width <- c(region[2]-region[1], region[4]-region[3])
            
            if (crosshairs)
            {
                halfVoxelWidth <- 0.5 / (dims[inPlaneAxes] - 1)
                lines(rep(voxelCentre[inPlaneAxes[1]],2), c(-halfVoxelWidth[2],1+halfVoxelWidth[2]), col="red")
                lines(c(-halfVoxelWidth[1],1+halfVoxelWidth[1]), rep(voxelCentre[inPlaneAxes[2]],2), col="red")
            }
            
            if (labels)
            {
                # Order is left, right, bottom, top
                currentLabels <- c(negativeLabels[inPlaneAxes[1]], positiveLabels[inPlaneAxes[1]], negativeLabels[inPlaneAxes[2]], positiveLabels[inPlaneAxes[2]])
                text(c(0.1*width[1]+region[1],0.9*width[1]+region[1],0.5*width[2]+region[3],0.5*width[2]+region[3]), c(0.5*width[1]+region[1],0.5*width[1]+region[1],0.1*width[2]+region[3],0.9*width[2]+region[3]), labels=currentLabels)
            }
        }
        
        if (!interactive)
            break
        
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

#' @export
layer <- function (image, scale = "grey", min = NULL, max = NULL)
{
    label <- deparse(substitute(image))
    image <- as.array(retrieveNifti(image))
    
    if (inherits(image, "rgbArray"))
        colours <- window <- NULL
    else
    {
        if (is.character(scale) && length(scale) == 1)
            colours <- switch(scale, grey=gray(0:99/99), gray=gray(0:99/99), greyscale=gray(0:99/99), grayscale=gray(0:99/99), heat=heat.colors(100), rainbow=rainbow(100,start=0.7,end=0.1), shades::gradient(scale,100))
        else
            colours <- scale
    
        if (is.null(min))
            min <- image$cal_min
        if (is.null(max))
            max <- image$cal_max
    
        window <- c(min, max)
        if (any(is.na(window)) || (min == max))
        {
            window <- quantile(image[is.finite(image)], c(0.025,0.975), na.rm=TRUE)
            if (diff(window) > abs(mean(window)))
                window[which.min(abs(window))] <- 0
        }
    
        image[image < window[1]] <- window[1]
        image[image > window[2]] <- window[2]
    }
    
    return (structure(list(image=image, label=label, colours=colours, window=window), class="viewLayer"))
}

#' @export
defaultInfoPanel <- function (point, data, labels)
{
    escapeToQuit <- isTRUE(names(dev.cur()) %in% c("quartz","RStudioGD"))
    quitInstructions <- paste(ifelse(escapeToQuit,"Press Esc","Right click"), "to exit", sep=" ")
    
    plot(NA, xlim=c(0,1), ylim=c(0,1), xlab="", ylab="", xaxt="n", yaxt="n", bty="n", main=paste("Location: (", paste(point,collapse=","), ")", sep=""))
    nImages <- min(4, length(labels))
    yLocs <- 0.95 - cumsum(c(0,rep(c(0.1,0.13),nImages)))
    yLocs[length(yLocs)] <- -0.05
    text <- quitInstructions
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
