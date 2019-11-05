.plotLayer <- function (layer, loc, asp = NULL, add = FALSE)
{
    axis <- which(!is.na(loc))
    indices <- alist(i=, j=, k=)
    indices[axis] <- loc[axis]
    data <- do.call("[", c(layer["image"],indices))
    dims <- dim(data)
    
    if (add)
        data <- replace(data, which(data==0), NA)
    
    if (is.null(asp))
        asp <- dims[2] / dims[1]
    
    oldPars <- par(mai=c(0,0,0,0), bg=layer$colours[1])
    on.exit(par(oldPars))
    image(data, col=layer$colours, axes=FALSE, asp=asp, add=add, zlim=layer$window)
}

#' @export
view <- function (..., x = c(0,0,0), y = NA, z = NA, world = TRUE, orientation = "LAS", interactive = base::interactive(), crosshairs = TRUE, labels = TRUE, infoPanel = defaultInfoPanel)
{
    layers <- list(...)
    layerExpressions <- sapply(substitute(list(...)), as.character)[-1]
    nLayers <- length(layers)
    if (nLayers == 0)
        stop("At least one image must be specified")
    
    for (i in seq_len(nLayers))
    {
        if (!inherits(layers[[i]], "viewLayer"))
        {
            layers[[i]] <- layer(layers[[i]])
            layers[[i]]$label <- layerExpressions[i]
        }
        orientation(layers[[i]]$image) <- orientation
    }
    
    dims <- dim(layers[[1]]$image)
    fov <- dims * pixdim(layers[[1]]$image)
    
    if (length(x) > 1)
        point <- x
    else
        point <- c(x, y, z)
    
    if (world)
        point <- round(worldToVoxel(point, xform(layers[[1]]$image)))
    
    labelText <- list(c("P","A","I","S"), c("R","L","I","S"), c("R","L","P","A"))
    
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
        point[point < 1] <- 1
        point[point > dims] <- dims[point > dims]
        voxelCentre <- (point - 1) / (dims - 1)
        
        starts <- ends <- numeric(0)
        
        # Plot the info panel first so that we have some handle on the coordinate system when we use locator()
        layout(matrix(c(2,3,4,1),nrow=2,byrow=TRUE))
        
        data <- lapply(layers, function(layer) do.call("[", c(list(layer$image),as.list(point))))
        imageNames <- sapply(layers, "[[", "label")
        infoPanel(point, data, imageNames)
        
        for (i in 1:3)
        {
            inPlaneAxes <- setdiff(1:3, i)
            loc <- replace(rep(NA,3), i, point[i])
            
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
                text(c(0.1*width[1]+region[1],0.9*width[1]+region[1],0.5*width[2]+region[3],0.5*width[2]+region[3]), c(0.5*width[1]+region[1],0.5*width[1]+region[1],0.1*width[2]+region[3],0.9*width[2]+region[3]), labels=labelText[[i]])
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
            point[2:3] <- round(adjustedPoint * (dims[2:3] - 1)) + 1
        }
        else if (nextPoint[1] > ends[5] && nextPoint[2] > ends[6])
        {
            adjustedPoint <- (nextPoint-ends[5:6]) / (ends[5:6]-starts[5:6]) * (ends[3:4]-starts[3:4]) + starts[3:4]
            point[c(1,3)] <- round(adjustedPoint * (dims[c(1,3)] - 1)) + 1
        }
        else
            point[1:2] <- round(nextPoint * (dims[1:2] - 1)) + 1
    }
    
    invisible(NULL)
}

#' @export
layer <- function (image, scale = "grey", min = NULL, max = NULL)
{
    label <- deparse(substitute(image))
    image <- as.array(retrieveNifti(image))
    
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
    
    return (structure(list(image=image, label=label, colours=colours, window=window), class="viewLayer"))
}

#' @export
defaultInfoPanel <- function (point, data, imageNames)
{
    usingQuartz <- isTRUE(names(dev.cur()) == "quartz")
    quitInstructions <- paste(ifelse(usingQuartz,"Press Esc","Right click"), "to exit", sep=" ")
    
    plot(NA, xlim=c(0,1), ylim=c(0,1), xlab="", ylab="", xaxt="n", yaxt="n", bty="n", main=paste("Location: (", paste(point,collapse=","), ")", sep=""))
    nImages <- min(4, length(imageNames))
    yLocs <- 0.95 - cumsum(c(0,rep(c(0.1,0.13),nImages)))
    yLocs[length(yLocs)] <- -0.05
    labels <- quitInstructions
    for (i in seq_len(nImages))
    {
        labels <- c(labels, {
            if (is.numeric(data[[i]]))
                as.character(signif(mean(data[[i]]),6))
            else
                data[[i]]
        }, imageNames[i])
    }
    text(0.5, yLocs, rev(labels), col=c(rep(c("white","red"),nImages),"grey70"), cex=pmin(1,1/strwidth(rev(labels))), xpd=TRUE)
}
