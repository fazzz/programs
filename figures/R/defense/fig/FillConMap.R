FillConMap <- function (x = seq(0, 1, length.out = nrow(z)), y = seq(0, 1, 
    length.out = ncol(z)), z, xlim = range(x, finite = TRUE), 
    ylim = range(y, finite = TRUE), zlim = range(z, finite = TRUE), 
    levels = pretty(zlim, nlevels), nlevels = 20, color.palette = cm.colors, 
    col = color.palette(length(levels) - 1), plot.title, plot.axes,
    key.title, key.axes, asp = NA, xaxs = "i", yaxs = "i", las = 1,
    axes = TRUE,wxaxis = "T",wyaxis = "T",xrange=xrange,yrange=yrange,
    xrange.axis=xrange.axis,yrange.axis=yrange.axis,labelflag="T", frame.plot = axes, ...) 
{
    if (missing(z)) {
        if (!missing(x)) {
            if (is.list(x)) {
                z <- x$z
                y <- x$y
                x <- x$x
            }
            else {
                z <- x
                x <- seq.int(0, 1, length.out = nrow(z))
            }
        }
        else stop("no 'z' matrix specified")
    }
    else if (is.list(x)) {
        y <- x$y
        x <- x$x
    }
    if (any(diff(x) <= 0) || any(diff(y) <= 0)) 
        stop("increasing 'x' and 'y' values expected")

    plot.new()
    plot.window(xlim, ylim, "", xaxs = "i", yaxs = "i",las=1)
    if (!is.matrix(z) || nrow(z) <= 1L || ncol(z) <= 1L) 
        stop("no proper 'z' matrix specified")
    if (!is.double(z)) 
        storage.mode(z) <- "double"

    .Internal(filledcontour(as.double(x), as.double(y), z, as.double(levels), 
        col = col))
    
    if (missing(plot.axes)) {                         
        if (axes) {                                   
            title(main = "", xlab = "", ylab = "")    
            cat(xrange,'\n')                          
            if (wxaxis=="T") {
              Axis(x, side = 1,xaxp=xrange.axis)
            }

            if (wyaxis=="T") {
              Axis(y, side = 2,las=1,yaxp=yrange.axis)
            }
        }                                             
    }                                                 
    else plot.axes                                    
    
   if (frame.plot)                                 
     box()                                        

    if (labelflag=="T") {
      if (missing(plot.title))
        title(...)
      else plot.title
    }
    
    invisible()
}
