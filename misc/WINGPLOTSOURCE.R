####################

#PrcCrds <- wings.LM[,4:99]

WingPlot <- function( PrcCrds, wingcol="blue", lines=TRUE, add=FALSE, wingcex=0.75, winglwd=1, winglty=1, wingpch=1, LM_highlit=FALSE, LM_pch=1 ){

	odds <- seq(1, 95, 2)
	evens <- seq(2, 96, 2)
	if ( class(PrcCrds) == "data.frame" ){
		LMx <- PrcCrds[,c(odds)]
		LMy <- PrcCrds[,c(evens)]
		LMx.mean <- colMeans( LMx, na.rm=TRUE)
		LMy.mean <- colMeans( LMy, na.rm=TRUE)
	} 	else {
			LMx.mean <- PrcCrds[c(odds)]
			LMy.mean <- PrcCrds[c(evens)]
	}

	if (add == FALSE){
	plot( LMx.mean, LMy.mean, xlim=c(-0.25,0.25), ylim=c(-0.25,0.25), cex=0.3, type="n", xaxt="n", yaxt="n", ann=FALSE)
	points( LMx.mean, LMy.mean, cex=wingcex, pch=wingpch)
		}
	from <- c( 6, 47, 5, 29, 5, 28, 27, 26, 25, 4, 24, 23, 3, 22, 2, 21, 20, 19, 1, 18, 17, 16, 15, 14, 12, 33, 32, 7, 31, 8, 48, 37, 9, 36, 8, 35, 34, 30, 11, 42, 10, 41, 40, 39, 38, 10, 11, 46, 45, 44, 43 )
	to <- c( 47, 5, 29, 30, 28, 27, 26, 25, 4, 24, 23, 3, 22, 2, 21, 20, 19, 1, 18, 17, 16, 15, 14, 13, 33, 32, 7, 31, 1, 48, 7, 9, 36, 8, 35, 34, 2, 11, 42, 10, 41, 40, 39, 38, 3, 9, 46, 45, 44, 43, 4 )

	if(lines==TRUE){
		for (i in 1:length(from)) {
    	lines( c(LMx.mean[from[i]], LMx.mean[to[i]]), c(LMy.mean[from[i]], LMy.mean[to[i]]), col=wingcol, lty=winglty, lwd=winglwd )
		}	}

	if( LM_highlit==TRUE ) {
	  lmco <- c(1:13, 30)
	  points( LMx.mean[ lmco ], LMy.mean[ lmco ], cex=wingcex, pch=LM_pch )
	}
	}

#WingPlot(PrcCrds, 48, wingcol="dark green")

comment(WingPlot) <- c("WingPlot takes a matrix of Proc Coords and 'joins the dots' to plot the wing outline and veins. It plots mean points for each LM and semi-LM and connects them with lines when lines=TRUE (default). Can accept a 'wingcol' argument if you don't like blue")


##################

WingBlur <- function ( PrcCrds, wingcol="#add8e632", winglwd=3, winglty=1, wingpch=1, wingcex=0.9, meanline=TRUE, meancol="blue", meanlty=1, meanlwd=2 ) {
	# set up required coord vectors
	from <- c( 6, 47, 5, 29, 5, 28, 27, 26, 25, 4, 24, 23, 3, 22, 2, 21, 20, 19, 1, 18, 17, 16, 15, 14, 12, 33, 32, 7, 31, 8, 48, 37, 9, 36, 8, 35, 34, 30, 11, 42, 10, 41, 40, 39, 38, 10, 11, 46, 45, 44, 43 )
	to <- c( 47, 5, 29, 30, 28, 27, 26, 25, 4, 24, 23, 3, 22, 2, 21, 20, 19, 1, 18, 17, 16, 15, 14, 13, 33, 32, 7, 31, 1, 48, 7, 9, 36, 8, 35, 34, 2, 11, 42, 10, 41, 40, 39, 38, 3, 9, 46, 45, 44, 43, 4 )
	x.vals <- seq(1, 95, 2)		#these pull out the odd...
	y.vals <- seq(2, 96, 2)		#...and even columns...
	shape.X <- PrcCrds[,x.vals]	# make x coord matrix
	shape.Y <- PrcCrds[,y.vals]	# and one for y
	LMx.mean <- colMeans(PrcCrds)[x.vals]	# make mean x vector
	LMy.mean <- colMeans(PrcCrds)[y.vals]	# and one for y

	# set up the plot
	plot( LMx.mean, LMy.mean, 
	      xlim = c(-0.25,0.25), ylim = c(-0.25,0.25), 
	      type = "n", xaxt = "n", yaxt = "n", ann = FALSE, asp = 1)

	# nested loops to draw each line segment for each row in shape.X matrices
	for (j in 1:nrow(shape.X)){
		for (i in 1:length(from)) {
		lines( c(shape.X[j,from[i]], shape.X[j,to[i]]), c(shape.Y[j,from[i]], shape.Y[j,to[i]]), col=wingcol, lty=winglty, lwd=3 )
		}	}

	# loop to draw the mean lines if requested
	if (meanline ==TRUE) {
		for (i in 1:length(from)) {
    	lines( c(LMx.mean[from[i]], LMx.mean[to[i]]), c(LMy.mean[from[i]], LMy.mean[to[i]]), col=meancol, lty=meanlty, lwd=meanlwd )
    	}	}

	# finally, this drop in the mean LM/semi-LM points
	points(LMx.mean, LMy.mean, cex = wingcex, pch = wingpch)
}


comment(WingBlur) <- c("This draws a wing outline&veins for each row of the matrix 'PrcCrds' - can be individuals, line means, whatever... the default setting are for translucent lines so that the coloured area represents that phenotypic space sampled by the wings measured. Mean shape overlaid as points + outline, but can be suppressed")

##################

WingEffect <- function ( meanshape, effectplus, effectminus, wingcol = c("black", "black", "red"), winglwd = c(2, 2, 2), winglty = c(1, 1, 1), scale.factor=2, meanline=FALSE, wingcex = c(0,0,0), winglabel = "", add = FALSE, wingframe = T, scale.display = T, wingpoints = F, wingpch = c(16,1,1) ) {

	# set up required coord vectors
		from <- c( 6, 47, 5, 29, 5, 28, 27, 26, 25, 4, 24, 23, 3, 22, 2, 21, 20, 19, 1, 18, 17, 16, 15, 14, 12, 33, 32, 7, 31, 8, 48, 37, 9, 36, 8, 35, 34, 30, 11, 42, 10, 41, 40, 39, 38, 10, 11, 46, 45, 44, 43 )	# from and to indices drawn from the 'segment' file
		to <- c( 47, 5, 29, 30, 28, 27, 26, 25, 4, 24, 23, 3, 22, 2, 21, 20, 19, 1, 18, 17, 16, 15, 14, 13, 33, 32, 7, 31, 1, 48, 7, 9, 36, 8, 35, 34, 2, 11, 42, 10, 41, 40, 39, 38, 3, 9, 46, 45, 44, 43, 4 )
		x.vals <- seq(1, 95, 2)	# indices to specify x and y coords
		y.vals <- seq(2, 96, 2)
		LMx.mean <- meanshape[x.vals]	# apply x & y indices to mean data
		LMy.mean <- meanshape[y.vals]
		# x & y indices to effectplus vector
		x.plus <- LMx.mean + (effectplus[x.vals] * scale.factor)
		y.plus <- LMy.mean + (effectplus[y.vals] * scale.factor)
		# x & y indices to effectminus vector
		x.minus <- LMx.mean - (effectminus[x.vals] * scale.factor)
		y.minus <- LMy.mean - (effectminus[y.vals] * scale.factor)

	# plot mean shape
		if ( add == FALSE ) { # opens blank plot
		plot( LMx.mean, LMy.mean, xlim = c(-0.22, 0.27), ylim = c(-0.3, 0.27), 
		      type = "n", xaxt = "n", yaxt = "n", ann = FALSE, frame = wingframe)
		}
		if (meanline == TRUE){
		# this for loop draws line segments for the mean shape
		for (i in 1:length(from)) {
	    	lines( c(LMx.mean[from[i]], LMx.mean[to[i]]), 
	    	       c(LMy.mean[from[i]], LMy.mean[to[i]]), 
	    	       col=wingcol[1], lty=winglty[1], lwd=winglwd[2])
	    	}	}
	    if (wingpoints == TRUE) {
	    	points( x.plus, y.plus, cex = wingcex[2], pch = wingpch[2], col = wingcol[2])
	    	points( x.minus, y.minus, cex = wingcex[3], pch = wingpch[3], col = wingcol[3])
	    }

		# this drop in the mean LM/semi-LM points
		points(LMx.mean, LMy.mean, cex = wingcex[1], pch = wingpch[1], col = wingcol[1])

	# this for loop draws the mean+effect line segments
	for (i in 1:length(from)) {
		lines( c(x.plus[from[i]], x.plus[to[i]]), c(y.plus[from[i]], y.plus[to[i]]), col=wingcol[2], lty=winglty[2], lwd=winglwd[2] )
		lines( c(x.minus[from[i]], x.minus[to[i]]), c(y.minus[from[i]], y.minus[to[i]]), col=wingcol[3], lty=winglty[3], lwd=winglwd[3] )
		}
	# adding an annotation of the effect scaling factor
	if (scale.display == T) {
		text( 0.025, -0.215, paste("scale factor =", as.character(scale.factor)) )
		}
	text(0, 0.25, winglabel) # optional centre-top label
}


comment(WingEffect) <- c("This makes mean ± an effect of choice wingplots - mean LMs/semi-LMs connected by outline&veins lines, with superimposed outline&veins lines for a specifed modifier effect. e.g. specify a vector of effect coefficients from a model to plot lines for mean+effect and mean-effect")

####################
require(animation)

WingMovie <- function ( meanshape, effect, wingcol = c("black", "red"), winglwd = c(2, 2), winglty = c(1, 1), wingcex = 1, scale.factor = 2, effectlines = FALSE, meanline = TRUE, winglabel="" ) {

	# set up required coord vectors
		from <- c( 6, 47, 5, 29, 5, 28, 27, 26, 25, 4, 24, 23, 3, 22, 2, 21, 20, 19, 1, 18, 17, 16, 15, 14, 12, 33, 32, 7, 31, 8, 48, 37, 9, 36, 8, 35, 34, 30, 11, 42, 10, 41, 40, 39, 38, 10, 11, 46, 45, 44, 43 )	# from and to indices drawn from the 'segment' file
		to <- c( 47, 5, 29, 30, 28, 27, 26, 25, 4, 24, 23, 3, 22, 2, 21, 20, 19, 1, 18, 17, 16, 15, 14, 13, 33, 32, 7, 31, 1, 48, 7, 9, 36, 8, 35, 34, 2, 11, 42, 10, 41, 40, 39, 38, 3, 9, 46, 45, 44, 43, 4 )
		x.vals <- seq(1, 95, 2)	# indices to specify x and y coords
		y.vals <- seq(2, 96, 2)
		LMx.mean <- meanshape[x.vals]	# apply x & y indices to mean data
		LMy.mean <- meanshape[y.vals]
	# set up X & Y coords from effect vector
		x.loadings <- effect[c(x.vals)]		#...to become the X...
		y.loadings <- effect[c(y.vals)]		#...and Y coordinates
	# set up stepping vectors
		multiplier <- c(seq(-scale.factor, scale.factor, (scale.factor/10)), seq(scale.factor, -scale.factor, -(scale.factor/10)))	# multiplier is just a sequence from -2 to 2 and back to scale the frames
		x.plus <- (LMx.mean + scale.factor *(x.loadings))  #vector of X mean + loadings
		x.minus <- (LMx.mean - scale.factor *(x.loadings))  #vector of X mean - loadings
		y.plus <- (LMy.mean + scale.factor *(y.loadings))  #vector of Y mean + loadings
		y.minus <- (LMy.mean - scale.factor *(y.loadings))  #vector of Y mean - loadings

	# set up animation 'frames'
	for (j in 1:length(multiplier)) {
		plotX <- as.numeric(LMx.mean + (x.loadings * multiplier[j]))
		plotY <- as.numeric(LMy.mean + (y.loadings * multiplier[j]))
		# plotX&Y are vectors of mean coords + (loadings * scale)

		plot( LMx.mean, LMy.mean, xlim = c(-0.22,0.27), ylim = c(-0.22,0.27), cex=wingcex, col=wingcol[1], xaxt="n", yaxt="n", ann=FALSE) 	# sets up plot, adds mean points
		text( 0.025, -0.215, paste("scale factor =", as.character(scale.factor)) ) #labels with scale factor
		text(0, 0.25, winglabel) # optional centre-top label

		for (i in 1:length(from)) {
			if (meanline==TRUE) {
				lines( c(LMx.mean[from[i]], LMx.mean[to[i]]), c(LMy.mean[from[i]], LMy.mean[to[i]]), col=wingcol[1], lty=winglty[1], lwd=winglwd[1] )
				}
		    	lines( c(plotX[from[i]], plotX[to[i]]), c(plotY[from[i]], plotY[to[i]]), col=wingcol[2], lty=winglty[2], lwd=winglwd[2] )	# draws the outline segments for each 'from-to' pair
    		if (effectlines==TRUE){
   				for (j in 1:length(LMx.mean)) {
					lines( c(x.plus[j], x.minus[j]), c(y.plus[j], y.minus[j]), col="black", lwd=1)
				}	}		# this loop draws in the vectors of change across the mean for each LM/sLM
    	}
	}
}


comment(WingMovie) <- c("This makes mean ± an effect of choice wingplots - mean LMs/semi-LMs connected by outline&veins lines, with an animated line tracking from +effect to -effect [*scale factor]. Colours and line types may be specified as a vector e.g. wingcols=c('blue', 'red')... Call with saveGIF(expr=WingMovie( meanshape, effect ), interval=0.1, filename='name.gif') ")


####################
require(animation)

WingMovieSNP <- function ( SNP, meanline=FALSE ) {

	LMx.mean <- c( 0.052, 0.221, 0.236, 0.149, -0.121, -0.195, -0.004, 0.003, -0.088, -0.086, -0.150, -0.170, -0.197, -0.182, -0.147, -0.099, -0.049, 0.001, 0.096, 0.141, 0.185, 0.234, 0.218, 0.187, 0.101, 0.053, 0.004, -0.058, -0.155, -0.196, 0.028, -0.054, -0.111, 0.151, 0.078, -0.041, -0.140, 0.175, 0.118, 0.052, -0.018, -0.115, 0.095, 0.043, -0.017, -0.076, -0.155, 0.002 )
	LMy.mean <- c( -0.112, -0.028, 0.021, 0.092, 0.060, 0.049, -0.059, -0.017, 0.000, 0.017, 0.020, -0.012, -0.022, -0.070, -0.103, -0.121, -0.126, -0.122, -0.103, -0.086, -0.062, -0.004, 0.054, 0.077, 0.100, 0.100, 0.095, 0.084, 0.035, 0.025, -0.081, -0.040, -0.024, -0.021, -0.017, -0.010, 0.002, 0.029, 0.031, 0.031, 0.028, 0.018, 0.078, 0.069, 0.058, 0.043, 0.062, -0.039 )

	x.vals <- seq(1, 95, 2)		#these pull out the odd...
	y.vals <- seq(2, 96, 2)		#...and even columns...
	x.loadings <- loadings[,c(x.vals)]		#...to become the X...
	y.loadings <- loadings[,c(y.vals)]		#...and Y coordinates

multiplier <- c(seq(-2, 2, 0.2), seq(2, -2, -0.2))
#multiplier is just a sequence from -2 to 2 and back to scale the frames

	x.plus <- (LMx.mean + x.loadings[SNP,])  #vector of X mean + loadings
	x.minus <- (LMx.mean - x.loadings[SNP,])  #vector of X mean - loadings
	y.plus <- (LMy.mean + y.loadings[SNP,])  #vector of Y mean + loadings
	y.minus <- (LMy.mean - y.loadings[SNP,])  #vector of Y mean - loadings

	#'from' and 'to' are paired vectors - the wing outline is made up of lines that are plotted FROM from[n] TO to[n]...
		from <- c( 6, 47, 5, 29, 5, 28, 27, 26, 25, 4, 24, 23, 3, 22, 2, 21, 20, 19, 1, 18, 17, 16, 15, 14, 12, 33, 32, 7, 31, 8, 48, 37, 9, 36, 8, 35, 34, 30, 11, 42, 10, 41, 40, 39, 38, 10, 11, 46, 45, 44, 43 )
		to <- c( 47, 5, 29, 30, 28, 27, 26, 25, 4, 24, 23, 3, 22, 2, 21, 20, 19, 1, 18, 17, 16, 15, 14, 13, 33, 32, 7, 31, 1, 48, 7, 9, 36, 8, 35, 34, 2, 11, 42, 10, 41, 40, 39, 38, 3, 9, 46, 45, 44, 43, 4 )

	par(mar=c(1, 1, 1, 1))	#margins not needed as we have no labels

for (j in 1:length(multiplier)) {	 #i.e. for each frame of the animation
	plotX <- as.numeric(LMx.mean + (x.loadings[SNP,] * multiplier[j]))
	plotY <- as.numeric(LMy.mean + (y.loadings[SNP,] * multiplier[j]))
	#plotX&Y are vectors of mean coords + (loadings * scale)

	plot( LMx.mean, LMy.mean, xlim=c(-0.22,0.27), ylim=c(-0.22,0.27), cex=1, col="black", xaxt="n", yaxt="n", ann=FALSE) 	#sets up the plot and adds mean points
	#text( 0, 0.25, "brk")		#labels the plot with 'SNP'

		for (i in 1:length(from)) {
			if(meanline==TRUE){
		lines( c(LMx.mean[from[i]], LMx.mean[to[i]]), c(LMy.mean[from[i]], LMy.mean[to[i]]), col="blue", lty=1, lwd=3 )
					}
    	lines( c(plotX[from[i]], plotX[to[i]]), c(plotY[from[i]], plotY[to[i]]), col="red", lty=1, lwd=3 )		#draws the outline segments for each 'from-to' pair

#    	for (j in 1:length(LMx.mean)) {
#	lines( c(x.plus[j], x.minus[j]), c(y.plus[j], y.minus[j]), col="blue", lwd=2)
#	}		#this loop draws in the vectors of change across the mean for each LM/sLM
    	}
	}
}

comment(WingMovieSNP) <- "the function requires 'loadings', which should be a 88 column dataframe (eg: loadings <- chrsX1[,10:97] ), comprising the 44x&44y coordinates part of the output from CPReaderer. 'SNP' is provided as a column index."


#wingmovie(85, meanline=TRUE)

#saveGIF(expr=wingmovie(85, meanline=TRUE), interval=0.1, "test!.gif")



################
