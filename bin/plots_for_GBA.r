# Extra functions
##
get_permuted_scores <- function(genes.labels, network, nFold){
        n = dim(network)[1]
        roc.permutes = matrix(0,ncol=nFold, nrow=dim(genes.labels)[2])

        for ( i in 1:nFold){
                permute_test = sample(n,n)
                permute      = neighbor.voting.CV(genes.labels, network[permute_test,], 3)
                roc.permutes[,i] = permute[,1]
        }
        return (roc.permutes)
}


## Plots the distribution (densities) of different data against each other
plot_density_all_greyscale <- function(out, labels, xlab, b, ...)
{
        # Get first set of data
        args <- list(...)
        a = args[[1]]
        l = density(a, bw=b)
        all = cbind(l$x, l$y)
        m = mean(a, na.rm=T)

        # Get remaining data
        args = args[2:length(args)]
        for (a in args){
                l = density(a, bw=b)
                m = cbind(m, mean(a, na.rm=T))
                all = cbind(all, l$x, l$y)
        }

        # Setup plot variables
        odd = (1:(dim(all)[2]/2))*2 - 1
        even = (1:(dim(all)[2]/2))*2
        ymax=max( all[,even])
        ymin=min( all[,even])
        xmin=min( all[,odd])
        xmax=max( all[,odd])

        greys = colors()[grep("grey", colors())]
        greys =  greys[grep ("[0-9]", greys, invert=T) ]
        cols=c("black", greys )

        legcols = list()

        # Plot
        png( paste(out, "png", sep="."))
        plot(-xmax,-ymax, type="l", lwd=2, col="grey", ylim=c(ymin,ymax),xlim=c(xmin,xmax), ylab="Density", xlab=xlab)
        j = 1
	for (i in odd){
                lines( all[,c(i,i+1)], lwd=3, col=cols[j])
		legcols = append(legcols, cols[j])
                abline( v= m[j], lty=2, col=cols[j])
                text( m[j], max(all[,c(i+1)])+ymax*0.01, label=round(m[j],2))
                j = j + 1
        }
        print(legcols)
        legend("topleft", (labels), col=unlist(legcols), bty="n", lwd=3)
        dev.off()
}

plot_hist_all_greyscale <- function(out, labels, xlab, b, ...)
{
        # Get first set of data
        args <- list(...)
        a = args[[1]]
        hist = hist(a, breaks=b, plot=F)
        m = mean(a, na.rm=T)
        all =  get_density(hist)
        b = hist$breaks

        # Get remaining data
        args = args[2:length(args)]
        for (a in args){
                hist = hist(as.numeric(a), breaks=b, plot=F)
                m = cbind(m, mean(a, na.rm=T))
                temp = get_density(hist)
                all = cbind(all, temp)
        }

        # Setup plot variables
        odd = (1:(dim(all)[2]/2))*2 - 1
        even = (1:(dim(all)[2]/2))*2
        ymax=max( all[,even])
        ymin=min( all[,even])
        xmin=min( all[,odd])
        xmax=max( all[,odd])

        #n_cols = dim(all)[2]/2
        #greys = colors()[grep("grey", colors())]
        #cols=c("black", greys[sort(sample(length(greys), n_cols))] )

        greys = colors()[grep("grey", colors())]
        greys =  greys[grep ("[0-9]", greys, invert=T) ]
        cols=c("black", greys )

        legcols = list()
        newlabels = list()

        # Plot
        png( paste(out, "png", sep="."))
        plot(-xmax,-ymax, type="l", lwd=2, col="grey", ylim=c(ymin,ymax),xlim=c(xmin,xmax), ylab="Density", xlab=xlab )
        j = 1
	for (i in odd){
                lines( all[,c(i,i+1)], lwd=2, col=cols[j])
		legcols = append(legcols, cols[j])
                newlabels = append(newlabels, paste( labels[j], "avg=", round(m[j],2) ))
                j = j + 1
        }
        legend("topleft", unlist(newlabels), col=unlist(legcols), bty="n", lwd=2)
        dev.off()
}


plot_convolved_greyscale <- function(out, f, listA, listB, xlab, ylab){
        conv_listB = convolve_nets( listA, listB, f)
        png(paste(out, "png",sep="."))
        plot( listA, listB, pch=19, cex=0.25, col="grey", ylab=ylab, xlab=xlab, sub=round(cor(listA, listB, method="s"),3))
        lines( conv_listB, lwd=3)
        dev.off()
}









# Tic toc functions
tic <- function(gcFirst = TRUE, type=c("elapsed", "user.self", "sys.self"))
{
	type <- match.arg(type)
	assign(".type", type, envir=baseenv())
	if(gcFirst) gc(FALSE)
	tic <- proc.time()[type]
	assign(".tic", tic, envir=baseenv())
	invisible(tic)
}

toc <- function()
{
	type <- get(".type", envir=baseenv())
	toc <- proc.time()[type]
	tic <- get(".tic", envir=baseenv())
	print(toc - tic)
	invisible(toc)
}



# Transparent colors
makeTransparent<-function(someColor, alpha=100)
{
	newColor<-col2rgb(someColor)
	apply(newColor, 2, function(curcoldata){rgb(red=curcoldata[1], green=curcoldata[2],
	blue=curcoldata[3],alpha=alpha, maxColorValue=255)})
}

