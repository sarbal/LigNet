#########################################
#  Network analysis helper functions    #
#########################################

##  Written: Sara Ballouz
##  Date: January 19th, 2014


## Necessary libraries
library(MASS)
library(Matrix)
library(scales)
library(gplots)

# Network manipulation

## Returns network filtered on row labels
get_subnet_rows <- function(network, row_labels)
{
        m = match(row_labels, rownames(network))
        f1 = !is.na(m)
        f2 = m[f1]

        subnet = network[f2,]
        return(subnet)
}


## Returns network filtered on column labels
get_subnet_cols <- function(network, col_labels)
{

        m = match(col_labels, colnames(network))
        f1 = !is.na(m)
        f2 = m[f1]

        subnet = network[,f2]
        return(subnet)
}


## Returns network filtered on both given row and column labels
get_subnet <- function(network, row_labels, col_labels)
{
        m = match(row_labels, rownames(network))
	f1 = !is.na(m)
	f2 = m[f1]

        m = match(col_labels, colnames(network))
	f3 = !is.na(m)
	f4 = m[f3]

        subnet = network[f2,f4]
        return(subnet)
}

## Returns network filtered on both given row and column labels
filter_network <- function(network, min, max)
{
	# Filter terms between min and max
	colsums <- apply(network, 2, sum)
	col.filter <- which( (colsums > min & colsums < max) )

        # Filter for genes with terms
	rowsums <- apply(network, 1, sum)
	row.filter <- which( rowsums > 0 )
	network <- network[row.filter,col.filter]

        return(network)
}


## Returns network filtered on both given row and column labels
filter_network_cols <- function(network, min, max)
{
	# Filter terms between min and max
	colsums <- apply(network, 2, sum)
	col.filter <- which( (colsums > min & colsums < max) )

	network <- network[,col.filter]

        return(network)
}

## Formats the density distribution from the histogram function
get_density <- function(hist)
{
        x = sort(rep(hist$breaks,2))
        y = matrix(rbind(hist$density, hist$density))
        y = c(0,y,0)

        return(cbind(x,y))
}


## Formats the counts distribution from the histogram function
get_counts <- function(hist)
{
        x = sort(rep(hist$breaks,2))
        y = matrix(rbind(hist$counts, hist$counts))
        y = c(0,y,0)

        return(cbind(x,y))
}


## Thresholds the network to top p% of genes
threshold_network_top_genes <- function(network, p)
{

        # Setup network
        diag(network) <- 0
        n = dim(network)[1]
        x = n-1
        top = ceiling(x*p/100)

        # Rank network, remove edges below top ranked
        ranks <-  apply( network, 2, rank,na.last="keep",ties.method="first")
        x = which(ranks < (n-top))
        network[x]=0

        # Make symmetrical
        a = network != 0
        b =  t(network) != 0
        c = a+b
        d = c==0
        c[d] = 1
	network = (network + t(network))/c

        return(network)
}


## Calculates the ROC for a given list of scores(or ranks) and the known labels
roc_score <- function(scores,labels)
{
	negatives = which(labels == 0, arr.ind=T)
	scores[negatives] <- 0

	p  = sum(scores)
	nL = length(labels)
	np = sum(labels)
        nn = nL - np

	roc  <-  (p/np - (np+1)/2)/nn

	return(roc)
}
# Calculates the ROC for a given list of scores(or ranks) and the known labels
plot_roc <- function(scores,labels, file)
{
        h1 = scores * labels
        h2 = scores * !labels

	auroc  <-  roc_score(scores, labels)
	tpr = cumsum(h1)/sum(h1)
	fpr = cumsum(h2)/sum(h2)
	png(file)
	plot(fpr,tpr, xlab="FPR", ylab="TPR")
	lines( c(0,1), c(0,1), col="grey")
	text(roc, 0.5, 1 )
	dev.off()
        pval = wilcox.test(scores[labels], scores[!labels])
	return(cbind(fpr,tpr))
}

plot_roc <- function(scores,labels, file)
{
        h1 = scores * labels
        h2 = scores * !labels

	auroc  <-  roc_score(scores, labels)
	tpr = cumsum(h1)/sum(h1)
	fpr = cumsum(h2)/sum(h2)
	plot(fpr,tpr, xlab="FPR", ylab="TPR")
	lines( c(0,1), c(0,1), col="grey")
	text(labels=auroc, 0.5, 1 )
        pval = wilcox.test(scores[labels], scores[!labels])
	return(cbind(fpr,tpr))
}




# Plotting functions

## Plots the heatmap of a given correlation matrix. Warning, for large matrices, this is slow/unusable
plot_heatmap <-function(out, corr)
{
	png(paste( out, "heatmap.png", sep="."))
	heatmap.2(corr, trace="none", density.info="none" )
	dev.off()
}


## Plots list of points versus another
plot_comp <-function(out, list1, list2, xlab,ylab)
{
	png(paste( out, "vs.png", sep="."))
	plot(list1, list2, xlab=xlab, ylab=ylab, pch=19, cex=0.5)
	lines(stats::lowess(cbind(list1,list2)),col="purple")
	lines( c(0,1),c(0,1),col="grey")
	dev.off()
}

## Plots the distribution (densities) of different data against each other
plot_hist_all <- function(out, labels, xlab, b, ...)
{
        # Get first set of data
        args <- list(...)
        a = args[[1]]
        hist = hist(a, breaks=b, plot=F)
        all =  get_density(hist)
        b = hist$breaks
        
        # Get remaining data
        args = args[2:length(args)]
        for (a in args){
                hist = hist(as.numeric(a), breaks=b, plot=F)
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
        cols=rainbow(dim(all)[2])

        # Plot
        png( paste(out, "png", sep="."))
        plot(-xmax,-ymax, type="l", lwd=2, col="grey", ylim=c(ymin,ymax),xlim=c(xmin,xmax), ylab="Density", xlab=xlab )
        for (i in odd){
                lines( all[,c(i,i+1)], lwd=2, col=cols[i] )
        }
        legend("topleft", labels, col= cols[odd], bty="n", lwd=2)
        dev.off()
}


## Plots the distribution (densities) of different data against each other
plot_hist_all_greyscale <- function(out, labels, xlab, b, ...)
{
        # Get first set of data
        args <- list(...)
        a = args[[1]]
        hist = hist(a, breaks=b, plot=F)
        all =  get_density(hist)
        b = hist$breaks
        
        # Get remaining data
        args = args[2:length(args)]
        for (a in args){
                hist = hist(as.numeric(a), breaks=b, plot=F)
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
        #cols=rainbow(dim(all)[2])
        cols=c("grey","black")
        leglty = list()
	legcols = list()

        # Plot
        png( paste(out, "png", sep="."))
        plot(-xmax,-ymax, type="l", lwd=2, col="grey", ylim=c(ymin,ymax),xlim=c(xmin,xmax), ylab="Density", xlab=xlab )
        j = 1
	for (i in odd){
                lines( all[,c(i,i+1)], lwd=2, col=cols[(j%%2) + 1 ], lty=j )   
		legcols = append(legcols, cols[(j%%2) + 1 ])
		leglty = append(leglty, j)

                j = j + 1
        }
        legend("topleft", labels, col=unlist(legcols) , lty=unlist(leglty), bty="n", lwd=2)
        dev.off()
}


## Plots the distribution (densities) of different data against each other
plot_hist_all_greyscale <- function(out, labels, xlab, b, ...)
{
        # Get first set of data
        args <- list(...)
        a = args[[1]]
        hist = hist(a, breaks=b, plot=F)
        all =  get_density(hist)
        b = hist$breaks

        # Get remaining data
        args = args[2:length(args)]
        for (a in args){
                hist = hist(as.numeric(a), breaks=b, plot=F)
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
        #cols=rainbow(dim(all)[2])
        cols=rev(c("black", "darkslategrey", "darkgrey"))
        legcols = list()

        # Plot
        png( paste(out, "png", sep="."))
        plot(-xmax,-ymax, type="l", lwd=2, col="grey", ylim=c(ymin,ymax),xlim=c(xmin,xmax), ylab="Density", xlab=xlab )
        j = 1
	for (i in rev(odd)){
                lines( all[,c(i,i+1)], lwd=2, col=cols[j])
		legcols = append(legcols, cols[j])

                j = j + 1
        }
        legend("topleft", (labels), col=unlist(rev(legcols)), bty="n", lwd=2)
        dev.off()
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
        #cols=rainbow(dim(all)[2])
        cols=(c("black", "darkslategrey", "darkgrey"))
        legcols = list()

        # Plot
        png( paste(out, "png", sep="."))
        plot(-xmax,-ymax, type="l", lwd=2, col="grey", ylim=c(ymin,ymax),xlim=c(xmin,xmax), ylab="Density", xlab=xlab)
        j = 1
	for (i in odd){
                lines( all[,c(i,i+1)], lwd=3, col=cols[j])
		legcols = append(legcols, cols[j])
                abline( v= m[j], lty=2, col=cols[j])
                text( m[j], max(all[,c(i+1)])+ymax*0.01, label=round(m[j],3))
                j = j + 1
        }
        legend("topleft", (labels), col=unlist(legcols), bty="n", lwd=3)
        dev.off()
}

# Gene ID mapping functions

## Returns gene properties from given gene list
### requires gene_attr
get_gene_size <- function(genes)
{

        m = match( genes, gene_attr$entrez_gene_id)
        f = !is.na(m)
        z = m[f]

        return( cbind(genes[f],gene_attr$gene_length[z]))
}


## Returns gene Entrez IDs given gene symbol
get_gene_id_from_name <- function(genes)
{
        m = match( genes, gene_attr$gene_symbol)
        f = !is.na(m)
        z = m[f]

        return( cbind(genes[f],gene_attr$entrez_gene_id[z]))
}


## Returns gene symbol given gene Entrez IDs
get_gene_name_from_id <- function(genes)
{
        m = match( genes, gene_attr$entrez_gene_id)
        f = !is.na(m)
        z = m[f]

        return( cbind(genes[f],gene_attr$gene_symbol[z] ))
}

