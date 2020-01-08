#########################################
#    Guilt-by-association functions     #
#########################################

##  Written: Sara Ballouz
##  Date: January 20th, 2014
##  Updated: June 18th, 2015

# Given gene annotations (genes.labels) and matching network,
# returns list of predictions (scored) for each annotation using
# simple neighbor voting

predictions <- function(genes.labels, network)
{
        # genes.label : needs to be in 1s and 0s
        l <- dim(genes.labels)[2]
        g <- dim(genes.labels)[1]
        ab <- which(genes.labels != 0, arr.ind=T)
        n  <- length(ab[,1])

	# max
        diag(network) <- max(network, na.rm=T)

        # Get sums - mat. mul.
        sumin    = ( t(network) %*% genes.labels)

        # Sum of all edge in network
        sumall   = matrix(apply(network,2,sum), ncol = dim(sumin)[2], nrow=dim(sumin)[1])

        # Predictions
        predicts = sumin/sumall
        return(predicts)
}

# Given gene annotations (genes.labels) and matching network,
# returns AUROC for each annotation using
# simple neighbor voting

neighbor.voting.CV <- function(genes.labels,network,nFold){

        # genes.label : needs to be in 1s and 0s
        l <- dim(genes.labels)[2]
        g <- dim(genes.labels)[1]
        ab <- which(genes.labels != 0, arr.ind=T)
        n  <- length(ab[,1])

	# max
        diag(network) <- max(network, na.rm=T)

        #print("Make genes label CV matrix")
        test.genes.labels = matrix(genes.labels, nrow=g, ncol=nFold*l)

        # For each fold in each GO group, remove 1/nth of the values of the genes.label
        for (j in 1:l){
                d <- which(ab[,2]==j)  # Which indices the genes are in this particular GO group
                t <- length(d)         # Total number of genes in the GO group
                r <- sample(1:t, replace=F)
                f <- t/nFold
                for (i in 1:nFold){
                        e = c( (1:f) + f*(i-1) )
                        e = sort(r[e])
                        c = j + l*(i-1)        # GO group to look at (ie column)
                        test.genes.labels[ab[d],c][e] <- 0
                }
        }


        #print("Get sums - mat. mul.")
        sumin    = ( t(network) %*% test.genes.labels)

        #print("Get sums - calc sumall")
        sumall   = matrix(apply(network,2,sum), ncol = dim(sumin)[2], nrow=dim(sumin)[1])

        #print("Get sums - calc predicts")
        predicts = sumin/sumall

        predicts0 = predicts
        predicts = get_subnet_rows(predicts, rownames(genes.labels))

        #print("Hide training data")
        nans = which(test.genes.labels == 1, arr.ind=T)
        predicts[nans] <- NA

        #print("Rank test data")
        predicts = apply(abs(predicts), 2, rank,na.last="keep",ties.method="average")


        filter = matrix(genes.labels, nrow=g, ncol=nFold*l)
        negatives = which(filter == 0, arr.ind=T)
        positives = which(filter == 1, arr.ind=T)

        predicts[negatives] <- 0

        #print("Calculate ROC - np")
        np = apply(filter,2,sum) - apply(test.genes.labels,2,sum) # Postives

        #print("Calculate ROC - nn")
        nn = dim(test.genes.labels)[1] - apply(filter,2,sum)      # Negatives

        #print("Calculate ROC - p")
        p =  apply(predicts,2,sum,na.rm=T)

        #print("Calculate ROC - rocN")
        rocN = (p/np - (np+1)/2)/nn
        rocN = matrix(rocN, ncol=nFold, nrow=l)
        rocN = rowMeans(rocN)

        #print("Calculate node degree")
        node_degree <- apply(network, 1,sum)
        colsums <- apply(genes.labels, 2, sum)

        #print("Calculate node degree - sum across gene labels")
        node_degree <- matrix(node_degree)
        temp <- t(node_degree) %*% genes.labels


        #print("Calculate node degree - average")
        average_node_degree <- t(temp)/colsums

        #print("Calculate node degree roc - rank node degree")
        ranks <- apply(abs(node_degree), 2, rank,na.last="keep",ties.method="average")
        ranks <- matrix(ranks, nrow=length(ranks), ncol=dim(genes.labels)[2])

        #print("Calculate node degree roc - remove negatives")
        negatives = which(genes.labels == 0, arr.ind=T)
        ranks[negatives] = 0

        #print("Calculate node degree roc - np")
        np = apply(genes.labels,2,sum)

        #print("Calculate node degree roc - nn")
        nn = dim(genes.labels)[1] - np

        #print("Calculate node degree roc - p")
        p =  apply(ranks,2,sum,na.rm=T)

        #print("Calculate node degree roc - roc")
        roc  <-  (p/np - (np+1)/2)/nn

        scores = cbind(rocN,matrix(average_node_degree)[,1],roc)

}

