#########################################
#    Guilt-by-Association: how to       #
#########################################

##  Written: Sara Ballouz
##  Date: January 21st, 2014

# Load data
load("data/GO_semantic_sim.Rdata")
load("data/KEGG.Rdata")
load("data/lig_data.Rdata")
load("data/extPPIN.Rdata")
load("data/ppin.Rdata")
load("data/coexp.Rdata")
load("data/9606.gene_attr.Rdata")

# Load functions
source("code/GBA_helper.r")
source("code/GBA.r")

#########################################
# Data descriptions
## GO_semantic_sim.Rdata
### GO.filt.prop : GO annotations without IEA, propogated
### GO.filt      : GO annotations without IEA, not propogated
### GO           : All GO annotations
### GO.IEA.prop  : Only IEA GO annotations, propogated
### GO.IEA       : ONLY IEA GO annotations, not propogated
### semantic_similarity_GO.filt.prop : Semantic similarity based on GO.filt.prop
### semantic_similarity_GO.filt      : Semantic similarity based on GO.filt
### semantic_similarity_GO           : Semantic similarity based on GO
### semantic_similarity_GO.IEA.prop  : Semantic similarity based on GO.IEA.prop
### semantic_similarity_GO.IEA       : Semantic similarity based on GO.IEA

## KEGG.Rdata
### KEGG : KEGG annotations, not filtered

## lig_data.Rdata
### all_genes : genes
### data.mat  : ligand data

## extPPIN.Rdata
### fullextPPI: full extended PPI network

## ppin.Rdata
### fullPPI: full PPI network from biogrid

## coexp.Rdata
### fullcoexpNetThres: full coexpression network, scores not ranked, thresholded for top 0.5% connections

#########################################

# To run on ligand data:
## ligand network
fullligNet = data.mat
## GO annotations, filtered for GO terms between minGOsize and maxGOsize genes
minGOsize = 20
maxGOsize = 1000

GO.filtered = filter_network(GO.filt.prop, minGOsize, maxGOsize000)

# Filter network and annotations
## make sure to check the labels match, should be sorted according to genes
genes =  rownames(fullligNet)
genes.labels = get_subnet_rows( GO.filtered, genes)

# Filter network symmetrically
network = get_subnet(fullligNet, rownames(genes.labels), rownames(genes.labels) )

# OR filter by rows (leaving columns, for extra predictions of unannotated genes)
network = get_subnet_rows(fullligNet, rownames(genes.labels) )

# Run GBA cross validation to get AUROC scores
## Number of folds
nFold = 3
scores = neighbor.voting.CV(genes.labels, network, nFold)
## every row has the AUROC for the GO term in the first column, second column is the average node degree of the term, and the last column is the node degree AUROC.

averageAUROC = mean(scores[,1], na.rm=T)

# Run GBA to get scores for gene predictions
predictions =  predictions(genes.labels, network)
## every column has the GO term, and the rows are the score for each gene
## rerank accordingly


# To run on KEGG annotations, eg:
KEGG.filtered = filter_network(KEGG, 2,10000)
genes =  rownames(fullligNet)
genes.labels = get_subnet_rows( KEGG.filtered, genes)
network = get_subnet_rows(fullligNet, rownames(genes.labels) )
scores = neighbor.voting.CV(genes.labels, network, 3)
averageAUROC = mean(scores[,1], na.rm=T)
predictions =  predictions(genes.labels, network)


# Word of caution:
# For dense and large networks, this takes a little time ~ but should be on the order of 
# a few seconds for the smaller ligand networks
