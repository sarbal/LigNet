data = read.table("/home/sballouz/data/SEA/9606.scores_best.parsed.tsv", header=T)
biogrid =  read.table("/home/sballouz/data/bioGrid/BIOGRID-ALL-3.2.106/9606.biogrid.parsed.tsv", sep="\t", header=T)
labels =  read.table("/home/sballouz/data/negatome/v2/labels", sep="\t")
negatome = read.table("/home/sballouz/data/negatome/v2/9606.negatome.parsed.tsv", header=T)

all_genes = sort(unique( c(data[,1], data[,3], biogrid[,1], biogrid[,2], negatome[,1],negatome[,2])))

all_genes = sort(unique( c(data[,1], data[,3])))
data.mat = matrix(0, ncol=length(all_genes), nrow=length(all_genes))
rownames(data.mat) = all_genes
colnames(data.mat) = all_genes

bio.mat = data.mat
neg.mat = data.mat


m.a = match(negatome[,1],all_genes)
f.na = !is.na(m.a)
f.na2 = m.a[f.na]

m.b = match(negatome[,2], all_genes)
f.nb = !is.na(m.b)
f.nb2 = m.b[f.nb]

neg.mat[cbind(f.na2, f.nb2)] = 1

m.a = match(biogrid[,1],all_genes)
f.na = !is.na(m.a)
f.na2 = m.a[f.na]

m.b = match(biogrid[,2], all_genes)
f.nb = !is.na(m.b)
f.nb2 = m.b[f.nb]

bio.mat[cbind(f.na, f.nb)] = 1


m.a = match(data[,1],all_genes)
f.na = !is.na(m.a)
f.na2 = m.a[f.na]

m.b = match(data[,3], all_genes)
f.nb = !is.na(m.b)
f.nb2 = m.b[f.nb]

data.mat[cbind(f.na2, f.nb2)] = 1

d = dim(data.mat)
A =(rep(all_genes, each=d[2]))
B = rep(rep(all_genes), d[2])
C = cbind(A,B,array(data.mat),array(bio.mat),array(neg.mat))



filt = C[,3] == 1

labelsfile =  "~/data/GO_March/final/goa_human.split.parsed.filt.desc"
GO = get_labels(labelsfile,"GOID")



labelsfile = "/KEGG/uid2KEGG.all"
KEGG = get_labels(labelsfile,"KEGG")
labelsfile = "~/sballouz/human/reactome/human.reactome.parsed"
reactome = get_labels(labelsfile,"reactome")


load("~/data/mat_files/human_RNAseq_aggregate_ranked.RData")

genes.coexp = rownames(network)
coexp.mat = matrix(0, ncol=length(all_genes), nrow=length(all_genes))

m = match(all_genes,genes.coexp)
f.n = !is.na(m)
f.n2 = m[f.n]

coexp.mat[f.n, f.n] = 1
subnet = network[f.n2,f.n2]
n = dim(subnet)[1]
ranked = rankM(subnet, n)

diag(subnet) = 0
ranked2 = rankM(subnet,n)


temp = array(data.mat[f.n, f.n])
png("test.png")
conv_smoother( array(subnet), temp, 1000, "coexp", "ligand")
dev.off()

png("test2.png")
conv_smoother( array(ranked), temp, 1000, "coexp", "ligand")
dev.off()

png("test10000.png")
conv_smoother( array(subnet), temp, 10000, "coexp", "ligand")
dev.off()


png("test10000-2.png")
conv_smoother( array(ranked), temp, 10000, "coexp", "ligand")
dev.off()


png("test10000-3.png")
conv_smoother( array(ranked2), temp, 10000, "coexp", "ligand")
dev.off()


png("test10000-4.png")
conv_smoother4( array(ranked2), temp, 10000, "coexp", "ligand")
dev.off()

ranked3 = ranked2
diag(ranked3) = NA

png("test10000-5.png")
conv_smoother4( array(ranked3), temp, 10000, "coexp", "ligand")
dev.off()

png("test5000-5.png")
conv_smoother4( array(ranked3), temp, 5000, "coexp", "ligand")
dev.off()

png("test1000-5.png")
conv_smoother4( array(ranked3), temp, 1000, "coexp", "ligand")
dev.off()

repmat <- function(X,m,n){
	##R equivalent of repmat (matlab)
	mx = dim(X)[1]
	nx = dim(X)[2]
	return (matrix(t(matrix(X,mx,nx*n)),mx*m,nx*n,byrow=T))
}



load("~/data/mat_files/human_RNAseq_aggregate.RData")











all_genes2 = sort( unique( c(biogrid[,1], biogrid[,2])))
bio.mat = matrix(0, ncol=length(all_genes2), nrow=length(all_genes2))
rownames(bio.mat) = all_genes2
colnames(bio.mat) = all_genes2

m.a = match(biogrid[,1],all_genes2)
f.na = !is.na(m.a)
f.na2 = m.a[f.na]

m.b = match(biogrid[,2], all_genes2)
f.nb = !is.na(m.b)
f.nb2 = m.b[f.nb]

cbind(f.na2,f.nb2)

 print_mtx(bio.mat, "all", "biogrid")




J = get_network("GO.jaccard.mtx.gz","GO.jaccard.keys")


m.j = match(rownames(J), all_genes)
f.nj = !is.na(m.j)
f.nj2 = m.j[f.nj]

subJ = J[f.nj,f.nj]
 diag(subJ) = 0
 temp = rankM(subJ,dim(subJ)[1])
 diag(temp) = NA
 rankedJ = array(temp)
ligJ = array(data.mat[f.nj2, f.nj2])

png("semsim10000.png")
conv_smoother4( rankedJ, ligJ, 10000, "semantic sim. GO", "ligand")
dev.off()



### sem version 2
GO = get_labels("~/data/GO_March/final/goa_human.split.parsed.filt.desc", "GO")
num = GO %*% t(GO)
sum = apply(GO, 1, sum)
sum_rep = repmat( as.matrix(sum), 1, length(sum))
denom = t(sum_rep) + sum_rep - num
J = num/denom

t = which(is.na(J))
J[t] = 0

m.j = match(rownames(J), all_genes)
f.nj = !is.na(m.j)
f.nj2 = m.j[f.nj]

subJ = J[f.nj,f.nj]
 diag(subJ) = 0
 temp = rankM(subJ,dim(subJ)[1])
 diag(temp) = NA
 rankedJ = array(temp)
ligJ = array(data.mat[f.nj2, f.nj2])

png("semsim10000-2.png")
conv_smoother4( rankedJ, ligJ, 10000, "semantic sim. GO", "ligand")
dev.off()


m.j = match(rownames(J), all_genes[f.n])
f.nj = !is.na(m.j)
f.nj2 = m.j[f.nj]

subJ = J[f.nj,f.nj]
 diag(subJ) = 0
 temp = rankM(subJ,dim(subJ)[1])
 diag(temp) = NA
 rankedJ = array(temp)
ligJ = array(data.mat[f.nj2, f.nj2])

png("semsim10000-3.png")
conv_smoother4( rankedJ, ligJ, 10000, "semantic sim. GO", "ligand")
dev.off()




m.g = match (all_genes, unlist(gene_list))
f.ng = !is.na(m.g)
f.ng2 = m.g[f.ng]


m.j = match(rownames(J), all_genes[f.ng] )
f.nj = !is.na(m.j)
f.nj2 = m.j[f.nj]

subJ = J[f.nj,f.nj]
 diag(subJ) = 0
 temp = rankM(subJ,dim(subJ)[1])
 diag(temp) = NA
 rankedJ = array(temp)
ligJ = array(data.mat[f.ng, f.ng])

diag(data.mat) = 0

png("semsim10000-4.png")
conv_smoother4( rankedJ, ligJ, 10000, "semantic sim. GO", "ligand")
dev.off()


### sem version 3
GOsum = apply(GO, 2, sum)
filt = which(GOsum > 10 & GOsum < 300) 

num = GO[,filt] %*% t(GO[,filt])
sum = apply(GO[,filt], 1, sum)
sum_rep = repmat( as.matrix(sum), 1, length(sum))
denom = t(sum_rep) + sum_rep - num
J = num/denom

t = which(is.na(J))
J[t] = 0

m.j = match(rownames(J), all_genes)
f.nj = !is.na(m.j)
f.nj2 = m.j[f.nj]

subJ = J[f.nj,f.nj]
 diag(subJ) = 0
 temp = rankM(subJ,dim(subJ)[1])
 diag(temp) = NA
 rankedJ = array(temp)
ligJ = array(data.mat[f.nj2, f.nj2])

png("semsim10000-2a.png")
conv_smoother4( rankedJ, ligJ, 10000, "semantic sim. GO", "ligand")
dev.off()


m.j = match(rownames(J), all_genes[f.n])
f.nj = !is.na(m.j)
f.nj2 = m.j[f.nj]

subJ = J[f.nj,f.nj]
 diag(subJ) = 0
 temp = rankM(subJ,dim(subJ)[1])
 diag(temp) = NA
 rankedJ = array(temp)
ligJ = array(data.mat[f.nj2, f.nj2])

png("semsim10000-3a.png")
conv_smoother4( rankedJ, ligJ, 10000, "semantic sim. GO", "ligand")
dev.off()



m.g = match (all_genes, unlist(gene_list))
f.ng = !is.na(m.g)
f.ng2 = m.g[f.ng]

m.j = match(rownames(J), all_genes[f.ng] )
f.nj = !is.na(m.j)
f.nj2 = m.j[f.nj]

subJ = J[f.nj,f.nj]
diag(subJ) = 0
temp = rankM(subJ,dim(subJ)[1])
diag(temp) = NA
rankedJ = array(temp)
ligJ = array(data.mat[f.ng, f.ng])

png("semsim10000-4a.png")
conv_smoother4( rankedJ, ligJ, 10000, "semantic sim. GO", "ligand")
dev.off()



### sem version 4
GOsum = apply(GO, 2, sum)
filt = which(GOsum > 20 & GOsum < 1000)


m.go = match(rownames(GO[,filt]), all_genes)
f.go = !is.na(m.go)
f.go2 = m.go[f.go]

num = GO[f.go,filt] %*% t(GO[f.go,filt])
sum = apply(GO[f.go,filt], 1, sum)
sum_rep = repmat( as.matrix(sum), 1, length(sum))
denom = t(sum_rep) + sum_rep - num
J = num/denom

t = which(is.na(J))
J[t] = 0


m.j = match(rownames(J), all_genes)
f.nj = !is.na(m.j)
f.nj2 = m.j[f.nj]

subJ = J[f.nj,f.nj]
 diag(subJ) = 0
 temp = rankM(subJ,dim(subJ)[1])
 diag(temp) = NA
 rankedJ = array(temp)
ligJ = array(data.mat[f.nj2, f.nj2])

png("semsim10000-2b.png")
conv_smoother4( rankedJ, ligJ, 10000, "semantic sim. GO", "ligand")
dev.off()


m.j = match(rownames(J), all_genes[f.n])
f.nj = !is.na(m.j)
f.nj2 = m.j[f.nj]

subJ = J[f.nj,f.nj]
 diag(subJ) = 0
 temp = rankM(subJ,dim(subJ)[1])
 diag(temp) = NA
 rankedJ = array(temp)
ligJ = array(data.mat[f.nj2, f.nj2])

png("semsim10000-3b.png")
conv_smoother4( rankedJ, ligJ, 10000, "semantic sim. GO", "ligand")
dev.off()



m.g = match (all_genes, unlist(gene_list))
f.ng = !is.na(m.g)
f.ng2 = m.g[f.ng]


m.j = match(rownames(J), all_genes[f.ng] )
f.nj = !is.na(m.j)
f.nj2 = m.j[f.nj]

subJ = J[f.nj,f.nj]
 diag(subJ) = 0
 temp = rankM(subJ,dim(subJ)[1])
 diag(temp) = NA
 rankedJ = array(temp)
ligJ = array(data.mat[f.ng, f.ng])

diag(data.mat) = 0

png("semsim10000-4b.png")
conv_smoother4( rankedJ, ligJ, 10000, "semantic sim. GO", "ligand")
dev.off()


### sem version 5
filt = which(GOsum > 10 & GOsum < 300)
m.go = match(rownames(GO[,filt]), all_genes)
f.go = !is.na(m.go)
f.go2 = m.go[f.go]

num = GO[f.go,filt] %*% t(GO[f.go,filt])
sum = apply(GO[f.go,filt], 1, sum)
sum_rep = repmat( as.matrix(sum), 1, length(sum))
denom = t(sum_rep) + sum_rep - num
J = num/denom

t = which(is.na(J))
J[t] = 0

m.j = match(rownames(J), all_genes)
f.nj = !is.na(m.j)
f.nj2 = m.j[f.nj]

subJ = J[f.nj,f.nj]
diag(subJ) = 0
temp = rankM(subJ,dim(subJ)[1])
diag(temp) = NA
rankedJ = array(temp)
ligJ = array(data.mat[f.nj2, f.nj2])

png("semsim10000-2c.png")
conv_smoother4( rankedJ, ligJ, 10000, "semantic sim. GO", "ligand")
dev.off()


m.j = match(rownames(J), all_genes[f.n])
f.nj = !is.na(m.j)
f.nj2 = m.j[f.nj]

subJ = J[f.nj,f.nj]
 diag(subJ) = 0
 temp = rankM(subJ,dim(subJ)[1])
 diag(temp) = NA
 rankedJ = array(temp)
ligJ = array(data.mat[f.nj2, f.nj2])

png("semsim10000-3b.png")
conv_smoother4( rankedJ, ligJ, 10000, "semantic sim. GO", "ligand")
dev.off()



m.g = match (all_genes, unlist(gene_list))
f.ng = !is.na(m.g)
f.ng2 = m.g[f.ng]


m.j = match(rownames(J), all_genes[f.ng] )
f.nj = !is.na(m.j)
f.nj2 = m.j[f.nj]

subJ = J[f.nj,f.nj]
 diag(subJ) = 0
 temp = rankM(subJ,dim(subJ)[1])
 diag(temp) = NA
 rankedJ = array(temp)
ligJ = array(data.mat[f.ng, f.ng])

diag(data.mat) = 0

png("semsim10000-4b.png")
conv_smoother4( rankedJ, ligJ, 10000, "semantic sim. GO", "ligand")
dev.off()







m.g = match (all_genes, unlist(gene_list))
f.ng = !is.na(m.g)
f.ng2 = m.g[f.ng]


m.j = match(rownames(J), all_genes[f.ng] )
f.nj = !is.na(m.j)
f.nj2 = m.j[f.nj]

subJ = J[f.nj,f.nj]
 diag(subJ) = 0
 temp = rankM(subJ,dim(subJ)[1])
 diag(temp) = NA
 rankedJ = array(temp)
ligJ = array(data.mat[f.ng, f.ng])

diag(data.mat) = 0

png("semsim10000-4b.png")
conv_smoother4( rankedJ, ligJ, 10000, "semantic sim. GO", "ligand")
dev.off()














### PPIN+ext
bio.mat.weights =  get_network("~/data/SEA/biogrid.all-paths-weights.mtx", "~/data/SEA/biogrid.all.keys")

m.ppi = match( rownames(bio.mat.weights),all_genes)
f.ppi = !is.na(m.ppi)
f.ppi2 = m.ppi[f.ppi]

subPPI = bio.mat.weights[f.ppi,f.ppi]
ligPPI = array(data.mat[f.ppi2, f.ppi2])
diag(subPPI) = 0
 temp = rankM(subPPI,dim(subPPI)[1])
 diag(temp) = NA
 rankedPPI = array(temp)

png("ppi10000.png")
conv_smoother4( rankedPPI, ligPPI, 10000, "PPIN (biogrid) + ext", "ligand")
dev.off()

 temp = matrix(rank(subPPI, na.last="keep", ties.method="first"), nrow=dim(subPPI)[1]))
 diag(temp) = NA
 rankedPPI = array(temp)


png("ppi10000b.png")
conv_smoother4( rankedPPI, ligPPI, 10000, "PPIN (biogrid) + ext", "ligand")
dev.off()


 temp = matrix(rank(subPPI, na.last="keep", ties.method="random"), nrow=dim(subPPI)[1])
 diag(temp) = NA
 rankedPPI = array(temp)

png("ppi10000c.png")
conv_smoother4( rankedPPI, ligPPI, 10000, "PPIN (biogrid) + ext", "ligand")
dev.off()


	n <- order(rankedPPI)
	m <- rankedPPI[n]
	i <- length(m)/f
	temp_b1 = m[ f/2 : (length(m)-f/2)]
	temp_b2 = convolve( ligPPI[n], rep(1,f),type="filter")
	temp_b2 = temp_b2/f




convolved = cbind(temp_b1,temp_b2)
colnames(convolved) = c("ppi","lig")
png("boxplot.png")
boxplot(lig~ppi, data=convolved,  xlab="PPIN (biogrid) + ext", ylab="ligand")
dev.off()



png("ppi10000.png")
conv_smoother_boxplot( rankedPPI, ligPPI, 10000, "PPIN (biogrid) + ext", "ligand")
dev.off()

ylab = "Ligand incidence"
ppilab = "PPIN (biogrid + ext) rank"
semlab = "Semantic similarity (GO) rank"
coexplab = "Coexpression rank"

genes.coexp = rownames(network)
genes.bio = rownames(bio.mat.weights)
genes.GO = rownames(GO)

all_genes = sort(unique( c(data[,1], data[,3])))
intersect_genes = sort(intersect( intersect( intersect(all_genes, genes.coexp), genes.bio), genes.GO))

# ligandI
data.mat = matrix(0, ncol=length(intersect_genes), nrow=length(intersect_genes))
rownames(data.mat) = intersect_genes
colnames(data.mat) = intersect_genes
m.a = match(data[,1],intersect_genes)
f.na = !is.na(m.a)
f.na2 = m.a[f.na]

m.b = match(data[,3], intersect_genes)
f.nb = !is.na(m.b)
f.nb2 = m.b[f.nb]

data.mat[cbind(f.na2, f.nb2)] = 1
diag(data.mat) = 0
ligI = array(data.mat)


# coexp
coexp.mat = matrix(0, ncol=length(intersect_genes), nrow=length(intersect_genes))
m = match(genes.coexp, intersect_genes)
f.n = !is.na(m)
f.n2 = m[f.n]
coexp.mat[f.n2, f.n2] = 1
subnet = network[f.n,f.n]
n = dim(subnet)[1]
diag(subnet) = 0
temp = rankM(subnet,n)
diag(temp)= NA
rankedCo = array(temp)

coexp.mat = matrix(0, ncol=length(all_genes), nrow=length(all_genes))
m = match(genes.coexp, all_genes)
f.n = !is.na(m)
f.n2 = m[f.n]
coexp.mat[f.n2, f.n2] = 1
subnet2 = network[f.n,f.n]

n = dim(subnet2)[1]
diag(subnet2) = 0
temp = rankM(subnet2,n)
diag(temp)= NA
rankedCo2 = array(temp)

data.mat = matrix(0, ncol=length(all_genes), nrow=length(all_genes))
rownames(data.mat) = all_genes
colnames(data.mat) = all_genes
m.a = match(data[,1],all_genes)
f.na = !is.na(m.a)
f.na2 = m.a[f.na]

m.b = match(data[,3], all_genes)
f.nb = !is.na(m.b)
f.nb2 = m.b[f.nb]

data.mat[cbind(f.na2, f.nb2)] = 1
diag(data.mat) = 0
ligICo = array(data.mat[f.n2,f.n2])



# semantic sim.
GOsum = apply(GO, 2, sum)
filt = which(GOsum > 20 & GOsum < 1000)
m.go = match(rownames(GO[,filt]), intersect_genes)
f.go = !is.na(m.go)
f.go2 = m.go[f.go]

num = GO[f.go,filt] %*% t(GO[f.go,filt])
sum = apply(GO[f.go,filt], 1, sum)
sum_rep = repmat( as.matrix(sum), 1, length(sum))
denom = t(sum_rep) + sum_rep - num
J = num/denom

t = which(is.na(J))
J[t] = 0
m.j = match(rownames(J), intersect_genes)
f.nj = !is.na(m.j)
f.nj2 = m.j[f.nj]
subJ = J[f.nj,f.nj]
 diag(subJ) = 0
 temp = rankM(subJ,dim(subJ)[1])
 diag(temp) = NA
 rankedJ = array(temp)

## PPI
m.ppi = match( genes.bio ,intersect_genes)
f.ppi = !is.na(m.ppi)
f.ppi2 = m.ppi[f.ppi]

subPPI = bio.mat.weights[f.ppi,f.ppi]
ligPPI = array(data.mat[f.ppi2, f.ppi2])
diag(subPPI) = 0
 temp = rankM(subPPI,dim(subPPI)[1])
 diag(temp) = NA
 rankedPPI = array(temp)


##
f=5000
png("LigandIncidence.png", width=480, height=480*3)
par(mfrow=c(3,1))
conv_smoother4( rankedCo, ligI, f, coexplab, ylab)
conv_smoother4( rankedJ, ligI, f, semlab, ylab)
conv_smoother_boxplot( rankedPPI, ligI, f, ppilab,ylab )
dev.off()

f = 1e4
zeroes = rep(0,f)


	n <- order(rankedCo)
	m <- rankedCo[n]
	i <- length(m)/f
	temp_b1 = m[ (f/2):(length(m)-f/2))]
	temp_b2 = convolve( ligI[n], rep(1,f),type="filter")
        temp_b2 = temp_b2/f

        temp_b3 = convolve( ligI[n], rep(1,f),type="open")
        temp_b3 = temp_b3/f
        temp_b4 = temp_b3[ (f/2):(length(temp_b3)-f/2)]

sample = sample(length(temp_b1), f)
png("coexp-a.png")
plot(temp_b1[sample], temp_b2[sample])
dev.off()

sample = sample(length(temp_b4), f*10)

png("coexp-b.png")
plot(m[sample], temp_b4[sample])
dev.off()


png("coexp.png")
conv_smoother5(rankedCo, ligI, f, "coexp", "ligand")
dev.off()

png("coexp2.png")
conv_smoother5( rankedCo2, ligI, f, "coexp", "ligand")
dev.off()




all_genes = sort(unique( c(data[,1], data[,3])))
data.mat = matrix(0, ncol=length(all_genes), nrow=length(all_genes))
rownames(data.mat) = all_genes
colnames(data.mat) = all_genes
m.a = match(data[,1],all_genes)
f.na = !is.na(m.a)
f.na2 = m.a[f.na]
m.b = match(data[,3], all_genes)
f.nb = !is.na(m.b)
f.nb2 = m.b[f.nb]
data.mat[cbind(f.na2, f.nb2)] = 1
diag(data.mat) = 0


all_genes = sort(unique( c(data[,1], data[,3])))
data.mat = matrix(0, ncol=length(intersect_genes), nrow=length(intersect_genes))
rownames(data.mat) = intersect_genes
colnames(data.mat) = intersect_genes
m.a = match(data[,1],intersect_genes)
f.na = !is.na(m.a)
f.na2 = m.a[f.na]
m.b = match(data[,3], intersect_genes)
f.nb = !is.na(m.b)
f.nb2 = m.b[f.nb]
data.mat[cbind(f.na2, f.nb2)] = 1
diag(data.mat) = 0
ligI = array(data.mat) 

png("coexp-testing.png")
conv_smoother4( rankedCo,ligI, 5000, "coexp", "ligand")
dev.off()

ligICo2 = array(subnetlig)
png("coexp-testing.png")
conv_smoother4( rankedCo2,ligICo2, 10000, "coexp", "ligand")
dev.off()


