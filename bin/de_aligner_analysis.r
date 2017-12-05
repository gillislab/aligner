# Load gene sets and attributes
load("gene_annotations_v22.Rdata")
load("DE_sex_diff_genesets.Rdata")

# Helper code and functions
source("de_aligner_helper.r")
require(EGAD)


# samples
# load geuvadis data set at default parameter 0.66
load("sample_counts.Rdata")
series = c("0.66")  # change if multiple parameters tested

group = as.numeric(sex)
nS = length(samples)
n_rep = 10    # number of repeats
n_sub = 200   # number of samples to sub sample to

# Filter for genes with entrezIDs and non zeroes
f.k = !is.na(attr$entrezID) & rowSums(counts_exp) > 0


deg.sub = list()
# for( i in 1:length(series)){
	X = counts_exp
	rownames(X) = attr$ensemblID
	colnames(X) = samples
	N = colSums(X,na.rm=T)
	# convert counts to CPM
	X.cpm = sapply(1:nS, function(k) 10^6*X[,k]/N[k] )
	        deg.sub[[i]] = list()
        	for( j in 1:n_rep){
        		sub = c(sample(which(group==1),n_sub ), sample(which(group==2),n_sub ))
			X = X.cpm[,sub]
                	m.X = rowMeans(X)

			m.X1 = rowMeans(X[,sex[sub] == "female"])
			m.X2 = rowMeans(X[,sex[sub] == "male"])
			fc = log2(m.X1/m.X2)

			X.ps = sapply(1:dim(X)[1], function(k) wilcox.test(X[k,sex[sub] == "female"], X[k,sex[sub] == "male" ])$p.val)
	       		X.padj = p.adjust(X.ps , method = "BH")

			de = cbind(m.X, fc, X.ps, X.padj)

			deg.sub[[i]][[j]]= list()
	        	deg.sub[[i]][[j]]$de = de
		}
# }

save(deg.sub, file="deg.sub.wilcox.Rdata")

# random DE (shuffle samples)
deg.sub.rand = list()
i = 1
# for( i in 1:length(series)){
	X = counts_exp
	rownames(X) = attr$ensemblID
	colnames(X) = samples
	N = colSums(X,na.rm=T)
	# convert counts to CPM
	X.cpm = sapply(1:nS, function(k) 10^6*X[,k]/N[k] )

                deg.sub.rand[[i]] = list()
        	j = 1
                for( j in 1:n_rep){
        		sub = c(sample(which(group==1),n_sub ), sample(which(group==2),n_sub ))
                        range = 1:(n_sub/2)
                        sub1 = c(sub[range],sub[ range + 2*(n_sub/2) ])
                        sub2 = c(sub[range +(n_sub/2)], sub[range + 3*(n_sub/2) ])

                        X = X.cpm[,sub]
                	m.X = rowMeans(X)

			m.X1 = rowMeans(X.cpm[,sub1 ])
			m.X2 = rowMeans(X.cpm[,sub2])
			fc = log2(m.X1/m.X2)

			X.ps = sapply(1:dim(X)[1], function(k) wilcox.test(X.cpm[k,sub1], X.cpm[k,sub2])$p.val)
	       		X.padj = p.adjust(X.ps , method = "BH")

			de = cbind(m.X, fc, X.ps, X.padj)

			deg.sub.rand[[i]][[j]]= list()
	        	deg.sub.rand[[i]][[j]]$de = de
		}
	deg.sub.rand.sub = deg.sub.rand[[i]]

        save(deg.sub.rand.sub,i, file=paste0(i,"deg.sub.rand.wilcox.Rdata") )

# }

save(deg.sub.rand, file="deg.sub.rand.wilcox.Rdata")



# Calculate enrichments, sub.labels is in the DE_sex_diff file 
rocs.all.r = list()
for( i in 1:length(series) ) {
        deg.adj = deg.sub.rand[[i]][[1]][[1]]
        rocs.all.r[[i]] = get_enrichment_rocs(deg.adj, sub.labels, f.k)
}
rocs.table.r = sapply( 1:length(series), function(i) unlist(rocs.all.r[[i]][[2]] ))


rocs.all= list()
for( i in 1:length(series) ) {
        deg.adj = deg.sub[[i]][[1]][[1]]
        rocs.all[[i]] = get_enrichment_rocs(deg.adj, sub.labels, f.k)
}
rocs.table= sapply( 1:length(series), function(i) unlist(rocs.all[[i]][[2]] ))




