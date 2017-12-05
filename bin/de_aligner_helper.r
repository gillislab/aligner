
get_roc <- function (scores, labels)
{
    o <- order(scores, decreasing = TRUE)
    h1 <- labels[o]
    h2 <- !labels[o]
    tpr <- cumsum(h1)/sum(h1)
    fpr <- cumsum(h2)/sum(h2)
    return(cbind(fpr, tpr))
}

 auroc_analytic <- function (scores, labels)
{
    negatives <- which(labels == 0, arr.ind = TRUE)
    scores[negatives] <- 0
    p <- sum(scores, na.rm = TRUE)
    nL <- length(labels)
    np <- sum(labels, na.rm = TRUE)
    nn <- nL - np
    auroc <- (p/np - (np + 1)/2)/nn
    return(auroc)
}



get_enrichment_rocs <-function(de, sub.labels, f.e){
	score = (rank(-log10(de[f.e,4]))+rank(abs(de[f.e,2])))/2
	curves = lapply(1:dim(sub.labels)[2], function(i) get_roc(score, sub.labels[f.e,i]) )
	aurocs = lapply(1:dim(sub.labels)[2], function(i) auroc_analytic(score, sub.labels[f.e,i]) )
	names(aurocs) = colnames(sub.labels)
	names(curves) = colnames(sub.labels)
	return( list(curves, aurocs) )
}



calc_DE <- function(X, group){

	m.X1 = rowMeans(X[,group==1])
	m.X2 = rowMeans(X[,group==2])
        m.X = rowMeans(X)
	fc = log2(m.X1/m.X2)

	X.ps = sapply(1:dim(X)[1], function(k) wilcox.test(X[k,group==1], X[k,group==2])$p.val)
	X.padj = p.adjust(X.ps , method = "BH")

	de = cbind(m.X, fc, X.ps, X.padj)

	return(de)
}


calc_DE_balanced <- function(X, group){

	nN = count(group)[which.min(count(group)[,2]),2]
	subsample1 = sample(which(group == 1), nN)
        subsample2 = sample(which(group == 2), nN)

	m.X1 = rowMeans(X[,subsample1])
	m.X2 = rowMeans(X[,subsample2])
        m.X = rowMeans(X[,c(subsample1,subsample2) ])
	fc = log2(m.X1/m.X2)

	X.ps = sapply(1:dim(X)[1], function(k) wilcox.test(X[k,subsample1], X[k,subsample2])$p.val)
	X.padj = p.adjust(X.ps , method = "BH")

	de = cbind(m.X, fc, X.ps, X.padj)

	return(de)
}


calc_DE_prop <- function(X, group, prop){

	X = X[,prop]

	nN = count(group[prop])[which.min(count(group[prop])[,2]),2]
	subsample1 = sample(which(group[prop] == 1), nN)
        subsample2 = sample(which(group[prop] == 2), nN)

	m.X1 = rowMeans(X[,subsample1] )
	m.X2 = rowMeans(X[,subsample2])
        m.X = rowMeans(X[,c(subsample1,subsample2)])
	fc = log2(m.X1/m.X2)

	X.ps = sapply(1:dim(X)[1], function(k) wilcox.test(X[k,subsample1], X[k,subsample2])$p.val)
	X.padj = p.adjust(X.ps , method = "BH")

	de = cbind(m.X, fc, X.ps, X.padj)

	return(de)
}


get_enrichment_rocs <-function(de, sub.labels, f.e){
	score = (rank(-log10(de[f.e,4]))+rank(abs(de[f.e,2])))/2
	curves = lapply(1:dim(sub.labels)[2], function(i) get_roc(score, sub.labels[f.e,i]) )
	aurocs = lapply(1:dim(sub.labels)[2], function(i) auroc_analytic(score, sub.labels[f.e,i]) )
	names(aurocs) = colnames(sub.labels)
	names(curves) = colnames(sub.labels)
	return( list(curves, aurocs) )
}


