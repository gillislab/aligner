
# counts => count reads per gene/transcript for a single sample, matrix n_genes x 1, repeat for multiple samples

# Downsampling fractions
frac_down = 10:19  / 2  * 10

# Setting up variables
n_genes = dim(counts)[1]
total_reads_mapped = sum(counts)
reads_down = total_reads_mapped -  colSums(X) * frac_down/100
colnames(reads_down) = frac_down/100
rownames(reads_down) = colnames(X)

# Initialize downsampling matrix
X = matrix(counts, ncol= length(frac_down), nrow=n_genes )

# Remove counts at random from genes/transcriptsfor each downsampling fraction
f = counts > 0
i = 1
for( size in reads_down ){
    mask = sample( total_reads_mapped,  size )
    b = c(0, cumsum(counts[f,1]) )
    h = hist( mask, breaks=b )
    X[f,i] = counts[f,1]-h$counts
    i = i + 1
}
