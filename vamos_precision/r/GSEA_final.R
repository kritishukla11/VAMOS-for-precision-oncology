#file names
inputge = "GOZGIT_ESR1_TARGETS_UP_TARGET_GE.csv"
input_gene_set = "GOZGIT_ESR1_TARGETS_UP.csv"
output_file = "GOZGIT_ESR1_TARGETS_UP_GSEA.csv"

#load packages
library(matrixStats)
library(circlize)
library(data.table)

#define function
ssgsea = function(X, gene_sets, alpha = 0.25, scale = T, norm = F, single = T) {
    row_names = rownames(X)
    num_genes = nrow(X)
    gene_sets = lapply(gene_sets, function(genes) {which(row_names %in% genes)})
    
    # Ranks for genes
    R = matrixStats::colRanks(X, preserveShape = T, ties.method = 'average')
    
    # Calculate enrichment score (es) for each sample (column)
    es = apply(R, 2, function(R_col) {
        gene_ranks = order(R_col, decreasing = TRUE)
        
        # Calc es for each gene set
        es_sample = sapply(gene_sets, function(gene_set_idx) {
            # pos: match (within the gene set)
            # neg: non-match (outside the gene set)
            indicator_pos = gene_ranks %in% gene_set_idx
            indicator_neg = !indicator_pos
            
            rank_alpha  = (R_col[gene_ranks] * indicator_pos) ^ alpha
            
            step_cdf_pos = cumsum(rank_alpha)    / sum(rank_alpha)
            step_cdf_neg = cumsum(indicator_neg) / sum(indicator_neg)
            
            step_cdf_diff = step_cdf_pos - step_cdf_neg
            
            # Normalize by gene number
            if (scale) step_cdf_diff = step_cdf_diff / num_genes
            
            # Use ssGSEA or not
            if (single) {
                sum(step_cdf_diff)
            } else {
                step_cdf_diff[which.max(abs(step_cdf_diff))]
            }
        })
        unlist(es_sample)
    })
    
    if (length(gene_sets) == 1) es = matrix(es, nrow = 1)
    
    # Normalize by absolute diff between max and min
    if (norm) es = es / diff(range(es))
    
    # Prepare output
    rownames(es) = names(gene_sets)
    colnames(es) = colnames(X)
    return(es)
}

#perform GSEA
data1 = read.csv(input_ge)
data2 = data1[ , !(names(data1) %in% c("Name"))]
data = as.matrix(data2)
rownames(data)<-paste0(data1[,1])
gene_set = read.csv(input_gene_set)
gene_sets = as.list(as.data.frame(gene_set))
system.time(assign('res', ssgsea(data, gene_sets, scale = TRUE, norm = FALSE)))
res1 = t(res)
head(res1)
write.csv(res1,output_file)
