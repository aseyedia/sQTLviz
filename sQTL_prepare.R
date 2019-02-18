library(dplyr)
library(data.table)
library(leafcutter)
library(stringr)
library(readr)
library(qvalue)

setwd("~/Desktop/sQTLviz-master/")

permutation_res <- "permutations_full.txt.gz"

# other files
VCF = "data/genotypes_MAF1.vcf.gz"
clusters_table <- "Ne-sQTL_perind_numers.counts.gz" # now just the counts - don't need to strip ratios

# PRE-REQUISITES:
# annotation
annotation_code <- "gencode_hg19"
exon_file <- paste0(annotation_code, "_all_exons.txt.gz")
all_introns <- paste0(annotation_code,"_all_introns.bed.gz" )
threeprime_file <- paste0( annotation_code,"_threeprime.bed.gz")
fiveprime_file <- paste0( annotation_code,"_fiveprime.bed.gz")

exons_table <- if (!is.null(exon_file)) {
  cat("Loading exons from", exon_file, "\n")
  as.data.frame(fread(paste("zless", exon_file)) )
} else {
  cat("No exon_file provided.\n")
  NULL
}

## CHECK VCFTOOLS IS INSTALLED
if( !file.exists(system("which vcftools", intern = TRUE) )){
  stop("vcftools is not in your $PATH - have you installed it?")
}
## BEDTOOLS
if( !file.exists(system("which bedtools", intern = TRUE) )){
  stop("bedtools is not in your $PATH - have you installed it?")
}

# read in clusters
print("reading in clusters")
clusters <- read.table(clusters_table, header = TRUE)

# find and harmonise sample names
samples <- names(clusters)[2:ncol(clusters)]
srr2gtex <- fread("tissue_table.txt")
gtex_samples <- srr2gtex[Run %in% samples]$submitted_subject_id

# write out samples
samples_file <- "used_samples.txt"
write.table(gtex_samples, samples_file, col.names = FALSE, row.names = FALSE, quote = FALSE)

# read in junction x snp results
print("reading in results")
res <- fread(permutation_res, header = FALSE) %>%
  setnames(., c("intron_cluster", "chrom", "pheno_start", "pheno_end", 
                "strand", "total_cis", "distance", "variant_id", "variant_chrom", 
                "var_start", "var_end", "df", "dummy", "param_1", "param_2",
                "p", "beta", "emp_p", "adj_p"))

res[, qval := qvalue(res$adj_p)$qvalues]

sig_fdr_10 <- unique(res[qval < 0.1]$variant_id)

sig_snp_file <- "sig_snps.txt"
write.table(sig_fdr_10, file = sig_snp_file, col.names = FALSE, row.names = FALSE, quote = FALSE)

######
# PREPARE GENOTYPES
######

# use vcftools to filter out just the snps and samples required
print("filtering VCF")
#cmd <- paste("vcftools --gzvcf", VCF,  "--snps", sig_snp_file, "--keep", samples_file, "--recode --stdout")
#vcf <- fread(cmd)
vcf <- fread("grep -v '^##' out.recode.vcf")


# round genotypes to 0,1,2
vcf <- vcf[complete.cases(vcf),]
vcf[, 10:ncol(vcf)] <- as.data.frame(apply( vcf[, 10:ncol(vcf)], MAR = c(1, 2), FUN = function(x) gsub(":.*", "", x)))
vcf[, 10:ncol(vcf)] <- as.data.frame(apply( vcf[, 10:ncol(vcf)], MAR = c(1, 2), FUN = function(x) if (x == "0/0") {return(0)} else if (x == "0/1") {return(1)} else if (x == "1/1") {return(2)} else if (x == "./.") {return(NA)} ))

get_vcf_meta <- function(vcf) {
  # take a VCF and split out the meta data 
  # to use for querying
  vcf_meta <- data.frame(SNP = vcf$ID, SNP_pos = paste(vcf$`#CHROM`, vcf$POS, sep = ":"), REF = vcf$REF, ALT = vcf$ALT, stringsAsFactors = FALSE)
  return(vcf_meta)
}

vcf_meta <- get_vcf_meta(vcf)

##################
# PREPARE CLUSTERS
##################

# from significant associations
sigClusters <- unique(sapply(strsplit(unique(res[qval < 0.1]$intron_cluster), ":"), function(x) x[4]))

introns <- get_intron_meta(row.names(clusters)) %>% as.data.table()
sig_index <- which(sapply(strsplit(rownames(clusters), ":"), function(x) x[4]) %in% sigClusters)

# remove non-significant clusters
introns <- introns[clu %in% sigClusters]
clusters <- clusters[sig_index,]

swap_id <- data.table(Run = colnames(clusters))
swap_id[, index := .I]
swap_id <- merge(swap_id, srr2gtex, "Run") %>%
  setorder(., index)
new_ids <- swap_id$submitted_subject_id
colnames(clusters) <- new_ids

# rearrange sample columns in clusters so they match the VCF
samples <- names(vcf)[10:ncol(vcf)]
clusters <- clusters[, samples]

introns_to_plot <- get_intron_meta(row.names(clusters))

# for each cluster work out mean proportion of each junction
# remove junctions < 1% contribution

juncProp <- function(cluster){
  cluster$prop <- cluster$meanCount / sum(cluster$meanCount) 
  return(cluster)
}

splitClusters <- introns_to_plot %>%
  mutate( 
    clu = factor(.$clu, levels = unique(.$clu)),
    meanCount = rowMeans(clusters) ) %>%
  split( .$clu ) %>%
  purrr::map_df( juncProp ) %>%
  mutate( clu = as.character(.$clu))

introns_to_plot <- introns_to_plot[splitClusters$prop >= 0.01,]
clusters <- clusters[splitClusters$prop >= 0.01,]
introns <- introns[splitClusters$prop >= 0.01,]

####################
# ANNOTATE JUNCTIONS
####################

# functions now in separate script
source("~/Desktop/sQTLviz-master/sQTL_annotation_functions.R")

intersects <- intersect_introns(introns)
threeprime_intersect <- intersects[[1]]
fiveprime_intersect <- intersects[[2]]
all.introns_intersect <- intersects[[3]]

print("Annotating junctions")

uniqueClusters <- unique( introns$clu ) 
#uniqueClusters <- uniqueClusters[7010:7020]
#annotatedClusters <- lapply( seq_along(uniqueClusters),
#                             FUN = function(i) annotate_single_cluster(introns, clu = uniqueClusters[i], cluIndex = i) )
#annotatedClusters <- do.call(rbind, annotatedClusters)

# replace with purrr?
annotatedClusters <- purrr::map_df( seq_along(uniqueClusters),
                                    ~annotate_single_cluster( introns, clu = uniqueClusters[.], cluIndex = .  ) )

annotatedClusters$gene[ is.na( annotatedClusters$gene) ] <- "."
annotatedClusters$ensemblID[ is.na( annotatedClusters$ensemblID) ] <- "."


#################
# PREPARE RESULTS
#################

get_snp_meta <- function(snps){
  snp_meta <-  as.data.frame(str_split_fixed(snps, "\\.", 3), stringsAsFactors = FALSE)
  colnames(snp_meta) <- c("snp_ID", "snp_chr", "snp_pos")
  # a few snps don't have any IDs - just the coordinates
  noID <- snp_meta[ grepl("\\.", snp_meta$snp_pos),]
  noID <- select(noID, snp_pos, snp_ID, snp_chr)
  names(noID) <- c("snp_ID", "snp_chr","snp_pos")
  # put back in
  snp_meta[ grepl("\\.", snp_meta$snp_pos),] <- noID
  
  snp_meta$snp_pos <- as.numeric(snp_meta$snp_pos)
  return(snp_meta)
}

## the results table should consist of a list of clusters with their most significant SNP
# bind intron_meta, intron, snp, p value, snp_meta
#sigJunctions <- cbind( get_intron_meta( res[,1]), res[, c(1,6,11)], get_snp_meta(res[,6]))

# use associations with q < 0.05 cut-off. 
# Bind together metadata with original results
sigJunctions <- cbind( get_intron_meta( res[,1]), 
                       res[, c(1,2,3)],
                       get_snp_meta(res[,2]))


# bind together - all will be accessible in same table
names(NallsJunctions) <- names(sigJunctions)
names(YangJunctions) <- names(sigJunctions)
sigJunctions <- rbind(sigJunctions, NallsJunctions)
sigJunctions <- rbind(sigJunctions, YangJunctions)
# sometimes there will be duplicates - remove!
sigJunctions <- dplyr::distinct(sigJunctions)

#names(sigJunctions)[8] <- "bpval"
# present most significant junction for each SNP?
# or most significant SNP for each junction?
resultsByCluster <- dplyr::group_by(sigJunctions[order(sigJunctions$bpval),], clu) %>% 
    dplyr::summarise( chr = first(chr),
                      start = min(start),
                      end = max(end),
                      snp = first(snp_ID),
                      snp_chr = first(snp_chr),
                      pos = first(snp_pos),
                      FDR = first(bpval) ) %>%
    dplyr::arrange(FDR)

# for Nalls and Yang SNPs don't thin at all - we want those SNPs!

####
## PREPARE FOR SHINY
####

code <- "test"
annotation_code <- "gencode_hg19"

resultsByCluster$gene <- annotatedClusters$gene[ match(resultsByCluster$clu, annotatedClusters$clusterID)]
resultsByCluster$SNP_pos <- paste0(resultsByCluster$snp_chr, ":", resultsByCluster$pos)

# fix coords without "chr"
if( all( !grepl("chr", sample(resultsByCluster$chr, 100)) ) ){
  resultsByCluster$chr <- paste0("chr", resultsByCluster$chr)
}

resultsByCluster$cluster_pos = paste0(resultsByCluster$chr,":", resultsByCluster$start,"-",resultsByCluster$end)

resultsToPlot <- as.data.frame( select( resultsByCluster,
                           SNP = snp,
                           SNP_pos, 
                           gene = gene,
                           cluster_pos,
                           q = FDR
) )

row.names(resultsToPlot) <- resultsByCluster$clu
resultsToPlot$q <- signif(resultsToPlot$q,  digits = 3)



# get the Betas and per-junction q values
# add in full permutation results to get Beta for each junction
permutation_full_res <- "data/permutations.all.CMC.txt.gz"
perm_full <- read_delim( permutation_full_res,
                         col_names = c("clusterID", "V2","V3","V4","V5","SNP","V7","V8","Beta","V10","FDR"),
                         delim = " "
)
perm_clean <- select(perm_full, clusterID, SNP, Beta, FDR)
perm_clean <- cbind( perm_clean,
                     get_intron_meta(perm_full$clusterID),
                     get_snp_meta(perm_full$SNP) )


# junction table - each junction with Beta, P value and annotation
junctionTable <- resultsToPlot %>%
  mutate( clu = row.names(resultsToPlot) ) %>%
  left_join(introns_to_plot, by = "clu" ) %>%
  rename(snp_ID = SNP) %>%
  left_join( perm_clean, 
             by = c("chr" = "snp_chr", "start", "end", "snp_ID", "clu", "middle")
  ) %>%
  mutate(coord = paste0( chr, ":", start, "-", end)) %>%
  left_join( annotatedClusters, 
             by = c("clu" = "clusterID", "coord" )
  ) %>%
  select(clu, coord, verdict, Beta, q = FDR) %>%
  mutate( Beta = signif(Beta, digits = 3),
          q = signif(q, digits = 3)) %>%
  mutate( Beta = ifelse(is.na(Beta), ".", Beta),
          q = ifelse(is.na(q), ".", q))



save.image("data/all_data.Rdata")
print("saving objects")
save( annotatedClusters, # every junction needed
      sigJunctions, # every junction x SNP interaction
      resultsToPlot, #significant clusters and the most significant SNP
      GWASresults, # associations with SNPs from a PD GWAS
      YangResults, # associations with Yang's TWAS hit SNPs
      clusters, # junction counts for each sample
      vcf,# the genotypes of each sample
      vcf_meta, # the vcf metadata
      introns_to_plot, # all the intron positions
      #counts, 
      #meta, 
      exons_table, # the annotation
      junctionTable, # the junctions to display for each cluster
      #pca, 
      #intron_summary, 
      #cluster_summary, 
      #introns_to_plot,
      #cluster_ids,
      #sample_table,
      annotation_code,
      code,
      file = paste0( "sQTLviz/sQTL_results.Rdata")
)

# to do - cut down size of exon table to increase speed of querying
#allGenes <- unique(resultsToPlot$gene)
#exons_table_cut <- exons_table[ exons_table$gene_name %in% allGenes ,]

