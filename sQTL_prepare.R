library(dplyr)
library(data.table)
library(leafcutter)
library(stringr)
library(readr)
library(qvalue)
library(pbmcapply)
library(pbapply)

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
  as.data.frame(fread(cmd = paste("zless", exon_file)) )
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
introns$chr <- paste0("chr", introns$chr)

####################
# ANNOTATE JUNCTIONS
####################

# functions now in separate script
#source("~/Desktop/sQTLviz-master/sQTL_annotation_functions.R")

intersect_introns <- function(introns) {	
  all.introns <- introns
  # for each splice site write out a bed file  	
  all.junctions <- dplyr::select(all.introns, chr, start, end, clusterID = clu)

  all.fiveprime <- data.frame(chr = all.introns$chr,	
                              start = all.introns$start,	
                              end = all.introns$start + 1,	
                              clusterID = all.introns$clu)

  all.threeprime <- data.frame(chr = all.introns$chr,	
                              start = all.introns$end,	
                              end = all.introns$end + 1,	
                              clusterID = all.introns$clu)

  all.file <- "all_junctions.bed"	
  all.fiveprime.file <- "all_fiveprime.bed"	
  all.threeprime.file <- "all_threeprime.bed"	
  
  write.table(all.junctions, all.file, col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t" )	
  write.table(all.threeprime, all.threeprime.file, col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t" )	
  write.table(all.fiveprime, all.fiveprime.file, col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t" )	
  
  print( "BedTools intersect junctions with list of known splice sites")	
  
  # first match junctions
  all.introns.cmd <- paste0("/anaconda3/bin/bedtools intersect -a ", all.file, " -b ", all_introns, " -wa -wb -loj -f 1 -nonamecheck")	
  all.introns_intersect <- fread(cmd = all.introns.cmd)
  
  # intersect with bedtools to find the annotations of each splice site	
  threeprime.cmd <- paste0("/anaconda3/bin/bedtools intersect -a ", all.threeprime.file, " -b ", threeprime_file, " -wa -wb -loj -f 1 -nonamecheck" )
  threeprime_intersect <- fread(cmd = threeprime.cmd)	
  
  fiveprime.cmd <- paste0("/anaconda3/bin/bedtools intersect -a ", all.fiveprime.file, " -b ", fiveprime_file, " -wa -wb -loj -f 1 -nonamecheck" )	
  fiveprime_intersect <- fread(cmd = fiveprime.cmd)	
  
  # remove temporary files	
  rm.cmd <- paste("rm ", all.file, all.fiveprime.file, all.threeprime.file) 	
  system(rm.cmd)
  
  return(list(threeprime_intersect, fiveprime_intersect, all.introns_intersect))	
}

intersects <- intersect_introns(introns)
threeprime_intersect <- intersects[[1]]
fiveprime_intersect <- intersects[[2]]
all.introns_intersect <- intersects[[3]]

print("Annotating junctions")
uniqueClusters <- unique(introns$clu) 

annotate_single_cluster <- function(introns, input_clu, cluIndex){	
  # for each intron in the cluster, check for coverage of both	
  # output a vector of string descriptions 	
  cluster <- introns[clu == input_clu]
  cluster$start <- as.integer(cluster$start)
  cluster$end <- as.integer(cluster$end)
  # subset intersects by clusterID (V4)	
  fprimeClu <- fiveprime_intersect %>%
    filter(., V4 == input_clu) %>%
    as.data.table()
  tprimeClu <- threeprime_intersect %>%
    filter(., V4 == input_clu) %>%
    as.data.table()
  bothSSClu <- all.introns_intersect %>%
    filter(., V4 == input_clu) %>%
    as.data.table()
  # for each intron in the cluster:	
  # create vector of overlapping splice sites, indexed by the row of the intersect	
  # five prime splice sites	
  fprime <- apply(cluster, MAR = 1, FUN = function(x) {	
    fprimeClu[paste(V1, V2, sep = "_") %in% paste(cluster$chr, cluster$start, sep = "_")]
  })	
  # three prime splice sites	
  tprime <- apply(cluster, MAR = 1, FUN = function(x) {	
    tprimeClu[paste(V1, V2, sep = "_") %in% paste(cluster$chr, cluster$end, sep = "_")]
  })
  # both splice sites	
  bothSS <- apply(cluster, MAR = 1, FUN = function(x) {	
    bothSSClu[paste(V1, V6, V7, sep = "_") %in% paste(cluster$chr, cluster$start, cluster$end, sep = "_")]
  })
  
  # find gene and ensemblID by the most represented gene among all the splice sites	
  cluster_genes <- names(sort(table(do.call(rbind, c(fprime, tprime, bothSS))$V8), decreasing = TRUE))
  
  cluster_gene <- cluster_genes[cluster_genes != "." ][1]	
  # if no cluster gene found then leave as "."	
  if( length(cluster_gene) == 0){	
    cluster_gene == "."	
  }	
  # do the same for EnsemblID	
  cluster_ensemblIDs <- names(sort(table(do.call(rbind, c(fprime, tprime, bothSS))$V9), decreasing = TRUE))	
  cluster_ensemblID <- cluster_ensemblIDs[ cluster_ensemblIDs != "." ][1]	
  if( length( cluster_ensemblID ) == 0 ){	
    cluster_ensemblID == "."	
  }	
  
  verdict <- c()
  coord <- c()
  gene <- c()
  ensemblID <- c()
  transcripts <- list()
  
  for(intron in 1:nrow(cluster)) {	
    coord[intron] <- paste0(cluster[intron,]$chr, ":", cluster[intron,]$start, "-", cluster[intron,]$end)
    gene[intron] <- cluster_gene
    ensemblID[intron] <- cluster_ensemblID
    
    # for each intron create vector of all transcripts that contain both splice sites	
    # TO DO : these do not cover all situations!
    # EDIT COMPLETE: change all(fprime[[intron]]$V5 != ".") to !all(fprime[[intron]]$V5 == ".")
    # check whether 3' or 5' splice sites are present among annotated introns
    # first part of table is the leafcutter intron, second part is the gencode transcript
    
    verdict[intron] <- "error"
    if (all(tprime[[intron]]$V5 == ".") & all(fprime[[intron]]$V5 == ".")) { 
      verdict[intron] <- "cryptic_unanchored"	# neither annotated
    } 
    if (all(tprime[[intron]]$V5 == ".") & !all(fprime[[intron]]$V5 == ".")) { 
      verdict[intron] <- "cryptic_threeprime"	# one annotated
    }	
    if (!all(tprime[[intron]]$V5 == ".") & all(fprime[[intron]]$V5 == ".")) { 
      verdict[intron] <- "cryptic_fiveprime" # one annotated
    }	
    if (!all(tprime[[intron]]$V5 == ".") & !all(fprime[[intron]]$V5 == "."))	{ 	
      # test if the splice sites are paired in a known intron	
      if (!all(bothSS[[intron]]$V5 == ".")) {	
        verdict[intron] <- "annotated" # both annotated
      } else { # both are annotated but never in the same junction	
        verdict[intron] <- "novel annotated pair"	
      }	
    }	
  }
  if (cluIndex %% 1 == 100) {	
    print(paste("processed", cluIndex, "clusters" ))	
  }	
  return(	
    data.table(	
      clusterID = input_clu,	
      coord = coord,	
      gene = gene,	
      ensemblID = ensemblID,	
      verdict = verdict)	
  )	
}

annotatedClusters <- do.call(rbind, 
        pblapply(1:length(uniqueClusters), 
               function(x) annotate_single_cluster(introns, input_clu = uniqueClusters[x], cluIndex = x))
        )

annotatedClusters[is.na(gene)]$gene <- "."
annotatedClusters[is.na(ensemblID)]$ensemblID <- "."


#################
# PREPARE RESULTS
#################

get_snp_meta <- function(snps) {
  snp_meta <-  data.table(str_split_fixed(snps, "_", 5), snps) %>%
    setnames(., c("snp_chr", "snp_pos", "snp_ref", "snp_alt", "genome_build", "snp_ID"))
  snp_meta$snp_pos <- as.numeric(snp_meta$snp_pos)
  return(snp_meta)
}

## the results table should consist of a list of clusters with their most significant SNP
# bind intron_meta, intron, snp, p value, snp_meta

# use associations with q < 0.10 cut-off
res <- res[qval < 0.1]
# Bind together metadata with original results
sigJunctions <- cbind(get_intron_meta(res$intron_cluster), 
                      res[, c("intron_cluster", "variant_id", "adj_p")],
                      get_snp_meta(res$variant_id),
                      res[, "qval"]) %>%
  as.data.table()

# sometimes there will be duplicates - remove!
sigJunctions <- dplyr::distinct(sigJunctions) %>%
  setorder(., adj_p)

#names(sigJunctions)[8] <- "bpval"
# present most significant junction for each SNP?
# or most significant SNP for each junction?
resultsByCluster <- dplyr::group_by(sigJunctions, clu) %>% 
  dplyr::summarise(chr = first(chr),
                   start = min(start),
                   end = max(end),
                   snp = first(snp_ID),
                   snp_chr = first(snp_chr),
                   pos = first(snp_pos),
                   pval = first(adj_p),
                   qval = first(qval)) %>%
  as.data.table() %>%
  setorder(., pval)

####
## PREPARE FOR SHINY
####

code <- "test"
annotation_code <- "gencode_hg19"

resultsByCluster$gene <- annotatedClusters$gene[match(resultsByCluster$clu, annotatedClusters$clusterID)]
resultsByCluster$SNP_pos <- paste0(resultsByCluster$snp_chr, ":", resultsByCluster$pos)

# fix coords without "chr"
# TO DO: this is a weird way to check this...
if (all(!grepl("chr", sample(resultsByCluster$chr, 100)))) {
  resultsByCluster$chr <- paste0("chr", resultsByCluster$chr)
}

resultsByCluster[, cluster_pos := paste0(chr, ":", start, "-", end)]

resultsToPlot <- select(resultsByCluster,
  SNP = snp,
  SNP_pos, 
  gene = gene,
  cluster_pos,
  q = qval) %>%
  as.data.frame()

row.names(resultsToPlot) <- resultsByCluster$clu
resultsToPlot$q <- signif(resultsToPlot$q,  digits = 3)

#######

# get the Betas and per-junction q values
# add in full permutation results to get Beta for each junction
permutation_full_res <- "permutations_full.txt.gz"
perm_full <- res
setnames(perm_full, c("beta", "qval", "intron_cluster", "variant_id"), c("Beta", "FDR", "clusterID", "SNP"))

perm_clean <- select(perm_full, clusterID, SNP, Beta, FDR)
perm_clean <- cbind(perm_clean,
                    get_intron_meta(perm_full$clusterID),
                    get_snp_meta(perm_full$SNP))


# junction table - each junction with Beta, P value and annotation
junctionTable <- resultsToPlot %>%
  mutate( clu = row.names(resultsToPlot) ) %>%
  left_join(introns_to_plot, by = "clu" ) %>%
  rename(snp_ID = SNP) %>%
  left_join( perm_clean, 
             by = c("chr" = "snp_chr", "start", "end", "snp_ID", "clu", "middle")
  ) %>%
  mutate(coord = paste0( chr, ":", start, "-", end)) %>%
  left_join(annotatedClusters, 
            by = c("clu" = "clusterID", "coord" )
  ) %>%
  select(clu, coord, verdict, Beta, q = FDR) %>%
  mutate( Beta = signif(Beta, digits = 3),
          q = signif(q, digits = 3)) %>%
  mutate( Beta = ifelse(is.na(Beta), ".", Beta),
          q = ifelse(is.na(q), ".", q))

  introns_to_plot$chr <- paste0("chr", introns_to_plot$chr)
  
save.image("all_data.Rdata")
print("saving objects")
save(annotatedClusters, # every junction needed
     sigJunctions, # every junction x SNP interaction
     resultsToPlot, #significant clusters and the most significant SNP
     clusters, # junction counts for each sample
     vcf, # the genotypes of each sample
     vcf_meta, # the vcf metadata
     introns_to_plot, # all the intron positions
     exons_table, # the annotation
     junctionTable, # the junctions to display for each cluster
     annotation_code,
     code,
     file = paste0("sQTL_results.Rdata")
)

