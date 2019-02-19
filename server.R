library(shiny)
library(dplyr)
library(ggplot2)
library(DT)
library(leafcutter)
library(reshape2)
library(gridExtra)
library(intervals)
library(foreach)
library(shinycssloaders)
library(grid)
library(gtable)
library(ggrepel)
library(ggbeeswarm)
library(stringr)

if (!exists("introns")){
  load("sQTL_results.Rdata")
  defaultValue <- 1
} else {
  defaultValue <- NULL
}

source("make_sQTL_cluster_plot.R")
source("make_sQTL_gene_plot.R")
source("make_sQTL_box_plot.R")

shinyServer(function(input, output) {
  
  # ALL CLUSTER X SNP TABLE
  output$all_clusters <- DT::renderDataTable({
    # clicked tab gives you the thing to subset
    subsetChoice <- eval(parse(text = input$navBarPage))
    neand <- fread("tag_snps.neand.EUR.bed") %>%
      mutate(., var_id_1 = paste(V1, V3, V4, V5, "b37", sep = "_")) %>%
      mutate(., var_id_2 = paste(V1, V3, V5, V4, "b37", sep = "_")) %>%
      as.data.table()
    neand_list <- c(neand$var_id_1, neand$var_id_2)
    resultsToPlot$is_neand <- (resultsToPlot$SNP %in% neand_list)
    df <- subset(resultsToPlot, row.names(resultsToPlot) %in% row.names(subsetChoice) )
    datatable(df,
              rownames = FALSE,
              selection = 'single',
              fillContainer = FALSE ,
              options = list(language = list(searchPlaceholder = "for a SNP or gene..."))
    )
  })
  
  # JUNCTION TABLE
  output$junctionTable <- DT::renderDataTable({
    jtable <- dplyr::filter(junctionTable, clu == mydata()$clusterID) %>%
      dplyr::select(-clu)
    datatable(jtable,
              escape = FALSE,
              rownames = FALSE,
              fillContainer = FALSE,
              options <- list( searching = FALSE, paging = FALSE, info = FALSE ))
  })
  
  # REACTIVITY
  values <- reactiveValues(default = defaultValue) 
  # REACTIVE VALUE IS UPDATED BY INPUT
  observeEvent(input$all_clusters_rows_selected,{
    print("new row selected!")
    values$default <- input$all_clusters_rows_selected # if all_clusters_rows_selected changes then update value - this sets everything!
    print(paste0("VALUE: ", values$default ))
  })
  
  # USE REACTIVE VALUE TO GENERATE ALL VARIABLES NEEDED
  mydata <- eventReactive(values$default, {
    subsetChoice <- eval(parse(text = input$navBarPage) )
    df <- subset(resultsToPlot, row.names(resultsToPlot) %in% row.names(subsetChoice))
    sel <- values$default
    print(sel)
    gene <- df[sel,]$gene
    SNP <- df[sel,]$SNP
    SNP_pos <- df[sel,]$SNP_pos
    clusterID <- row.names(df)[sel]
    print(clusterID)
    cluster_pos <- df[sel,]$cluster_pos
    # get the most significant junction in the selected cluster
    junction <- sigJunctions[variant_id == SNP & cluster_pos == paste0("chr", chr, ":", start, "-", end)]
    sigJunction <- junction$intron_cluster
    return(list(gene = gene, SNP = SNP, SNP_pos = SNP_pos, cluster_pos = cluster_pos, clusterID = clusterID, width = "auto", junction = junction))
  })
  
  # PLOTTING
  output$select_cluster_plot <- renderPlot({
    suppressWarnings( print(
      make_sQTL_cluster_plot(
                         cluster_to_plot =  mydata()$clusterID,
                         main_title = NA,
                         vcf = vcf,
                         vcf_meta = vcf_meta,
                         exons_table = exons_table,
                         counts = clusters,
                         introns = annotatedClusters,
                         cluster_ids = annotatedClusters$clusterID,
                         snp_pos = mydata()$SNP_pos,
                         snp = mydata()$SNP,
                         sigJunction = mydata()$junction)
    ))
  }, width = "auto", height = "auto",  res = 60
  )
  
  # WHOLE GENE PLOTTING
  observeEvent(values$default,{ 
    output$select_gene_plot <- renderPlot({
      suppressWarnings(print( 
        make_gene_plot(mydata()$gene,
                       counts = clusters,
                       introns = annotatedClusters,
                       exons_table = exons_table,
                       cluster_list = NULL,
                       clusterID = mydata()$clusterID,
                       introns_to_plot = introns_to_plot,
                       snp_pos = mydata()$SNP_pos,
                       snp = mydata()$SNP)
      )
      )
    }, width = mydata()$width, height = "auto", res = 90 # try changing height param
    )
  })

  # BOX PLOTS OF GENOTYPE AGAINST NORMALISED JUNCTION COUNTS
  output$select_box_plot <- renderPlot({
    #plotTitle <- c(mydata()$gene, as.character(mydata()$clusterID) )
    suppressWarnings( print(
      make_sQTL_box_plot(
        cluster_to_plot =  mydata()$clusterID,
        junction_to_plot = mydata()$junction,
        all_junctions = all_junctions,
        main_title = NA,
        vcf = vcf,
        vcf_meta = vcf_meta,
        exons_table = exons_table,
        counts = clusters,
        introns = annotatedClusters,
        cluster_ids = annotatedClusters$clusterID,
        junctionTable = junctionTable,
        snp_pos = mydata()$SNP_pos,
        snp = mydata()$SNP )
    ))
  }, width = "auto", height = "auto",  res = 90
  )
  
  # VIEW CLUSTER IN UCSC
  output$view_cluster_UCSC <- renderUI({
    coord <- mydata()$cluster_pos
    print("coord:")
    print(coord)
    snp <- mydata()$SNP
    url <- paste0("http://genome.ucsc.edu/cgi-bin/hgTracks?&org=human&db=hg19&position=", 
                  coord,"&hgFind.matches=", snp  )
  return(tags$a(href = url, "view on UCSC", target = "_blank", class = "btn btn_default", id = "UCSC" ) )
  })
  
})
