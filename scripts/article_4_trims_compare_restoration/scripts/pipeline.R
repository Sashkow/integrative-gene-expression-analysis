
my_libs = .libPaths()
my_libs
my_libs = c("/home/sashkoah/soft/rpackages", my_libs)
.libPaths(my_libs)


# BiocManager::install("SMITE")
library(plotrix)
library(erer)
library(useful)
library(gtools)
library(RCurl)
library(igraph)
library(stringr)
library(ggplot2)
library(enrichR)
library(writexl)
library(plyr)
library(SMITE)
library(hpar)
library(dplyr)


# BiocManager::install("hpar")
# remove.packages("hpar",lib =  "/home/sashkoah/R/x86_64-pc-linux-gnu-library/4.0")

source('~/a/r/igea-r/article_4_trims_compare/scripts/cluster_functional_enrichment.R')
source('~/a/r/igea-r/article_4_trims_compare/scripts/plots.R')
source('~/a/r/igea-r/article_4_trims_compare/scripts/centrality.R')
source('~/a/r/igea-r/article_4_trims_compare/scripts/gene_metadata.R')


### R code from vignette source 'STRINGdb.Rnw'

###################################################
### code chunk number 1: initialization
###################################################
# library(STRINGdb)

detach(STRINGdb)

library(easypackages)

libraries(getDependencies("STRINGdb"))

source("/home/sashkoah/a/r/igea-r/article_4_trims_compare/scripts/r_utils.R")
source("/home/sashkoah/a/r/igea-r/article_4_trims_compare/scripts/rstring.R")

# difexp_folder = '2_3_separate'

# for (confidence in seq(0, 900, by=100)){
  # print(confidence)
  # if (confidence==0){
    # confidence = 42
  # }
  # summary_pipeline(difexp_folder, confidence)
# }


difexp_folder = "1_2_no_55439"
# difexp_folder = "stas3"
confidence = 100


pipeline = function(difexp_folder, confidence){
  
  base_folder = "/home/sashkoah/a/r/igea-r/article_4_trims_compare/placenta"
  comparison = paste(difexp_folder,as.character(confidence),sep = "_")
  trim_folder = paste(difexp_folder,"/", sep = "")
  trim_folder_noslash = difexp_folder
  
  base_path = file.path(base_folder,trim_folder_noslash)
  base_path
  target_path = file.path(base_folder,comparison)

  if (!(dir.exists(target_path))){
    dir.create(target_path)
  }
  setwd(target_path)
  earlier_path = getwd()
  earlier_path
  

  comparison_summary = read.csv(file.path(base_folder, "comparison_summary.csv"))
  
  if (!(comparison %in% comparison_summary$comparison)){
    comparison_summary[nrow(comparison_summary)+1,]$comparison = comparison
  } else {
    comparion_to_delete = comparison_summary$comparison == comparison
    class(comparison_summary$comparison)
    
    comparison_summary = comparison_summary[!comparion_to_delete,]
    comparison_summary[nrow(comparison_summary)+1,]$comparison = comparison
  }
  
  cmp = comparison_summary$comparison==comparison

  
  sub_exprs = read.csv(file.path(base_path,"sub_exprs.tsv"), row.names = 1, sep='\t')
  # write.csv(pdata,file.path(base_path,"sub_pdata.csv"))
  # sub_exprs = read.csv(file.path(base_path,"sub_exprs.csv"))
  nrow(sub_exprs)
   
  comparison_summary[cmp,]$genes_in_background = nrow(sub_exprs)
  nrow(sub_exprs)
  # difexp = read.csv('/home/sashkoah/a/r/igea-r/article_3/placenta/map up down on clusters/difexpdifexp_1_2_term_global.tsv', sep = '\t')
  if ("difexp.tsv" %in% list.files(base_path)){
    difexp = read.csv(file.path(base_path,"difexp.tsv"), sep = '\t')
  } else {
    difexp = read.csv(file.path(base_path,"difexp.csv"))
    difexp = difexp[!duplicated(difexp$SYMBOL), ]
  }
  
  
  # difexp$X
  # rownames(difexp) = difexp$X
  rownames(difexp) = difexp$SYMBOL
  nrow(difexp)
  comparison_summary[cmp,]$difexp_genes = nrow(difexp)
  
  
  difexp = difexp[which(abs(difexp$logFC)>=1),]
  difexp
  
  
  comparison_summary[cmp,]$difexp_genes_logfc1 = nrow(difexp)
  nrow(difexp)
  
  
  string_db <- STRINGdb$new(version="11", species=9606, 
                             score_threshold=confidence, input_directory="")
  
  
  ###################################################
  ### code chunk number 4: map
  ###################################################
  sub_exprs
  
  if ("X" %in% rownames(sub_exprs)){
    sub_exprs$SYMBOL = sub_exprs$X  
  } else {
    sub_exprs$SYMBOL = rownames(sub_exprs)
  }
  
  
  sub_exprs_mapped <- string_db$map(sub_exprs, "SYMBOL", removeUnmappedRows = TRUE )
  
  nrow(sub_exprs)
  nrow(sub_exprs_mapped)
  
  comparison_summary[cmp,]$genes_in_background_mapped_stringdb = nrow(sub_exprs_mapped)
  nrow(sub_exprs_mapped)
  
  backgroundV = sub_exprs_mapped$STRING_id
  
  #16765
  string_db$set_background(backgroundV)
  #9839772
  #980394
  

  
  difexp_mapped <- string_db$map(difexp, "SYMBOL", removeUnmappedRows = TRUE, takeFirst = FALSE)
  nrow(difexp_mapped)
  # difexp_mapped = difexp_mapped[!duplicated(difexp_mapped$SYMBOL),]
  nrow(difexp_mapped)
  difexp_mapped = difexp_mapped[!duplicated(difexp_mapped$STRING_id ),]
  nrow(difexp_mapped)
  difexp_mapped_but_not_necessarily_on_graph = difexp_mapped
  # comparison_summary = rename(comparison_summary,c('difexp_genes_mapped'='difexp_genes_mapped_stringdb'))
  
  comparison_summary[cmp,]$difexp_genes_mapped_stringdb = nrow(difexp_mapped)
  nrow(difexp_mapped)
  
  hits <- difexp_mapped$STRING_id  
  leftout_genes = setdiff(as.character(difexp_mapped$STRING_id),as.character(difexp_mapped$STRING_id))
  leftout_genes
  
  
  # cluster_enrichment <- string_db$get_enrichment(difexp_mapped$STRING_id, category = "All")
  
  
  ###################################################
  ### code chunk number 15: clustering
  ###################################################
  
  E(string_db$graph)
  
  G = simplify(string_db$get_subnetwork(difexp_mapped$STRING_id))
  get.edgelist(G)[1]
  
  length(names(V(G)))
  
  Isolated = which(degree(G)==0)
  
  G = delete.vertices(G, Isolated)
  
  weight.community=function(row,membership,weigth.within,weight.between){
    if(as.numeric(membership[which(names(membership)==row[1])])==as.numeric(membership[which(names(membership)==row[2])])){
      weight=weigth.within
    }else{
      weight=weight.between
    }
    return(weight)
  }
  get.edgelist(G)
  E(G)$weight
  
  E(G)$weight=apply(get.edgelist(G),1,weight.community,membership(sub_fgreedy),2,1)
  E(G)$old_weight=apply(get.edgelist(sub_g),1,weight.community,membership(sub_fgreedy),1,1)
  
  
  
  
  difexp_mapped$on_graph = 1
  difexp_mapped$on_graph = ifelse(difexp_mapped$STRING_id %in% names(V(G)),1,0)
  
  
  nrow(difexp_mapped[difexp_mapped$on_graph==1,]) 
  length(names(V(G)))
  
  write.table(difexp_mapped,"difexp_with_on_graph.csv", row.names = FALSE, sep = ',')
  
  difexp_mapped = read.csv("difexp_with_on_graph.csv")
  
  nrow(difexp_mapped)
  
  
  difexp_mapped = add_cluster_column(difexp_mapped[difexp_mapped$on_graph==1,], string_db)
  # difexp_mapped = difexp
  nrow(difexp_mapped)
  table(difexp_mapped$cluster)
  
  comparison_summary[cmp,]$difexp_clustered = nrow(difexp_mapped)
  write.table(difexp_mapped,"difexp_with_clust.csv", row.names = FALSE, sep = ',')


  ###################################################
  ### code chunk number 11: enrichment and plotting
  ###################################################

  enrich = custom_get_enrichmment(difexp_mapped, earlier_path, string_db)

  comparison_summary[cmp,]$biological_process = nrow(enrich)
  difexp_mapped = boyanize_difexp(difexp_mapped, enrich)

  print(c("amount of genes, categories:", nrow(difexp_mapped),nrow(enrich)))
  
  enrich$description
  
  sub_g = G
  
  
  # clusters_2_3 = c(1,2,3,5,6,7,8)
  
  
  library(igraph)
  sub_fgreedy <- fastgreedy.community(sub_g, merges=TRUE, modularity=TRUE)
  comparison_summary[cmp,]$modularity = modularity(sub_fgreedy)
  modularity(sub_fgreedy)
  sizes(sub_fgreedy)
  


  weight.community=function(row,membership,weigth.within,weight.between){
    if(as.numeric(membership[which(names(membership)==row[1])])==as.numeric(membership[which(names(membership)==row[2])])){
      weight=weigth.within
    }else{
      weight=weight.between
    }
    return(weight)
  }

  E(sub_g)$weight=apply(get.edgelist(sub_g),1,weight.community,membership(sub_fgreedy),2,1)
  E(sub_g)$old_weight=apply(get.edgelist(sub_g),1,weight.community,membership(sub_fgreedy),1,1)

  # sub_g$layout = layout_with_fr(sub_g)
  sub_g$layout=layout.fruchterman.reingold(sub_g,weights=E(sub_g)$weight)

  sav = difexp_mapped
  difexp_mapped = add_is_local_abs_max_column(difexp_mapped)
  colnames(difexp_mapped)
  
  cluster_max_difexp = difexp_mapped[difexp_mapped$is_local_max==TRUE, c("logFC","cluster","SYMBOL")] # + "GENENAME"
  cluster_max_difexp = cluster_max_difexp[order(cluster_max_difexp$cluster),]
  cluster_max_difexp$logFC = round(cluster_max_difexp$logFC,2)
  cluster_max_difexp = cluster_max_difexp[order(cluster_max_difexp$logFC),]
  write_xlsx(cluster_max_difexp,file.path(target_path, paste(difexp_folder, "max_abs_locfc_per_cluster.xlsx", sep = "_")))
  
  # difexp_mapped$important_symbol = ifelse(difexp_mapped$is_local_max==TRUE,difexp_mapped$cluster, "")
  difexp_mapped$updown = ifelse(difexp_mapped$logFC>0,"green","red")
  
  
  
  
  #################################
  # additionally describe genes
  length(names(V(sub_g)))
  nrow(difexp_mapped)
  sav = difexp_mapped
  difexp_mapped = add_cluster_max_btw_column(difexp_mapped, sub_g)
  
  difexp_mapped = add_is_secreted_placenta_column(difexp_mapped, target_path)
  
  secreted_genes = difexp_mapped[difexp_mapped$is_secreted,]
  secteted_enrich = custom_get_enrichmment(difexp = secreted_genes, string_db = string_db, earlier_path = earlier_path)
  
  
  if (!dir.exists(file.path(target_path,"insights"))){
    dir.create(file.path(target_path,"insights"))
  }
  secreted_genes = secreted_genes[order(secreted_genes$SYMBOL),]
  write_xlsx(secreted_genes, file.path(target_path, "insights", "secreted_genes.xlsx"))
  write_xlsx(secteted_enrich, file.path(target_path, "insights", "secreted_enrich.xlsx"))
  
  # write.table(secteted_enrich[c("revigo_term","p_value")],
  #             file.path(target_path, "insights", "cluster_to_revigo.tsv"),
  #             sep = ' ', row.names = FALSE, quote = FALSE)
  

  top_btw_genes = difexp_mapped[difexp_mapped$btw >= quantile(difexp_mapped$btw, .9),]
  top_btw_genes = top_btw_genes[order(-top_btw_genes$btw),]
  top_btw_genes[,c("SYMBOL","logFC","btw", "boyan","cluster")]
  top_btw_genes = top_btw_genes[order(top_btw_genes$SYMBOL),]
  write_xlsx(top_btw_genes, file.path(target_path, "insights", "top_90_percent_btw_genes.xlsx"))
  
  # difexp_2_3 = read.csv("/home/sashkoah/a/r/igea-r/article_4_trims_compare/placenta/2_3_no_55439_100/difexp_final.csv")
  
  # difexp_mapped = difexp_2_3
  top_logfc_genes = difexp_mapped[order(-abs(difexp_mapped$logFC)),]
  top_logfc_genes = difexp_mapped[abs(difexp_mapped$logFC) >= quantile(abs(top_logfc_genes$logFC), .9),]
  top_logfc_genes = top_logfc_genes[order(top_logfc_genes$cluster),]
  top_logfc_genes[,c("SYMBOL","logFC","btw", "boyan", "cluster")]
  write_xlsx(top_logfc_genes, file.path(target_path, "insights", paste(difexp_folder, "top_90_percent_abs_logfc_genes.xlsx", sep = '_')))
  
  difexp_mapped[which(difexp_mapped$btw == max(difexp_mapped$btw)),c("SYMBOL", "logFC", "btw")]
  # sav = difexp_mapped
  difexp_mapped = sav
  difexp_mapped = add_hpa_columns(difexp_mapped)
  
  #################################
  difexp_mapped$STRING_id
  names(V(sub_g))
  difexp_mapped = difexp_mapped[match(names(V(sub_g)),difexp_mapped$STRING_id ),] # !!!
  names(V(sub_g)) == difexp_mapped$STRING_id
  difexp_mapped$STRING_id
  
  
  difexptest = difexp_mapped
  rownames(difexptest) = difexp$STRING_id
  difexptest = difexptest[names(V(sub_g)),]
  !(FALSE %in% (rownames(difexptest) == names(V(sub_g))))
  
  
  plot_the_graphs(sub_g, difexptest, sub_fgreedy)
  
  
  
  for (i in 1:50){
    print(i)
    dev.off()
  }
  
  


  if (nrow(enrich)!=0){

    l = add_category_sorting_criteria(difexp_mapped, enrich)
    enrich = l[[1]]
    difexp_mapped = l[[2]]
    #####################
    # coverage
    #####################
    enrichment_coverage = round(nrow(difexp_mapped[difexp_mapped$boyan!=0,])/nrow(difexp_mapped),2)
    print(c("enrichment_coverage", enrichment_coverage))
  }

  write.table(difexp_mapped,"difexp_with_clust_no_isolated.csv", row.names = FALSE, sep = ',')
  write.table(difexp_mapped,"difexp_final.csv", row.names = FALSE, sep = ',')

  write.table(enrich[,!names(enrich) %in% c("list_genes")],"enrichment.csv", sep = ',', row.names = FALSE)


  write.table(enrich[c("revigo_term","p_value")],"cluster_to_revigo.tsv", sep = ' ', row.names = FALSE, quote = FALSE)
  
  getwd()

  if ("average_logfc" %in% colnames(enrich)){
    write.table(enrich[c("revigo_term","average_logfc")],"logfc_to_revigo.tsv", sep = ' ', row.names = FALSE, quote = FALSE)
  }


  gene_union = NULL

  for (i in 1:nrow(enrich)){
    gene_union = union(gene_union, enrich[i,]$list_genes[[1]])
  }

  coverage = length(gene_union)/nrow(difexp_mapped)
  coverage
  comparison_summary[cmp,]$coverage = coverage
  t(comparison_summary)
  write.csv(comparison_summary, file.path(base_folder,"comparison_summary.csv"), row.names = FALSE)
  write.csv(t(comparison_summary), file.path(base_folder,"comparison_summary_t.csv"))

  
  ##################################################################3333333
  difexp_mapped = difexp_mapped[!duplicated(difexp_mapped$SYMBOL),]
  
  re = difexp_mapped
  
  g = sub_g

  nrow(re)
  V(g)
  
  
  earlier_path = target_path
  compute_summary(re, earlier_path)

  print("Computing clusters...")

  cluster_recurse(difexp = re, earlier_path = earlier_path, string_db = string_db)
  ################################################################################
  # path must have summary.csv in it
  earlier_path
  genes_covered_in_clusters = compute_genes_covered_in_clusters(earlier_path = earlier_path)
  logfc_covered_in_clusters = compute_logfc_covered_in_clusters(earlier_path = earlier_path)
  comparison_summary[cmp,]$genes_covered_in_clusters = genes_covered_in_clusters
  comparison_summary[cmp,]$cluster_coverage = genes_covered_in_clusters/comparison_summary[cmp,]$difexp_genes_mapped_stringdb
  comparison_summary[cmp,]$logfc_covered_in_clusters = logfc_covered_in_clusters
  comparison_summary[cmp,]$logfc_cluster_coverage = logfc_covered_in_clusters/sum(abs(difexp_mapped_but_not_necessarily_on_graph$logFC))

  write.csv(comparison_summary, file.path(base_folder,"comparison_summary.csv"), row.names = FALSE)
  write.csv(t(comparison_summary), file.path(base_folder,"comparison_summary_t.csv"))
}



# difexp_folder = "stas"
confidence = 300

# pipeline(difexp_folder, confidence)
# to_revigo.csv
earlier_path = "/home/sashkoah/a/r/igea-r/article_4_trims_compare/placenta/stas_norm_preeclampsia_1_study_100"
enrichment_after_revigo(earlier_path = earlier_path)


clusters_folder = "/home/sashkoah/a/r/igea-r/article_4_trims_compare/placenta/stas_norm_preeclampsia_1_study_100/clusters"
# emergency(clusters_folder)

create_xlsx_files(earlier_path = earlier_path)


for (confidence in seq(200, 900, by=100)){
  print(confidence)
  if (confidence==0){
    confidence = 42
  }
  pipeline(difexp_folder, confidence)
}


 # summary_pipeline(difexp_folder, 250)




base_folder = "/home/sashkoah/a/r/igea-r/article_3/placenta"
comparison = paste(difexp_folder,as.character(confidence),sep = "_")
trim_folder = paste(difexp_folder,"/", sep = "")
trim_folder_noslash = difexp_folder


target_path = file.path(base_folder,comparison)

if (!(dir.exists(target_path))){
  dir.create(target_path)
}
setwd(target_path)
earlier_path = getwd()
earlier_path



re = difexp_mapped
g = sub_g

nrow(re)
V(g)

# earlier_path = file.path(base_folder,trim_folder_noslash,"clusters", current_cluster)
earlier_path

compute_summary(re, earlier_path)
#compute united_summary

# cluster_recurse(difexp = re, earlier_path = earlier_path)

emergency = function(clusters_folder){
  paths = list.dirs(clusters_folder, recursive = FALSE)
  cluster_enrichment_combined = data.frame()
  cluster_didfexp_combined = data.frame()
  path = paths[1]
  for (path in paths){
    split = strsplit(path,'/')[[1]]
    cluster_number = split[length(split)]
    if (file.exists(file.path(path,"cluster_enrichment_united.csv"))){
      cluster_enrichment = read.csv(file.path(path,"cluster_enrichment_united.csv"))
      cluster_enrichment = cluster_enrichment[which(cluster_enrichment$eliminated==0),]
      if (nrow(cluster_enrichment)>0){
        cluster_enrichment$cluster = cluster_number
        cluster_enrichment_combined = rbind(cluster_enrichment_combined, cluster_enrichment)
      }
    }
    
    if (file.exists(file.path(path,"cluster_difexp.csv"))){
      cluster_enrichment = read.csv(file.path(path,"cluster_difexp.csv"))
      cluster_enrichment = cluster_enrichment[which(cluster_enrichment$eliminated==0),]
      if (nrow(cluster_enrichment)>0){
        cluster_enrichment$cluster = cluster_number
        cluster_enrichment_combined = rbind(cluster_enrichment_combined, cluster_enrichment)
      }
    }
    
    
  }
  cluster_enrichment_combined
  # as.data.frame(table(cluster_enrichment_combined[,c("cluster")]))
  write.csv(cluster_enrichment_combined,file.path(clusters_folder,"cluster_enrichment_combined.csv"), row.names = FALSE)
}

clusters_folder = "/home/sashkoah/a/r/igea-r/article_3/placenta/stas_200/clusters"
emergency(clusters_folder)


difexp_folder = '1_2_separate_new'
confidence = 700





summary_pipeline = function(difexp_folder, confidence){
  base_folder = "/home/sashkoah/a/r/igea-r/article_3/placenta"
  comparison = paste(difexp_folder,as.character(confidence),sep = "_")
  trim_folder = paste(difexp_folder,"/", sep = "")
  trim_folder_noslash = difexp_folder
  
  base_path = file.path(base_folder,trim_folder)
  target_path = file.path(base_folder,comparison)
  
  if (!(dir.exists(target_path))){
    dir.create(target_path)
  }
  setwd(target_path)
  earlier_path = getwd()
  earlier_path
  
  comparison_summary = read.csv(file.path(base_folder, "comparison_summary.csv"))
  
  if (!(comparison %in% comparison_summary$comparison)){
    comparison_summary[nrow(comparison_summary)+1,]$comparison = comparison
  } else {
    # comparion_to_delete = comparison_summary$comparison == comparison
    class(comparison_summary$comparison)
    
    # comparison_summary = comparison_summary[!comparion_to_delete,]
    # comparison_summary[nrow(comparison_summary)+1,]$comparison = comparison
  }
  
  cmp = comparison_summary$comparison==comparison
  
  if (!("cluster_coverage" %in% colnames(comparison_summary))){
    comparison_summary$cluster_coverage = NA
  }
  
  if (!("cluster_count" %in% colnames(comparison_summary))){
    comparison_summary$cluster_count= NA
  }
  
  if (!("coverage_including_isolated" %in% colnames(comparison_summary))){
    comparison_summary$coverage_including_isolated= NA
  }
  
  
  # cluster_coverage = compute_cluster_coverage(target_path)
  # comparison_summary[cmp,]$cluster_coverage = cluster_coverage
  if (file.exists(file.path(target_path,"summary.csv"))){
    summary = read.csv(file.path(target_path,"summary.csv"))
    cluster_count = nrow(summary)
  } else {
    cluster_count = 0
  }
  
  comparison_summary[cmp,]$cluster_count = cluster_count
  
  
  enrich = read.csv(file.path(target_path, "enrichment.csv"))
  enrich$list_genes = strsplit(enrich$inputGenes,",")
  
  for (i in 1:nrow(enrich)){
    gene_union = union(gene_union, enrich[i,]$list_genes[[1]])
  }
  
  difexp_with_isolated = read.csv(file.path(target_path,"difexp_with_isolated.csv"))
  enrich = custom_get_enrichmment(difexp_with_isolated, target_path)
  
  write.csv(comparison_summary, file.path(target_path,"enrichment_with_isolated.csv"), row.names = FALSE)
  
  
  
  coverage_including_isolated = length(gene_union)/comparison_summary[cmp,]$difexp_genes_mapped_stringdb
  comparison_summary[cmp,]$coverage_including_isolated = coverage_including_isolated
  
  cluster_coverage = sum(summary$genes_total * summary$enrichment_coverage) / comparison_summary[cmp,]$difexp_genes_mapped_stringdb
  comparison_summary[cmp,]$cluster_coverage = cluster_coverage
  
  print("here")
  
  print(comparison_summary[cmp,])
  
  write.csv(comparison_summary, file.path(base_folder,"comparison_summary.csv"), row.names = FALSE)
  write.csv(t(comparison_summary), file.path(base_folder,"comparison_summary_t.csv"))
  Sys.sleep(1)
}




