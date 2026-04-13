
library(writexl)
# difexp = difexp_mapped[difexp_mapped$on_graph==1,]
add_cluster_column = function(difexp, string_db){
  #remove cluster column
  difexp = difexp[,!(names(difexp) %in% c("cluster"))]
  the_g = simplify(string_db$get_subnetwork(difexp$STRING_id))
  V(the_g)
  fgreedy<-fastgreedy.community(the_g, merges=TRUE, modularity=TRUE)
  com <- cbind(V(the_g)$name,membership(fgreedy)) 
  com = as.data.frame(com)
  colnames(com)[colnames(com) == 'V2'] <- 'cluster'
  difexp = merge(difexp,com,by.x = "STRING_id", by.y = "V1")
  difexp$cluster = as.numeric(difexp$cluster)
  return(difexp)
}

# difexp =read.csv(file.path(target_path,"difexp_with_clust_no_isolated.csv"), header = TRUE)

compute_summary <- function(difexp, earlier_path){
  summary = as.data.frame(table(difexp[,"cluster"]))
  # summary$Var1 = as.numeric(as.character(summary$Var1))
  # summary$Freq = as.numeric(as.character(summary$Freq))
  # summary = summary[order(summary$Var1),]
  # row.names(summary) <- NULL
  
  colnames(summary) = c("cluster", "genes_total")
  summary$cluster = as.numeric(as.character(summary$cluster))
  updowntable = table(difexp[,c("cluster", "updown")])
  summary$up = 0
  summary$down = 0
  # summary$sum_positive_logfc = 0
  # summary$sum_negative_logfc = 0
  summary
  
  get_measure_logfc <- function(x,measure) round(measure(difexp[difexp$cluster==x['cluster'],]$logFC),2)
  
  for ( i in 1:nrow(summary)){
    summary[i,]$up = updowntable[i]
    summary[i,]$down = updowntable[i+nrow(summary)]
    # summary[i,]$sum_positive_logfc = round(sum(difexp[which(difexp$logFC>0 & difexp$cluster==i),]$logFC),2)
    # summary[i,]$sum_negative_logfc = round(sum(difexp[which(difexp$logFC<0 & difexp$cluster==i),]$logFC),2)
  }
  summary
  
  summary$logFC_min = apply(summary, 1, get_measure_logfc, measure=min)
  summary$logFC_max = apply(summary, 1, get_measure_logfc, measure=max)
  summary$logFC_mean = apply(summary, 1, get_measure_logfc, measure=mean)
  meanabs <- function(x) mean(abs(x))
  summary$logFC_mean_abs = apply(summary, 1, get_measure_logfc, measure=meanabs)
  rbind(summary)
  summary$cluster = as.numeric(summary$cluster)
  nrow(summary)
  # total = c("total", 
  #           nrow(difexp),
  #           nrow(difexp[difexp$updown=="green",]),
  #           nrow(difexp[difexp$updown=="red",]),
  #           # round(sum(difexp[which(difexp$logFC>0),]$logFC)/nrow(difexp[which(difexp$logFC>0),]),2),
  #           # round(sum(difexp[which(difexp$logFC<0),]$logFC)/nrow(difexp[which(difexp$logFC<0),],2),
  #           # round(mean(difexp$logFC),2),
  #           # round(mean(abs(difexp$logFC)),2)
  #           )
  summary = summary[order(summary$cluster),]
  write.csv(summary,file.path(earlier_path,"summary.csv"), quote = FALSE, row.names = FALSE)
} 




custom_get_enrichmment = function(difexp, earlier_path, string_db, max_pvalue=.05){
  # print(c("enriching", nrow(difexp), "genes", ",", min(abs(difexp$logFC)),max(abs(difexp$logFC)) ))
  
  cluster = difexp$STRING_id
  
  cluster_enrichment_all <- string_db$get_enrichment(cluster)
  write.table(cluster_enrichment_all, file.path(earlier_path,"cluster_enrichment_all.csv"), sep = ',', row.names = FALSE)
  cluster_enrichment = cluster_enrichment_all[cluster_enrichment_all$category=="Process",]
  
  if (nrow(cluster_enrichment)==0){
    return(cluster_enrichment)
  }
  #filter by pvalue
  cluster_enrichment = cluster_enrichment[cluster_enrichment$p_value<max_pvalue,]
  
  cluster_enrichment$revigo_term = str_replace(cluster_enrichment$term,"\\.",":")
  
  if (grepl(",",cluster_enrichment$inputGenes[1])){
    cluster_enrichment$list_genes = strsplit(cluster_enrichment$inputGenes,",")
  } else {
    cluster_enrichment$list_genes = strsplit(cluster_enrichment$inputGenes,";")
  }
  
  return(cluster_enrichment[cluster_enrichment$category=="Process",])
}


boyan = function(gene, enrich){
  count = 0
  for (i in 1:nrow(enrich)){
    if (gene %in% enrich[i,]$list_genes[[1]]){
      count = count + 1
    }
  }
  return(count)
}


# difexp = cluster_difexp
# enrich = upper_categories


# difexp = cluster_difexp
# enrich = cluster_enrichment
boyanize_difexp = function(difexp, enrich){
  difexp$boyan = 0
  if (nrow(enrich)==0){
    return(difexp)
  }
  enrich$list_genes = strsplit(enrich$inputGenes, ",")
  for (gene in difexp$STRING_id){
    gene
    difexp[difexp$STRING_id==gene,]$boyan = boyan(gene, enrich)
  }
  return(difexp)
}


# difexp = re[re$cluster==2,]
# nrow(difexp)
# enrich <- custom_get_enrichmment(difexp)
# nrow(enrich)
# invisible_cluster = difexp
# anno$string_ids

# anno$list_genes = strsplit(anno$string_ids,",")
# invisible_cluster_b = boyanize_difexp(invisible_cluster, anno)
# invisible_cluster_b$boyan
# l = add_category_sorting_criteria(invisible_cluster_b, anno)
# sorted_anno = l[[1]]
# sorted_anno


# write.csv(sorted_anno[,!names(sorted_anno) %in% c("list_genes")], "soted_anno_invisible_cluster_2.csv", row.names = FALSE)

add_category_sorting_criteria = function(difexp, enrich){
  enrich$boyan = 0
  enrich$boyan_sum = 0
  enrich$discounted_genes_amount = 0
  enrich$uniqueness_and_logfc = 0
  enrich$average_logfc = 0
  enrich$sum_logfc = 0
  
  row = 1
  for (row in 1:nrow(enrich))
  {
    boyan_sum = 0
    discounted_genes_amount = 0
    uniqueness_and_logfc = 0 
    just_logfc = 0
    for (gene in enrich[row,]$list_genes[[1]]){
      # print(c(gene, difexp[difexp$STRING_id==gene,]$boyan))
      # print(c("here",gene,row))
      # print(re[re$STRING_id==gene,]$boyan)
      boyan_sum = boyan_sum + difexp[difexp$STRING_id==gene,]$boyan
      discounted_genes_amount = discounted_genes_amount + 1/difexp[difexp$STRING_id==gene,]$boyan
      uniqueness_and_logfc = uniqueness_and_logfc + (1/difexp[difexp$STRING_id==gene,]$boyan) * abs(difexp[difexp$STRING_id==gene,]$logFC)
      just_logfc = just_logfc + abs(difexp[difexp$STRING_id==gene,]$logFC)
    }
    
    
    enrich[row,]$boyan_sum = boyan_sum
    enrich[row,]$boyan = boyan_sum/length(enrich[row,]$list_genes[[1]])
    enrich[row,]$discounted_genes_amount = discounted_genes_amount
    enrich[row,]$uniqueness_and_logfc = uniqueness_and_logfc
    enrich[row,]$average_logfc = just_logfc/length(enrich[row,]$list_genes[[1]])
    enrich[row,]$sum_logfc = just_logfc
    
    # print(c(boyan_sum,enrich[row,]$boyan))
  }
  
  return(list(enrich,difexp))
}



########################################################################################################################################

# difexp = re
# earlier_path = earlier_path
# string_db = string_db
# string_db
# earlier_path = "/home/sashkoah/a/r/igea-r/article_3/placenta/1_2_no_55439_100"
# list.files(earlier_path)
# difexp = read.csv(file.path(earlier_path, "difexp_with_clust_no_isolated.csv"))

cluster_recurse = function(difexp, earlier_path, string_db){
  
  # re is difexp table with cluster and string_id mappings
  setwd(earlier_path)
  
  # summary = read.csv("summary.csv")
  
  difexp = difexp[,!(names(difexp) %in% c("cluster"))]
  the_g = simplify(string_db$get_subnetwork(difexp$STRING_id))
  fgreedy<-fastgreedy.community(the_g, merges=TRUE, modularity=TRUE)
  com <- cbind(V(the_g)$name,membership(fgreedy)) 
  com = as.data.frame(com)
  colnames(com)[colnames(com) == 'V2'] <- 'cluster'
  difexp = merge(difexp,com,by.x = "STRING_id", by.y = "V1")
  getwd()
  if (!dir.exists("clusters")){
    dir.create("clusters")
  }
  
  
  summary = read.csv("summary.csv",header = TRUE)
  # summary = read.csv("united_summary.csv",header = TRUE)
  summary$enrichment_coverage = 0
  # sum abs logfc in genes covered by enrichment
  summary$logfc_covered = 0
  
  
  cluster_number = 1
  for (cluster_number in 1:length(fgreedy)){ 
    p = file.path(getwd(),"clusters", as.character(cluster_number))
    
    if (!dir.exists(p)){
      dir.create(p)
    }
    
    setwd(as.character(p))
    print(getwd())
    cluster_difexp = difexp[difexp$cluster==cluster_number,]
    
    cluster = cluster_difexp$STRING_id
    
    print(c("cluster number, length", cluster_number, length(cluster)))
    
    print("plotting graphs")
    sub_g = simplify(string_db$get_subnetwork(cluster))
    sub_g$layout=layout_with_fr(sub_g)
    sub_fgreedy<-fastgreedy.community(sub_g, merges=TRUE, modularity=TRUE)
    
    cluster_enrichment = custom_get_enrichmment(cluster_difexp, earlier_path, string_db)
    nrow(cluster_enrichment)
    
    cluster_difexp = boyanize_difexp(cluster_difexp, cluster_enrichment)
    cluster_difexp$boyan
    print(c("amount of genes, categories:", nrow(cluster_difexp),nrow(cluster_enrichment)))
    #################add_cluster_column
  
    cluster_difexp = add_cluster_column(cluster_difexp, string_db)
    
    
    
    plot_the_graphs(sub_g, cluster_difexp, sub_fgreedy)
    
    
    
    if (nrow(cluster_enrichment)!=0){
      l = add_category_sorting_criteria(cluster_difexp, cluster_enrichment)
      cluster_enrichment = l[[1]]
      cluster_difexp = l[[2]]

      #####################
      # coverage
      #####################
      gene_union = NULL
      
      enrichment_coverage = round(nrow(cluster_difexp[cluster_difexp$boyan!=0,])/nrow(cluster_difexp),2)
      
      logfc_covered = sum(abs(cluster_difexp[cluster_difexp$boyan!=0,]$logFC))
      
      
      summary[cluster_number,]$enrichment_coverage = enrichment_coverage
      summary[cluster_number,]$logfc_covered = logfc_covered
    } else {
      summary[cluster_number,]$enrichment_coverage = 0
      summary[cluster_number,]$logfc_covered = 0
    }
    print(c("coverage", summary[cluster_number,]$enrichment_coverage))
    write.table(cluster_difexp, "cluster_difexp.csv", row.names = FALSE, sep = ',')
    write.table(cluster_enrichment[,!names(cluster_enrichment) %in% c("list_genes")],"cluster_enrichment.csv", sep = ',', row.names = FALSE)
    if (("revigo_term" %in% colnames(cluster_enrichment)) &&  ("p_value" %in% colnames(cluster_enrichment)) ){
      write.table(cluster_enrichment[,c("revigo_term","p_value")],"cluster_to_revigo.tsv", sep = ' ', row.names = FALSE, quote = FALSE)
    }
    
    setwd(earlier_path)
  }
  # write.table(summary, file.path(earlier_path,"united_summary.csv"), sep = ',', row.names = FALSE)
  write.table(summary, file.path(earlier_path,"summary.csv"), sep = ',', row.names = FALSE)
}


list_enrichment_in_difexp = function(difexp, enrichment){
  enrichment$list_genes = strsplit(enrichment$inputGenes,",")
  difexp$categories = "_"
  for (i in 1:nrow(difexp)){
    gene = difexp[i,]$STRING_id
    categories = NA
    for (j in 1:nrow(enrichment)){
      if (gene %in% enrichment[j,]$list_genes[[1]]){
        if (is.na(categories)){
          categories = enrichment[j,]$description
        } else {
          categories = paste(categories,enrichment[j,]$description ,sep=", ")
        }
      }
    }
    difexp[i,]$categories = categories
  }
  return(difexp)
}


compute_genes_covered_in_clusters <- function(earlier_path) {
  summary = read.csv(file.path(earlier_path,"summary.csv"))
  genes_covered_in_clusters = sum(summary$genes_total * summary$enrichment_coverage)
  return(genes_covered_in_clusters)
}

compute_logfc_covered_in_clusters <- function(earlier_path) {
  summary = read.csv(file.path(earlier_path,"summary.csv"))
  logfc_covered_in_clusters = sum(summary$logfc_covered)
  return(logfc_covered_in_clusters)
}

reunite_cluster_difexp = function(earlier_path){
  
}

plot_category_logfc_distribution_hists = function(path, cluster_difexp, upper_categories){
  # use with enrichment_after_revigo
  #path - cluster path 
  unlink(file.path(path,"hists"), recursive=TRUE)
  dir.create(file.path(path,"hists"))
  
  logfcmin = min(cluster_difexp$logFC) - 1
  logfcmax = max(cluster_difexp$logFC) + 1
  breaks = round((ceiling(logfcmax)-floor(logfcmin))/.5)
  print(c(floor(logfcmin), ceiling(logfcmax)))
  for (i in 1:nrow(upper_categories)){
    # print(upper_categories[i,])
    upper_category_genes = strsplit(upper_categories[i,]$inputGenes,",")[[1]]
    logfcs = cluster_difexp[which(cluster_difexp$STRING_id %in% upper_category_genes),]$logFC
    filename = paste(str_replace_all(upper_categories[i,]$description," ","_"), ".png",sep = "")
    filename = str_replace_all(filename,"/","__")
    png(file.path(path,"hists",filename), height = 800, width = 600)
    hist(logfcs, main=upper_categories[i,]$description, breaks=breaks, xlim=c(floor(logfcmin), ceiling(logfcmax)))
    dev.off()
  }
}




# re is difexp table with cluster and string_id mappings
# g is igraph of a difexp graph



# re = difexp_mapped
# g = sub_g
# 
# 
# nrow(re)
# V(g)
# 
# earlier_path = file.path(base_folder,trim_folder_noslash,"clusters", current_cluster)
# earlier_path = "/home/sashkoah/a/r/igea-r/article_4_trims_compare/placenta/stas2_100/clusters/5"
# earlier_path
#
# difexp = read.csv(file.path(earlier_path, "cluster_difexp.csv"))
#
# difexp = add_cluster_column(difexp, string_db)
#
# compute_summary(difexp, earlier_path)
#
# cluster_recurse(difexp, earlier_path, string_db)
# 
# # to_revigo.py
# 
# enrichment_after_revigo(earlier_path = earlier_path)
# 
# create_xlsx_files(earlier_path = earlier_path)





# merge enrichment with revigo results, plot hists
enrichment_after_revigo = function(earlier_path){
  paths = list.dirs(file.path(earlier_path, "clusters"), recursive = FALSE)
  summary = read.csv(file.path(earlier_path,"summary.csv"),header = TRUE)
  summary$coverage_revigo = 0
  summary$genes_covered_revigo = 0
  summary$logfc_covered_revigo = 0
  summary$logfc_sum = 0
  # path = paths[4]
  # path
  for (path in paths ){
    cluster_number = tail(strsplit(path,'/')[[1]], n=1)
    cluster_number
    if (file.exists(file.path(path,"revigo.csv"))){
      revigo = read.csv(file.path(path,"revigo.csv"))
      cluster_difexp = read.csv(file.path(path,"cluster_difexp.csv"))
      revigo = revigo[,c("TermID","Uniqueness","Dispensability","Representative","Eliminated")]
      names(revigo)[names(revigo) == "TermID"] <- "term_ID"
      names(revigo)[names(revigo) == "Uniqueness"] <- "uniqueness"
      names(revigo)[names(revigo) == "Dispensability"] <- "dispensability"
      names(revigo)[names(revigo) == "Representative"] <- "representative"
      names(revigo)[names(revigo) == "Eliminated"] <- "eliminated"
      revigo = revigo %>% 
        mutate(across(where(is.character), str_trim))
      
      if (file.exists(file.path(path,"cluster_enrichment.csv"))){
        cluster_enrichment = read.csv(file.path(path,"cluster_enrichment.csv"))
      } else {
        next
      }
      
      cluster_difexp = read.csv(file.path(path,"cluster_difexp.csv"))
      # print(c(nrow(cluster_difexp), nrow(cluster_enrichment),nrow(cluster_enrichment_m)))
      cluster_enrichment_m = merge(cluster_enrichment,revigo, by.x = "revigo_term", by.y = "term_ID")
      cluster_enrichment_m[cluster_enrichment_m$eliminated=="False",]
      cluster_difexp = list_enrichment_in_difexp(cluster_difexp, cluster_enrichment_m[cluster_enrichment_m$eliminated=="False",])
      revigo_cluster_enrichmnet = cluster_enrichment_m[cluster_enrichment_m$eliminated==0,]
      
      
      write.csv(cluster_difexp, file.path(path,"cluster_difexp.csv"), row.names = FALSE)
      write.csv(cluster_enrichment_m, file.path(path,"cluster_enrichment_united.csv"), row.names = FALSE)
      upper_categories = cluster_enrichment_m[cluster_enrichment_m$eliminated==0,]
      upper_categories = upper_categories[order(-upper_categories$number_of_genes),]
      
      cluster_difexp = boyanize_difexp(cluster_difexp, upper_categories)
      genes_covered_after_revigo = nrow(cluster_difexp[which(cluster_difexp$boyan!=0),])
      coverage = genes_covered_after_revigo/nrow(cluster_difexp)
      
      summary[which(summary$cluster==cluster_number),]$coverage_revigo = coverage
      summary[which(summary$cluster==cluster_number),]$genes_covered_revigo = genes_covered_after_revigo
      summary[which(summary$cluster==cluster_number),]$logfc_covered_revigo = sum(abs(cluster_difexp[which(cluster_difexp$boyan!=0),]$logFC))
      summary[which(summary$cluster==cluster_number),]$logfc_sum = sum(abs(cluster_difexp$logFC))
    }
  }
  write.csv(summary,file.path(earlier_path,"summary.csv"), quote = FALSE, row.names = FALSE)
}

# earlier_path = "/home/sashkoah/a/r/igea-r/article_3/placenta/1_2_no_55439_100/clusters/2"
# enrichment_after_revigo(earlier_path = earlier_path)



create_xlsx_files = function(earlier_path){
  # list.files(earlier_path)
  folder_name = tail(strsplit(earlier_path,'/')[[1]],1)
  
  if (file.exists(file.path(earlier_path,"summary.csv"))){
    file_name = paste(folder_name, "summary.xlsx", sep = "_")
    file = read.csv(file.path(earlier_path,"summary.csv"))
    write_xlsx(file, file.path(earlier_path,file_name))
  }
  
  if (file.exists(file.path(earlier_path,"clusters/cluster_enrichment_combined.csv"))){
    file_name = paste(folder_name, "cluster_enrichment_combined.xlsx", sep = "_")
    file = read.csv(file.path(earlier_path,"clusters/cluster_enrichment_combined.csv"))
    write_xlsx(file, file.path(earlier_path,file_name))
  }
  
  if (file.exists(file.path(earlier_path,"difexp_with_on_graph.csv"))){
    file_name = paste(folder_name, "difexp.xlsx", sep = "_")
    file = read.csv(file.path(earlier_path,"difexp_with_on_graph.csv"))
    write_xlsx(file, file.path(earlier_path,file_name))
  }
  
  if (file.exists(file.path(earlier_path,"enrichment.csv"))){
    file_name = paste(folder_name, "enrichment.xlsx", sep = "_")
    file = read.csv(file.path(earlier_path,"enrichment.csv"))
    write_xlsx(file, file.path(earlier_path,file_name))
  }
  
  paths = list.dirs(file.path(earlier_path, "clusters"), recursive = FALSE)
  # create short files
  path = paths[3]
  for (path in paths){
    cluster_name = tail(strsplit(path,'/')[[1]],1)
    pages = list()
    file.exists(file.path(path,"cluster_enrichment_united.csv"))
    if (file.exists(file.path(path,"cluster_enrichment_united.csv"))){
      print(path)
      cluster_enrichment_united = read.csv(file.path(path,"cluster_enrichment_united.csv"))
      cluster_enrichment_united$p_group = 1
      
      
      
      
      
      to_int = function(row) as.character(as.numeric(strsplit(row['term'], ".", fixed=TRUE)[[1]][2]))
      cluster_enrichment_united$int_term=apply(cluster_enrichment_united,1,to_int)
      i=1
      for (i in 1:nrow(cluster_enrichment_united)){
        if (cluster_enrichment_united[i,]$eliminated=="False"){
          cluster_enrichment_united[i,]$p_group = cluster_enrichment_united[i,]$p_value
        } else{
          cluster_enrichment_united[i,]$p_group = cluster_enrichment_united[which(cluster_enrichment_united$int_term==cluster_enrichment_united[i,]$representative & cluster_enrichment_united$eliminated=="False"),]$p_value
        }
      }
      
      # cluster_enrichment_united[order(cluster_enrichment_united[,"representative"]),c("p_group", "representative","p_value", "eliminated")]
      odd_columns = c("boyan", "boyan_sum", "discounted_genes_amount", "uniqueness_and_logfc")
      filename = paste(cluster_name,"enrich", sep = '_')
      pages[[filename]] = cluster_enrichment_united[,!names(cluster_enrichment_united) %in% odd_columns]
      # write_xlsx(pages, file.path(path,filename))
    }
    
    if (file.exists(file.path(path,"cluster_difexp.csv"))){
      print(path)
      cluster_difexp = read.csv(file.path(path,"cluster_difexp.csv"))
      
      odd_columns = c("AveExpr", "t",	"P.Value",	"adj.P.Val",	"B",	"First.Trimester.Average",
                      "Second.Trimester.Average", "isolated", "is_local_max","important_symbol",
                      "updown", "covered")
      short_difexp = cluster_difexp[,!names(cluster_difexp) %in% odd_columns]
      filename = paste(cluster_name, "difexp", sep = '_')
      pages[[filename]] = cluster_difexp[,!names(cluster_difexp) %in% odd_columns] 
    }
    
    save_path = file.path(path,paste(folder_name, cluster_name, ".xlsx", sep = '_'))
    names(pages)
    write_xlsx(pages, save_path)
  }
}

# earlier_path = "/home/sashkoah/a/r/igea-r/article_3/placenta/1_2_no_55439_100/clusters/1"
# create_xlsx_files(earlier_path = earlier_path)

