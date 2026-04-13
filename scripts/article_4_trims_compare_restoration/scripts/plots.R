# g = sub_g
# difexp = difexp_mapped
# plot_the_graphs(g,difexp)

# difexp_folder = "1_2_no_55439"
# confidence = 100
# base_folder = "/home/sashkoah/a/r/igea-r/article_3/placenta"
# comparison = paste(difexp_folder,as.character(confidence),sep = "_")
# trim_folder = paste(difexp_folder,"/", sep = "")
# trim_folder_noslash = difexp_folder
# base_path = file.path(base_folder,trim_folder_noslash)
# target_path = file.path(base_folder,comparison)

get_cluster = function(g, v){
  # get cluster of vertex v in graph g
  clusters_df = as.data.frame(V(g)$community)
  clusters_df$vertex_id = names(V(g))
  clusters_df[clusters_df$vertex_id==v,]$`V(g)$community`
}



# get_difexp = function(path){
#   filepath = file.path(path,"difexp_with_clust_no_isolated.csv")
#   if (file.exists(filepath)){
#     return(read.csv(filepath))
#   }
#   filepath = file.path(path,"cluster_difexp.csv")
#   if (file.exists(filepath)){
#     return(read.csv(filepath))
#   }
# }
# 
# 
# 
# map_data_to_string = function(difexp_folder, confidence){
#   base_folder = "/home/sashkoah/a/r/igea-r/article_3/placenta"
#   comparison = paste(difexp_folder,as.character(confidence),sep = "_")
#   trim_folder = paste(difexp_folder,"/", sep = "")
#   trim_folder_noslash = difexp_folder
#   base_path = file.path(base_folder,trim_folder_noslash)
#   target_path = file.path(base_folder,comparison)
#   
#   if (!(dir.exists(target_path))){
#     dir.create(target_path)
#   }
#   setwd(target_path)
#   earlier_path = getwd()
#   sub_exprs = read.csv(file.path(base_path,"sub_exprs.tsv"), row.names = 1, sep='\t')
#   string_db <- STRINGdb$new(version="11", species=9606, 
#                             score_threshold=confidence, input_directory="")
#   
#   sub_exprs$SYMBOL = row.names(sub_exprs)
#   sub_exprs_mapped <- string_db$map(sub_exprs, "SYMBOL", removeUnmappedRows = TRUE )
#   backgroundV = sub_exprs_mapped$STRING_id
#   string_db$set_background(backgroundV)
#   setwd(base_folder)
#   return(string_db)
# }
# 
# 
# get_g = function(path, difexp){
#   string_db = map_data_to_string(difexp_folder, confidence)
#   g = simplify(string_db$get_subnetwork(difexp_mapped$STRING_id))
#   return(g)
# }
# 
# difexp = get_difexp(target_path)
# g = get_g(target_path, difexp)
# plot_the_graphs(g, difexp, target_path)



# g = sub_g
# difexp = difexptest
# difexp = difexp_mapped
#
# length(V(g))
# nrow(difexp)
# difexp = difexp[!duplicated(difexp$STRING_id),]


# sub_fgreedy = sub_fgreedy_old


plot_the_graphs = function(g, difexp, sub_fgreedy){
  difexp = difexp[!duplicated(difexp$STRING_id),]
  # sub_fgreedy = fastgreedy.community(g, merges=TRUE, modularity=TRUE)
  # g$layout=layout.fruchterman.reingold(g,weights=E(g)$weight)
  ############################
  # count local max
  # difexp = add_is_local_abs_max_column(difexp)
  # difexp$is_local_max = FALSE
  # max(difexp$cluster)
  # for (i in 1:max(difexp$cluster)){
  #   cluster_max = max(abs(difexp[difexp$cluster == i,]$logFC))
  #   max_gene = difexp[abs(difexp$logFC) == cluster_max,]$STRING_id[1]
  #   difexp[difexp$STRING_id==max_gene,]$is_local_max = TRUE
  # }
  
  
  cluster_max_difexp = difexp[difexp$is_local_max==TRUE ,c("logFC","cluster","SYMBOL")] #,"GENENAME"
  cluster_max_difexp = cluster_max_difexp[order(cluster_max_difexp$cluster),]
  cluster_max_difexp$logFC = round(cluster_max_difexp$logFC,2)
  write_xlsx(cluster_max_difexp,file.path(target_path, paste(difexp_folder,"max_abs_locfc_per_cluster.xlsx", sep = "_")))
  
  difexp$important_symbol = ifelse(difexp$is_local_max==TRUE,difexp$cluster, "")
  difexp$updown = ifelse(difexp$logFC>0,"green","red")
  ###################################
  legend_dot_size = 0.025
  svg_width = 20
  svg_height = 10
  scale=1
  V(g)$community <- sub_fgreedy$membership
  nb.cols <- length(unique(V(g)$community))
  mycolors <- colorRampPalette(brewer.pal(8,"Dark2"))(nb.cols)
  colrs = adjustcolor(mycolors)
  
  table(difexp$max_tissue)
  
  print('plotting secreted graph')
  
  #################################################
  
  svg(file.path(getwd(), "tissue_specificity_graph.svg"), width=svg_width, height=svg_height)
  par(mar=c(0,0,0,0))
  
  # print(c(cluster_difexp$important_symbol,cluster_difexp$updown))
  current_g = g
  current_g$layout = norm_coords(current_g$layout,xmin = -scale*2,xmax = scale*2,ymin = -scale,ymax = scale)
  
  # data(hpaSecretome)
  # difexpm = merge(difexp,hpaSecretome,by.x="SYMBOL", by.y="Gene.name", all.x=TRUE)
  # nrow(difexpm)
  # difexp = difexpm
  property = "max_tissue"
  
  difexp[,property] = ifelse(is.na(difexp[,property]), "Unknown", difexp[,property])
  tissue_specificity_colors <- colorRampPalette(brewer.pal(8,"Set3"))(length(levels(as.factor(difexp[,property]))))
  tissue_specificity_colors = adjustcolor(tissue_specificity_colors)
  
  ##colors legend
  tissue_specificity_colors
  # labels legend
  levels(as.factor(difexp[,property]))
  
  
  
  V(current_g)$color <- colrs[V(current_g)$community]
  E(current_g)$color <- apply(as.data.frame(get.edgelist(current_g)), 1, 
                              function(x) ifelse(get_cluster(current_g,x[1]) == get_cluster(current_g,x[2]), 
                                                 colrs[get_cluster(current_g,x[1])], '#AAAAAAAA'))
  
  plot(current_g, vertex.label.dist=0,
       vertex.label.color= "black",
       vertex.label.font=2,
       # vertex.label.cex = 3,
       vertex.label=NA,
       vertex.color=colrs[V(current_g)$community],
       vertex.size=abs(difexp$logFC),
       edge.color=E(current_g)$color,
       rescale = FALSE)
  
  
  # a = legend('bottomleft', legend=levels(as.factor(difexp[,property])))
  # x <- (a$text$x + a$rect$left) / 2
  # y <- a$text$y
  # symbols(x,y,
  #         circles=rep(legend_dot_size,length(tissue_specificity_colors)),
  #         inches=FALSE,add=TRUE,bg=tissue_specificity_colors)
  
  df_layout = as.data.frame(current_g$layout)
  
    
  addShadowText(x = df_layout$V1,
                y = df_layout$V2,
                ifelse(difexp$max_tissue!="Uknown",difexp$max_tissue,""),
                col="black",
                bg="white",
                cex = .5)
  
 
  
 
    
  
  
  
  print('plotting secreted graph')
  #################################################
  svg(file.path(getwd(), "secreted_betweennes_graph.svg"), width=svg_width, height=svg_height)
  par(mar=c(0,0,0,0))
  
  
  
  # print(c(cluster_difexp$important_symbol,cluster_difexp$updown))
  current_g = g
  current_g$layout = norm_coords(current_g$layout,xmin = -scale*2,xmax = scale*2,ymin = -scale,ymax = scale)
  # names(btw) == names(V(current_g))
  rownames(difexp) = difexp$STRING_id
  difexp = difexp[names(V(current_g)),]
  names(V(current_g))==difexp$STRING_id
  
  
  fine = 500000 # this will adjust the resolving power.
  pal = colorRampPalette(c('white','red'))
  graphCol = pal(fine)[as.numeric(cut(difexp$btw,breaks = fine))]
  
  
  plot(current_g,
       vertex.label=NA,
       vertex.color=graphCol,
       vertex.size=abs(difexp$logFC),
       rescale=FALSE)
  
  
  
  
  # a = legend('bottomleft',legend=c("is secreted in placenta","not secreted in placenta"))
  # x <- (a$text$x + a$rect$left) / 2
  # y <- a$text$y
  # symbols(x,y,circles=rep(legend_dot_size,2),inches=FALSE,add=TRUE,bg=c(categorical_pal(1),"white"))
  
  df_layout = as.data.frame(current_g$layout)
  addShadowText(x = df_layout$V1,
                y = df_layout$V2,
                ifelse(difexp$is_secreted ,difexp$SYMBOL,""),
                col="black",
                bg="white",
                cex = 1)

  
  dev.off()
  
  
  print('plotting betweenness graph')
  ############################################33
  
  svg(file.path(getwd(), "betweenness_graph.svg"), width=svg_width, height=svg_height)
  par(mar=c(0,0,0,0))
  
  # print(c(cluster_difexp$important_symbol,cluster_difexp$updown))
  current_g = g
  current_g$layout = norm_coords(current_g$layout,xmin = -scale*2,xmax = scale*2,ymin = -scale,ymax = scale)
  
  btw = betweenness(g, v = V(g), directed = FALSE,
                  weights = NULL, nobigint = TRUE)
  
  
  fine = 500000 # this will adjust the resolving power.
  pal = colorRampPalette(c('blue','green','red'))
  graphCol = pal(fine)[as.numeric(cut(difexp$btw ,breaks = fine))]
  
  # difexp = add_cluster_max_btw_column(difexp, g)
  
  
  
  plot(current_g,
       vertex.label=NA,
       vertex.color=graphCol,
       vertex.size=abs(difexp$logFC),
       rescale=FALSE)
  
  
  # a = legend('bottomleft',legend=c("Up","Down"))
  # x <- (a$text$x + a$rect$left) / 2
  # y <- a$text$y
  # symbols(x,y,circles=rep(legend_dot_size,2),inches=FALSE,add=TRUE,bg=c("green","red"))
  # quantile(difexp[!is.na(difexp$Tissue.RNA...placenta..NX.),]$Tissue.RNA...placenta..NX.,.9)
  
  
  df_layout = as.data.frame(current_g$layout)
  addShadowText(x = df_layout$V1,
                y = df_layout$V2,
                ifelse(difexp$btw>quantile(difexp$btw,.99)[[1]], difexp$SYMBOL ,""),
                col="black",
                bg="white",
                cex = .5)
  # addShadowText(x = df_layout$V1,
  #               y = df_layout$V2,
  #               ifelse(!is.na(difexp$Tissue.RNA...placenta..NX.) & difexp$Tissue.RNA...placenta..NX.>32.62, difexp$SYMBOL ,""),
  #               col="black",
  #               bg="white",
  #               cex = 1)
  dev.off()
  
  
  print('plotting difexp graph')
  ######################################################################################
  svg(file.path(getwd(), "difexp_graph.svg"), width=svg_width, height=svg_height)
  par(mar=c(0,0,0,0))
  
  # print(c(cluster_difexp$important_symbol,cluster_difexp$updown))
  
  g_difexp = g
  g_difexp$layout = norm_coords(g_difexp$layout,xmin = -scale*2,xmax = scale*2,ymin = -scale,ymax = scale)
  V(g_difexp)$color <- colrs[V(g_difexp)$community]
  E(g_difexp)$color <- apply(as.data.frame(get.edgelist(g_difexp)), 1, 
                              function(x) ifelse(get_cluster(g_difexp,x[1]) == get_cluster(g_difexp,x[2]), 
                                                 colrs[get_cluster(g_difexp,x[1])], '#AAAAAAAA'))
  
  plot(g_difexp,
       vertex.label=NA,
       vertex.color=difexp$updown,
       vertex.size=abs(difexp$logFC),
       rescale=FALSE)
  
  a = legend('bottomleft',legend=c("Up","Down"))
  x <- (a$text$x + a$rect$left) / 2
  y <- a$text$y
  symbols(x,y,circles=rep(legend_dot_size,2),inches=FALSE,add=TRUE,bg=c("green","red"))
  df_layout = as.data.frame(g_difexp$layout)
  addShadowText(x = df_layout$V1,
                y = df_layout$V2,
                ifelse(difexp$is_local_max, difexp$SYMBOL,""),
                col="black",
                bg="white",
                cex = 1.5)
  dev.off()
  
  
  print('plotting cluster graph')
  ################################################3
  svg(file.path(getwd(), "cluster_graph.svg"), width=svg_width, height=svg_height)
  
  par(mar=c(0,0,0,0))
  # print(c(cluster_difexp$important_symbol,cluster_difexp$updown))
  
  
  g_cluster = g
  g_cluster$layout = norm_coords(g_cluster$layout,xmin = -scale*2,xmax = scale*2,ymin = -scale,ymax = scale)
  V(g_cluster)$color <- colrs[V(g_cluster)$community]
  
  
  E(g_cluster)$color <- apply(as.data.frame(get.edgelist(g_cluster)), 1, 
                                  function(x) ifelse(get_cluster(g_cluster,x[1]) == get_cluster(g_cluster,x[2]), 
                                                     colrs[get_cluster(g_cluster,x[1])], '#AAAAAAAA'))

  
  table(difexp$cluster)
  
  class(difexp$cluster)
  
  
  plot(g_cluster, vertex.label.dist=0,
       vertex.label.color= "black",
       vertex.label.font=2,
       # vertex.label.cex = 3,
       vertex.label=NA,
       vertex.color=colrs[V(g_cluster)$community],
       vertex.size=abs(difexp$logFC)*2,
       edge.color=E(g_cluster)$color,
       rescale = FALSE)
 
  a = legend('topleft', legend=1:length(colrs), cex=1)
  x <- (a$text$x + a$rect$left) / 2
  y <- a$text$y
  # colrs_2_3 = colrs[c(1,2,3,5,6,7,8)]
  # colrs_2_3
  symbols(x,y,
          circles=rep(legend_dot_size,length(colrs)),
          inches=FALSE,add=TRUE,bg=colrs)
  
  
  # print(difexp[difexp$cluster_max_btw!="",c("SYMBOL","cluster")])
  
  
  difexptest = difexp
  rownames(difexptest) = difexp$STRING_id
  difexptest =  difexptest[names(V(g_cluster)),]
  !(FALSE %in% (rownames(difexptest) == names(V(g_cluster))))
  
  
  # map_clusters_2_3 = c("1","2","3","","4","5","6","7","","","","","","","")
  
  df_layout = as.data.frame(g_cluster$layout)
  
  #map g_cluster v to difexp symbol
  symbols = sapply(names(V(g_cluster)),function(x) as.character(difexp[difexp$STRING_id == x,]$SYMBOL))
  symbols
  
  addShadowText(x = df_layout$V1,
                y = df_layout$V2,
                symbols,
                # ifelse(difexptest$is_local_max, difexptest$cluster,""),
                col="black",
                bg="white",
                cex = .55)
  

  


  dev.off()
  ############################################33
  
  svg(file.path(getwd(), "enrichment_coverage_graph.svg"), width=svg_width, height=svg_height)
  par(mar=c(0,0,0,2.1))
  g_coverage = g
  g_coverage$layout = norm_coords(g_coverage$layout,xmin = -scale*2,xmax = scale*2,ymin = -scale,ymax = scale)
  V(g_coverage)$color <- colrs[V(g_coverage)$community]
  E(g_coverage)$color <- apply(as.data.frame(get.edgelist(g_coverage)), 1, 
                              function(x) ifelse(get_cluster(g_coverage,x[1]) == get_cluster(g_coverage,x[2]), 
                                                 colrs[get_cluster(g_coverage,x[1])], '#AAAAAAAA'))
  
  difexp$covered = ifelse(difexp$boyan==0,"black","green")
  
  plot(g_coverage,
       vertex.label=NA,
       vertex.color=difexp$covered,
       vertex.size=abs(difexp$logFC),
       rescale=FALSE)
  
  a = legend('bottomleft',legend=c("covered","not covered"))
  x <- (a$text$x + a$rect$left) / 2
  y <- a$text$y
  symbols(x,y,circles=rep(legend_dot_size,2),inches=FALSE,add=TRUE,bg=c("green","black"))
  
  df_layout = as.data.frame(g_cluster$layout)
  addShadowText(x = df_layout$V1,
                y = df_layout$V2,
                # symbols,
                ifelse(difexp$boyan!=0,difexp$SYMBOL, ""),
                col="black",
                bg="white",
                cex = .6)
  
  dev.off()
}

