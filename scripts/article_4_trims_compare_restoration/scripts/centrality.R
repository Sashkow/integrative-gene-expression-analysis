
add_btw_column = function(difexp, g){
  difexp = difexp[,!(names(difexp) %in% c("btw"))]
  btw = betweenness(g, v = V(g), directed = FALSE,
                    weights = NULL, nobigint = TRUE)
  
  btw_df = as.data.frame(btw)
  btw_df$string_id = rownames(btw_df)
  nrow(btw_df) == nrow(difexp)
  difexp_m = merge(difexp,btw_df, by.x = "STRING_id", by.y = "string_id")
  difexp_m[order(difexp_m$btw),]$STRING_id == btw_df[order(btw_df$btw),]$string_id
  
  difexp = difexp_m
  return(difexp)
}


# difexp = difexp_mapped
# g = sub_g
# nrow(difexp)
# length(names(V(g)))

add_cluster_max_btw_column = function(difexp, g, string_db = NULL){

  #adds btw column as well

  difexp = difexp[,!(names(difexp) %in% c("cluster_max_btw"))]
  length(V(g)) == nrow(difexp)
  difexp = add_btw_column(difexp, g)
  difexp$cluster_max_btw = ""
  difexp[which(difexp$btw == max(difexp$btw)),]$cluster_max_btw = paste(
    "global",
    difexp[which(difexp$btw == max(difexp$btw)),]$SYMBOL,
    difexp[which(difexp$btw == max(difexp$btw)),]$btw)
  i = 1

  for (i in 1:length(unique(difexp$cluster))){
    cluster_difexp = difexp[which(difexp$cluster==i),]
    sub_g = simplify(string_db$get_subnetwork(cluster_difexp$STRING_id))
    cluster_difexp = add_btw_column(cluster_difexp, sub_g)
    
    max_gene = head(cluster_difexp[which(cluster_difexp$btw == max(cluster_difexp$btw)),],1)
    max_gene$SYMBOL
    
    difexp[which(difexp$STRING_id == max_gene$STRING_id),]$cluster_max_btw = paste(
      max_gene$cluster,
      max_gene$SYMBOL,
      round(max_gene$btw), sep = ", ")
    difexp[which(difexp$STRING_id == max_gene$STRING_id),]$cluster_max_btw
  }
  return(difexp)
}

# difexp = difexp_mapped

add_is_local_abs_max_column = function(difexp){
  difexp$is_local_max = FALSE
  i = 1
  for (i in 1:max(difexp$cluster)){
    print(i)
    cluster_max = max(abs(difexp[difexp$cluster == i,]$logFC))
    cluster_max
    max_gene = difexp[abs(difexp$logFC) == cluster_max,]$STRING_id[1]
    max_gene
    difexp[difexp$STRING_id==max_gene,]$is_local_max = TRUE
  }
  return(difexp)
}
