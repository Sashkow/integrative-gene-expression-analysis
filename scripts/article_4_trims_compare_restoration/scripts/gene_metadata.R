add_is_secreted_placenta_column = function(difexp, target_path){
  list.files(target_path)
  secreted_genes = read.csv(file.path(base_dir, "data", "secreted_genes_placenta.csv"))
  nrow(secreted_genes)
  
  difexp$is_secreted = ifelse(difexp$SYMBOL %in% secreted_genes$Gene,TRUE, FALSE)
  return(difexp)
}


# difexp = difexp_mapped

add_hpa_columns = function(difexp, target_path){
  data("hpaNormalTissue")
  data(rnaGeneTissue)
  data(hpaSecretome)
  difexp$SYMBOL
  difexp = merge(difexp,hpaSecretome,by.x="SYMBOL", by.y="Gene.name", all.x=TRUE)
  difexp
  nrow(difexp)
  first_column = "Tissue.RNA...adipose.tissue..NX."
  last_column = "Tissue.RNA...total.PBMC..NX."
  first_column_i = which(colnames(difexp) == first_column)
  last_column_i = which(colnames(difexp) == last_column)
  first_column_i
  last_column_i
  df = difexp[,first_column_i:last_column_i]
  
  my_which_max = function(x,difexp){
    if (is.na(max(x))){
      return("Uknown")
    }
    tmp = colnames(difexp)[which.max(x)]
    tmp = strsplit(tmp,"\\.\\.\\.")[[1]][2]
    tmp = strsplit(tmp,"\\.\\.")[[1]][1]
    return(tmp)
  }
  difexp$max_tissue = apply(df,1,my_which_max, difexp=df)
  return(difexp)
}





# path = "/home/sashkoah/a/r/igea-r/article_4_trims_compare/placenta/1_2_no_55439_100/protein_class_Transcription.tsv"
# transcription = read.table(path, sep = '\t', header = TRUE)
# transcription$Gene
# difexp$is_transcription_factor = ifelse(difexp$SYMBOL %in% transcription$Gene, TRUE, FALSE)
# 
# difexp
# 
# difexp[which(difexp$is_secreted & difexp$is_transcription_factor),]
# 
# property = "is_secreted"
# numeric_property = "boyan"
# more = TRUE
# is_property_has_outstanding_numeric_property(difexp,property, numeric_property, more)




is_property_has_outstanding_numeric_property = function(difexp, property, numeric_property, more=TRUE){
  
  base_value = mean(abs(difexp[difexp[,property] == TRUE,][,numeric_property]))
  total_average = mean(abs(difexp[,numeric_property]))
  
  amount_of_genes_with_property = nrow(difexp[difexp[,property]==TRUE,])
  amount_of_genes_with_property
  
  cases_outstanding = 0
  trials = 1000
  for (i in 1:trials){
    random_genes = difexp[sample(nrow(difexp), amount_of_secreted_genes), ]
    super_value = mean(abs(random_genes[,numeric_property]))
    if (more){
      if (super_value > base_value){
        cases_outstanding = cases_outstanding + 1
      }
    }else{
      if (super_value < base_value){
        cases_outstanding = cases_outstanding + 1
      }
    }
  }
  
  alpha = .05
  p = cases_outstanding/trials
  print(paste(ifelse(p<alpha,"Passed:","Failed:"),
              "On average",property,"genes have",
              ifelse(more,"HIGHER","LOWER"),
              "abs",numeric_property,"(",round(base_value,2),")",
              "than abs average among all difexp genes (",round(total_average,2),").",
              ifelse(p==0,paste("p <",1/trials),paste("p =", p))))
}

