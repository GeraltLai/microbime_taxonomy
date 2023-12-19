##########################
####### function #########
##########################
# log-normal data
log_normal_data = function(group, sample, case_genus, control_genus, mu1, delta1){
  sigma1 = 1
  case_mu2 = mu1 + delta1
  control_mu2 = mu1
  sigma2 = sigma1
  genus= case_genus + control_genus
  genus_OTU = 3
  total_sample = sample * group
  total_OTU = genus * genus_OTU
  
  # log normal simulation data
  control_data = rlnorm(n = total_OTU * total_sample, meanlog = mu1, sdlog = sigma1)
  OTU = matrix(data = control_data, nrow = total_sample, ncol = total_OTU )
  case_data = rlnorm(case_genus * genus_OTU * sample, meanlog = case_mu2, sdlog = sigma2)
  case_matrix = matrix(data = case_data, nrow = sample, ncol = case_genus * genus_OTU )
  OTU[(sample + 1) : (2 *sample), 1: (case_genus * genus_OTU)] = case_matrix 
  row_name = c()
  for(i in 1:nrow(OTU)){
    row_name = c(row_name, paste0('sample',i))
  }
  OTU_name = c()
  for(i in 1:ncol(OTU)){
    OTU_name = c(OTU_name, paste0('OTU',i))
  }
  colnames(OTU) = OTU_name
  row.names(OTU) = row_name
  OTU = as.data.frame(round(OTU,0))
  index = OTU == 0
  OTU[index] = 1
  return(OTU)
}

# meta data simulation
meta_data = function(group, sample){
  factor = rep(letters[1:group], each = sample)
  meta = matrix(data = factor, nrow = group * sample)
  row_name = sapply(1:nrow(meta), function(x) paste0('sample',x))
  colnames(meta) = "factor"
  row.names(meta) = row_name
  meta = as.data.frame(meta)
  return(meta)
}

# tax data simulation
tax_data = function(case_genus, control_genus){
  genus = case_genus + control_genus
  total_OTU = 3 * genus
  tax_family_name = rep("family", total_OTU)
  tax_genus_name = sapply(1:genus, function(x) rep(paste0('genus', x), 3))
  tax_genus_name = as.vector(tax_genus_name)
  OTU_name = sapply(1:total_OTU, function(x) paste0('OTU', x))
  tax = data.frame(tax_family_name,tax_genus_name )
  row.names(tax) = OTU_name
  colnames(tax) = c("family", "genus")
  return(tax)
}

# weighted  bray-curits function
weighted_bray = function(data){
  matrix = diag(x = 1, nrow = nrow(data))
  for(ii in 1:(nrow(data)-1))
    for(jj in (ii+1):nrow(data)){
      A = data[ii,] ; B = data[jj,]
      weight = sapply(1:length(A), function(x) 1/((A[x] + B[x])/sum(A,B)))
      top = sapply(1:length(A), function(x) weight[x] * abs(A[x]-B[x]))
      down = sapply(1:length(A), function(x) weight[x] * (A[x] + B[x]))
      top = do.call(cbind, top) ; down = do.call(cbind, down)
      matrix[ii,jj] = sum(top)/sum(down)
    }
  label = t(matrix)[lower.tri(matrix)]
  # After transposing, remove the triangle and fill the lower triangle
  matrix[lower.tri(matrix)]=t(matrix)[lower.tri(matrix)]
  return(matrix)
}
# compare result
compare_result = function(all_result, pairwise_result){
  genus = all_result[-1,]
  n = nrow(genus)
  table = matrix(data = 0, nrow = n, ncol = n)
  for(i in 1:n){
    output = rep(0,n)
    step1 = genus[i,4]
    output[i] = step1
    name = row.names(genus)[i]
    id = grep(name, row.names(pairwise_result))
    step2 = pairwise_result[id, 4]
    output[-i] = step2
    table[i,] = output
  }
  row.names(table) = row.names(genus)
  colnames(table) = row.names(genus)
  return(table)
}
# count table
count_table = function(data, step1_result){
  step1 = step1_result[-1,]
  sig_id = step1[,4]<0.05 
  sig_gene = colnames(data)[sig_id]
  nonsig_id = step1[,4]>=0.05 
  nonsig_gene = colnames(data)[nonsig_id]
  individual_count = c()
  with_sig = c()
  with_nonsig = c()
  total_count = c()
  for(i in 1:nrow(data)){
    gene = row.names(data)[i]
    if(sum(str_detect(gene, sig_gene)) == 1 ){
      new_sig_gene = sig_gene[-grep(gene, sig_gene)]
      new_nonsig_gene = nonsig_gene
    }else{
      new_sig_gene = sig_gene
      new_nonsig_gene = nonsig_gene[-grep(gene, nonsig_gene)]
    }
    individual = ifelse(data[i,i] < 0.05, 1, 0)
    individual_count = c(individual_count, individual)
    count_sig = ifelse(data[i, new_sig_gene]<0.05, 1, 0)
    with_sig = c(with_sig, sum(count_sig))
    count_nonsig = ifelse(data[i, new_nonsig_gene]<0.05, 1, 0)
    with_nonsig = c(with_nonsig, sum(count_nonsig))
    total = sum(data[i,]< 0.05)
    total_count = c(total_count, total)
  }
  output = data.frame(individual_count, with_sig, with_nonsig, total_count)
  row.names(output) = row.names(data)
  return(output)
}
# merge_diff_mu_result
merge_diff_mu_result = function(data, mu){
  merge_table = list()
  all_list_name = names(data)
  for (i in seq_along(mu1_test)) {
    mu_name = paste0("mu=", mu1_test[i], "\\|")
    mu_index = grep(mu_name, all_list_name)
    mu_list = data[mu_index]
    
    # total_count table
    total_count = do.call(cbind, lapply(mu_list, function(x) x$total_count))
    col_name = gsub(mu_name, "", colnames(total_count))
    row.names(total_count) = unique(c(tax$family, tax$genus))
    
    # individual_count table
    individual_count = do.call(cbind, lapply(mu_list, function(x) x$individual_count))
    colnames(individual_count) = col_name
    row.names(individual_count) = unique(c(tax$family, tax$genus))
    
    # with_sig table
    with_sig_count = do.call(cbind, lapply(mu_list, function(x) x$with_sig))
    colnames(with_sig_count) = col_name
    row.names(with_sig_count) = unique(c(tax$family, tax$genus))
    
    # with_nonsig table
    with_nonsig_count = do.call(cbind, lapply(mu_list, function(x) x$with_nonsig))
    colnames(with_nonsig_count) = col_name
    row.names(with_nonsig_count) = unique(c(tax$family, tax$genus))
    
    # p-value table
    p_value_table = do.call(cbind, lapply(mu_list, function(x) x$`Pr(>F)`))
    colnames(p_value_table) = col_name
    row.names(p_value_table) = unique(c(tax$family, tax$genus))
    
    # size table
    size_table = do.call(cbind, lapply(mu_list, function(x) x$size))
    colnames(size_table) = col_name
    row.names(size_table) = unique(c(tax$family, tax$genus))
    
    ### merge to list
    merge_table[[mu_name]] = list(total = total_count, individual = individual_count,
                                  with_sig = with_sig_count,  with_nonsig = with_nonsig_count,
                                  p_value = p_value_table, size = size_table)
  }
  return(merge_table)
}
# distinguish different genus result
cluster_genus_result = function(loop_result){
  output = list()
  loop = length(loop_result)
  col_length = length(loop_result[[1]]) + 1
  for(i in unique(tax$genus)){
    pick_gene = i
    result_matrix = matrix(0, ncol = col_length, nrow = length(delta1_test) * loop)
    col_size = sapply(loop_result, function(x) x$size[pick_gene,])
    col_p = sapply(loop_result, function(x) x$p_value[pick_gene,])
    col_individual = sapply(loop_result, function(x) x$individual[pick_gene,])
    col_with_sig = sapply(loop_result, function(x) x$with_sig[pick_gene,])
    col_with_nonsig = sapply(loop_result, function(x) x$with_nonsig[pick_gene,])
    col_total = sapply(loop_result, function(x) x$total[pick_gene,])
    col_delta = rep(delta1_test,loop)
    result_matrix[,1] = col_size
    result_matrix[,2] = col_p 
    result_matrix[,3] = col_individual
    result_matrix[,4] = col_with_sig
    result_matrix[,5] = col_with_nonsig
    result_matrix[,6] = col_total
    result_matrix[,7] = col_delta
    result_matrix = as.data.frame(result_matrix)
    colnames(result_matrix) = c("size", "p-value", "individual", "with_sig", 
                                "with_nonsig", "total","delta1")
    result_matrix = arrange(result_matrix, delta1)
    output[[i]] = result_matrix
  }
  return(output)
}
# significant table
significant_table = function(data, tax_table, delta, case_genus, control_genus){
  output = list()
  total_genus = case_genus + control_genus
  for(ii in 1:length(data)){
    input = data[[ii]]
    id = which(input$individual == 1)
    output_table = matrix(data = 0, ncol = 3, nrow = length(delta))
    output_table[,1] = delta
    count_with_sig = c()
    count_with_nonsig = c()
    count_significant = c()
    new_data = input[id,]
    for(i in seq_along(delta)){
      pick = new_data[new_data$delta1 == delta[i],]
      output_table[i,2] = dim(pick)[1]
      output_table[i,3] = round(mean(pick$size),0)
      count_with_sig = rbind(count_with_sig, table(factor(pick$with_sig, levels = 0:(total_genus-1))))
      count_with_nonsig = rbind(count_with_nonsig,table(factor(pick$with_nonsig, levels = 0:(total_genus-1))))
      count_significant = rbind(count_significant, table(factor(pick$total, levels = 1:total_genus)))			
    }
    output_table = cbind(output_table, count_with_sig,count_with_nonsig, count_significant)
    output_table = as.data.frame(output_table)
    output[[ii]] = output_table
  }
  names(output) = unique(tax_table$genus)
  return(output)
}
# non-significant table
nonsignificant_table = function(data, tax_table, delta, case_genus, control_genus){
  output = list()
  total_genus = case_genus + control_genus
  for(ii in 1:length(data)){
    input = data[[ii]]
    id = which(input$individual == 0)
    output_table = matrix(data = 0, ncol = 3, nrow = length(delta))
    output_table[,1] = delta
    count_with_sig = c()
    count_with_nonsig = c()
    count_significant = c()
    new_data = input[id,]
    for(i in seq_along(delta)){
      pick = new_data[new_data$delta1 == delta[i],]
      output_table[i,2] = dim(pick)[1]
      output_table[i,3] = round(mean(pick$size),0)
      count_with_sig = rbind(count_with_sig, table(factor(pick$with_sig, levels = 0:(total_genus-1))))
      count_with_nonsig = rbind(count_with_nonsig,table(factor(pick$with_nonsig, levels = 0:(total_genus-1))))
      count_significant = rbind(count_significant, table(factor(pick$total, levels = 0:(total_genus-1))))			
    }
    output_table = cbind(output_table, count_with_sig,count_with_nonsig, count_significant)
    output_table = as.data.frame(output_table)
    output[[ii]] = output_table
  }
  names(output) = unique(tax_table$genus)
  return(output)
}
# percentage significant table
percentage_significant_table = function(data, tax_table, case_genus, control_genus){
  output = list()
  total_genus = case_genus + control_genus
  for(ii in 1:length(data)){
    new_data = data[[ii]]
    new_output = c()
    for(i in 1:nrow(new_data)){
      total = new_data[i,2]
      data1 = new_data[i,4: (4+total_genus-1)]
      data1 = round(data1/total, 4)
      data2 = new_data[i,(4+total_genus):(4+ 2*total_genus-1)]
      data2 = round(data2/total, 4)
      data3 = new_data[i,(4+ 2*total_genus):(4+ 3*total_genus -1)]
      data3 = round(data3/total, 4)
      new_output = rbind(new_output, cbind(data1,data2,data3))
    }
    output[[ii]] = new_output
  }
  names(output) = unique(tax_table$genus)
  return(output)
}
# percentage nonsignificant table
percentage_nonsignificant_table = function(data, tax_table, case_genus, control_genus){
  output = list()
  total_genus = case_genus + control_genus
  for(ii in 1:length(data)){
    new_data = data[[ii]]
    new_output = c()
    for(i in 1:nrow(new_data)){
      total = new_data[i,2]
      data1 = new_data[i,4: (4+total_genus-1)]
      data1 = round(data1/total, 4)
      data2 = new_data[i,(4+total_genus):(4+ 2*total_genus-1)]
      data2 = round(data2/total, 4)
      data3 = new_data[i,(4+ 2*total_genus):(4+ 3*total_genus -1)]
      data3 = round(data3/total, 4)
      new_output = rbind(new_output, cbind(data1,data2,data3))
    }
    output[[ii]] = new_output
  }
  names(output) = unique(tax_table$genus)
  return(output)
}



