######
rm(list = ls())

#### start
source("C:/Users/lab205/Desktop/lab205/Microbiota_Analysis_BRC/simulation_data_analysis/function.R")

# simulate data
mu1_test = c(1)
delta1_test = seq(from = 0.6, to = 2.0, by = 0.2)
group = 2 
sample = 4 
case_genus = 3
control_genus = 3

# loop parameter
mu1 = mu1_test[1]
loop = 10
loop_table = list()

# negative record
negative_count = 0
negative_delta = c()
negative_OTU = list()
negative_step2_output = list()

### basic
weighted_bray = function(data){
  matrix = diag(x = 1, nrow = nrow(data))
  for(ii in 1:(nrow(data)-1))
    for(jj in (ii+1):nrow(data)){
      A = data[ii,] ; B = data[jj,]
      value = sapply(1:length(A), function(x) abs(A[x]-B[x])/sum(A,B))
      value = do.call(cbind, value)
      matrix[ii,jj] = sum(value)
    }
  label = t(matrix)[lower.tri(matrix)]
  # After transposing, remove the triangle and fill the lower triangle
  matrix[lower.tri(matrix)]=t(matrix)[lower.tri(matrix)]
  return(matrix)
}

##########################
###### for loop ##########
##########################
for(l in 1:loop){
  final_result = list()
  for(aa in 1:length(delta1_test)){
    # set parameter
    delta1 = delta1_test[aa]
    print(paste("The loop is", l,". The delta1 is ", delta1))
    
    # simulate data
    OTU = log_normal_data(group = group, 
                          sample = sample, 
                          case_genus = case_genus, 
                          control_genus = control_genus,
                          mu1 = mu1,
                          delta1 = delta1 )
    meta = meta_data(group = group, sample = sample)
    tax = tax_data(case_genus = case_genus, control_genus = control_genus)
    
    ###############################################
    ###### step1 : individual genus output ########
    ###############################################
    step1_output = c()
    dist = bray0(OTU)
    result = adonis2(dist ~ factor, data = meta)[1,3:5]
    size = round(sum(OTU),0)
    result = cbind(size, result)
    step1_output = rbind(step1_output, result)
    genus = unique(tax$genus)
    for(i in genus){
      genus_group = subset(tax, tax[["genus"]] == i)
      genus_id = row.names(genus_group)
      new_OTU = select(OTU, all_of(genus_id))
      size = round(sum(new_OTU), 0)
      dist = bray0(new_OTU)
      result = adonis2(dist ~ factor, data = meta)[1,3:5]
      result = cbind(size, result)
      step1_output = rbind(step1_output, result)
    }
    row.names(step1_output) = c("total", genus)
    colnames(step1_output) = c("size", "R2", "F", "Pr(>F)")
    step1_output = as.data.frame(step1_output)
    
    #############################################
    ###### step2 : pairwise genus output ########
    #############################################
    step2_output = c()
    name = c()
    length = length(unique(tax$genus))
    for(i in 1:(length-1))
      for(j in (i+1):length){
        genus_id1 = unique(tax$genus)[i]
        genus_id2 = unique(tax$genus)[j]
        name = c(name, paste0(genus_id1,"-",genus_id2))
        genus_group1 = subset(tax, tax[["genus"]] == genus_id1)
        genus_group2 = subset(tax, tax[["genus"]] == genus_id2)
        genus_OTU = c(row.names(genus_group1),row.names(genus_group2))
        clean = select(OTU, all_of(genus_OTU))
        size = sum(clean)
        dist = weighted_bray(data = clean)
        ans = adonis2(dist ~ factor, data = meta)[1,3:5]
        output = round(c(size, ans[1,1], ans[1,2], ans[1,3]), 4)
        step2_output = rbind(step2_output, output)
      }
    colnames(step2_output) = c("size", "R2", "F", "Pr(>F)")
    row.names(step2_output) = name
    step2_output = as.data.frame(step2_output)
    if(sum(step2_output[,"F"]<0) > 0 ){
      negative_count = negative_count + 1
      negative_delta = c(negative_delta, delta1)
      negative_OTU[[negative_count]] = OTU
      negative_step2_output[[negative_count]] = step2_output
    }
    
    # compare result
    all_result = compare_result(all_result = step1_output, pairwise_result = step2_output) 
    
    # summary result
    individual_result = step1_output[,c(1,4)]
    my_result = count_table(data = all_result, step1_result = step1_output)
    if(individual_result[1,2] >= 0.05){
      my_result = rbind(rep(0,4), my_result)
    }else{
      my_result = rbind(rep(1,4), my_result)
    }
    summary_result = cbind(individual_result, my_result)
    row.names(summary_result)[1] = "family"
    list_name = paste0("mu=",mu1,"|","delta=",delta1)
    final_result[[list_name]] = summary_result
  }
  # merge to table
  merge_table = merge_diff_mu_result(data = final_result, mu = mu1_test)
  # loop_table
  loop_table[l] = merge_table
}
loop_table

# distinguish different genus result
result_list = cluster_genus_result(loop_result = loop_table)
result_list[[1]]

# significant genus table
significant_genus_result = significant_table(data = result_list, tax_table= tax, delta = delta1_test, case_genus = case_genus, control_genus = control_genus)
percentage_significant_genus_result = percentage_significant_table(data = significant_genus_result, tax_table= tax, case_genus = case_genus, control_genus = control_genus)

# nonsignificant genus table
nonsignificant_genus_result = nonsignificant_table(data = result_list, tax_table= tax, delta = delta1_test, case_genus = case_genus, control_genus = control_genus)
percentage_nonsignificant_genus_result = percentage_nonsignificant_table(data = nonsignificant_genus_result, tax_table= tax, case_genus = case_genus, control_genus = control_genus)

