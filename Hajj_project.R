library(foreign)
#library(Hmisc)
library(ivpack)
library(onehot)
library(survey)
library(dplyr)
library(FindIt)
library(dagitty)
g = dagitty("dag{
            Y<-X->Z
            X->D->Y
            Z->D")
plot(graphLayout(g))
hajj = read.dta('/Users/shashalin/Documents/CausalInference/Hajj/Hajj\ Distribution/hajj_public.dta')
vars = colnames(hajj)
#These are to construct the cell variable, controlling for differences in group
hajj$pod = as.factor(hajj$pod)
hajj$ptygrp = as.factor(hajj$ptygrp)
hajj$cat = as.factor(hajj$cat)
hajj$district = as.factor(hajj$district)
# pod_encoder = onehot(hajj['pod'])
# pod = predict(pod_encoder, hajj['pod'])
# ptygrp_encoder = onehot(hajj['ptygrp'])
# cat_encoder = onehot(hajj['cat'])

# ptygrp = predict(ptygrp_encoder, hajj['ptygrp'])
# cat = predict(cat_encoder, hajj['cat'])

other_vars = vars %in% c('pod', 'ptygrp', 'cat', 'hhid', 'persid', 'success', 'hajj2006', 'clusterid', 'district', 'subsample', 'female', 'age', 'literate', 'urban', 'smallpty', 'newcluster')
deps = vars[!other_vars]
library(standardize)
for(dep in deps){hajj[dep] = as.numeric(scale(hajj[dep]))}

###independent tests, not good coz too many tests!
# independent_models = list()
# coefs = list()
# i=1
# for(j in seq(1, length(deps))){
#   ivmodel = ivreg(deps_values[,j] ~ hajj2006, instruments=~success, data=hajj)
#   nonmissing = rep(TRUE, dim(hajj)[1])
#   nonmissing[unname(ivmodel$na.action)] = FALSE
#   coefs[[deps[j]]] = ivmodel$coefficients['hajj2006']
#   independent_models[[deps[j]]] = cluster.robust.se(ivmodel, clusterid=hajj$clusterid[nonmissing])}
########

##SUR
##results are somewhat different from those reported in original paper, not exactly sure why
# library(tidyr)
# hajj_orthodox=gather(hajj, key='question', value='value', x_s14aq1,x_s14aq3,x_s14aq4,x_s14aq5,x_s14aq6,x_s14bq2,x_s14aq8,x_s14aq9,x_s14aq12,x_s14aq13, factor_key = )
# 
# hajj_out = interaction(hajj_orthodox$hajj2006, hajj_orthodox$question)
# success_out = interaction(hajj_orthodox[,'success'], hajj_orthodox$question)
# covariates = interaction(hajj_orthodox$pod, hajj_orthodox$cat, hajj_orthodox$ptygrp, hajj_orthodox$question)
# orthodox_sur = ivreg(hajj_orthodox$value~hajj_out+covariates, ~success_out+covariates)
# orthodox_coefs = orthodox_sur$coefficients[2:20]
# mean(orthodox_coefs)

####################################
#AES as done in original paper
index_AES = function(cols, name){
  dependent_vars = hajj[, cols]
  if(length(cols) > 1){dependent_var = rowMeans(dependent_vars, na.rm=TRUE)}
  else{dependent_var = dependent_vars}
  model=ivreg(dependent_var ~ hajj$hajj2006 +hajj$ptygrp*hajj$cat*hajj$pod,
                           ~ hajj$success+ hajj$ptygrp*hajj$cat*hajj$pod, x=TRUE)
  model=ivreg(dependent_var ~ hajj$hajj2006 +hajj$ptygrp*hajj$cat*hajj$pod,
              ~ hajj$success+ hajj$ptygrp*hajj$cat*hajj$pod, x=TRUE)
  nonmissing = rep(TRUE, dim(hajj)[1])
  nonmissing[unname(model$na.action)] = FALSE
  model_result = cluster.robust.se(model, clusterid=hajj$clusterid[nonmissing])
  return(model_result)
}

AES_dic = list()
make_AES_dic = function(name, cols, all_measures=AES_dic, num_terms=37){
  all_measures[[name]][['cols']] = cols
  result = index_AES(cols, name)
  #all_measures[[name]][['result']] = result
  all_measures[[name]][['beta']] = result[2]
  all_measures[[name]][['std_err']] = result[2+num_terms]
  all_measures[[name]][['t']] = result[2+2*num_terms]
  all_measures[[name]][['p']] = result[2+3*num_terms]
  return(all_measures)
}
AES_dic = make_AES_dic('religious', c('x_s14aq10'), AES_dic)

AES_dic = make_AES_dic('belief_local', c('x_s14aq15d',
                               'x_s14aq16b',
                               'x_s14aq16a',
                               'x_s14aq16d',
                               'x_s14aq16e',
                               'x_s14aq16f', 
                               'x_s14dq6',
                               'x_s14cq6',
                               'x_s10dq7',
                               'x_s14dq7'), AES_dic)
AES_dic = make_AES_dic('belief_orthodox', c('x_s14aq1', 'x_s14aq3','x_s14aq4','x_s14aq5', 
             'x_s14aq6',
             'x_s14bq2', 
             'x_s14aq8',
             'x_s14aq9',
             'x_s14aq12',
             'x_s14aq13'), AES_dic)

AES_dic = make_AES_dic('participate_local', c('x_2_s14aq16b',
             'x_2_s14aq16a',
             'x_2_s14aq16d',
             'x_2_s14aq16f'), AES_dic)

AES_dic = make_AES_dic('view_other_countries', c('x_s15cq5', 
                                       'x_s15cq6',
                                       'x_s15cq7',
                                       'x_s15cq8',
                                       'x_s15cq10',
                                       'x_s15cq9' ), AES_dic)

AES_dic = make_AES_dic('view_other_group', c('x_s15b_b6',  
                                   'x_s15b_c6',
                                   'x_s15b_a6'), AES_dic)

AES_dic = make_AES_dic('harmony', c('x_s15cq2',
                          'x_s15cq3',
                          'x_s15cq1',
                          'x_s14aq7'), AES_dic)

AES_dic = make_AES_dic('peaceful', c('x_s10bq4',
                           'x_s10bq5',
                           'x_s10cq5', 
                           'x_s10cq4', 
                           'x_s10cq3',
                           'x_s10hq6',
                           'x_s10hq7',
                           'x_s10hq8'), AES_dic)

AES_dic =  make_AES_dic('political_islam', c('x_s7q12a',
                                  'x_s7q12f',
                                  'x_s7q12d',
                                  'x_s7q12c',
                                  'x_s7q11c'), AES_dic)

AES_dic = make_AES_dic('view_west', c('x_s10aq3',
                            'x_s10aq2',
                            'x_s10aq4',
                            'x_s10aq5'), AES_dic)

AES_dic = make_AES_dic('view_women', c('x_s10dq15a',
                             'x_s10dq15b',
                             'x_s10dq15c',
                             'x_s10dq18'), AES_dic)

AES_dic = make_AES_dic('women_QoL', c('xda_s10dq19b', 
                            'xda_s10dq19c', 
                            'xda_s10dq19d', 
                            'x_2_s10dq8', 
                            'xd_s10dq8'), AES_dic)

AES_dic = make_AES_dic('girl_education', c('x_s10eq2',
                                 'x_s10eq5', 
                                 'x_s10eq4', 
                                 'x_s10eq6',  
                                 'x_s10eq3' ), AES_dic)

AES_dic = make_AES_dic('women_in_work', c('x_s1_3q3',
                                'x_s1_3q4',
                                'x_s1_3q5r4'), AES_dic)

AES_dic = make_AES_dic('gender_authority', c('x_s10dq15d',
                                   'x_s10dq6',
                                   'x_s10hq3', 
                                   'x_s10hq11',
                                   'x_s10dq1',
                                   'x_s10dq14',
                                   'x_2_s14cq8_9'), AES_dic)

AES_dic = make_AES_dic('K6', c('s3q1a', 
                     's3q1b',
                     's3q1c',
                     's3q1d',
                     's3q1e', 
                     's3q1f'), AES_dic)

AES_dic = make_AES_dic('positive_feeling', c('x_s3q2a',
                                   'x_s3q2b', 
                                   'x_s3q2c',
                                   'x_s3q7',
                                   'x_s3q3'), AES_dic)

AES_dic = make_AES_dic('life_satisfaction', c('x_s3q4', 
                                    'x_s3q5',
                                    'x_s3q6'), AES_dic)

AES_dic = make_AES_dic('physical_health', c('x_s2q1', 
                                  'x_s2q2' ),AES_dic)

AES_dic = make_AES_dic('socio_economic_engagement', c('x_s5_5q1',
                                            'x_s5_5q2',
                                            'x_s5_5q3',
                                            'x_s5_5q4',
                                            's5_2q1a',
                                            's5_2q2a',
                                            's5_2q3a',
                                            's5_2q1b',
                                            's5_2q2b',
                                            's5_2q3b',
                                            'x_s5_4q1a',
                                            'x_s5_4q2a',
                                            'x_s5_4q3a', 
                                            'x_s12q2',
                                            'x_s12q1'),AES_dic)

AES_dic = make_AES_dic('engage_politics', c('x_s7q7',
                                  'x_s7q1',
                                  'x_s7q2',
                                  'x_s5_4q5a',
                                  'x_s5_4q4a',
                                  'x_s7q10',
                                  'x_s7q9'), AES_dic)

AES_dic = make_AES_dic('islam_knowledge', c('x_pillars',
                                  'x_s14bq5', 
                                  'x_s14bq7',
                                  'x_s14cq2',
                                  'x_s14cq1', 
                                  'x_s14dq2', 
                                  'x_s14dq1', 
                                  'x_s14dq5', 
                                  'x_s14bq3', 
                                  'x_s14cq3'),AES_dic)

AES_dic = make_AES_dic('diversity_knowledge', c('x_s14dq3', 
                                     'x_2_s14dq6', 
                                     'x_s14dq4'), AES_dic)

AES_dic = make_AES_dic('gender_knowledge', c('x_s14bq4',
                                   'x_s14bq6',
                                   'x_s14cq4',  
                                   'x_2_s14cq6', 
                                   'x_s10dq12',
                                   'xop_s10dq19b',
                                   'xop_s10dq19c',
                                   'xop_s10dq19d'), AES_dic)

AES_dic = make_AES_dic('global_knowledge', c('x_s8q2',
                                   'x_s8q3', 
                                   'x_s8q7',
                                   'x_s8q8',
                                   'x_s8q6',
                                   'x_s8q5'), AES_dic)
pvals = list()
for(x in names(AES_dic)){
  pvals[x] = AES_dic[[x]]['p']
}
holm_ps = p.adjust(pvals, method='holm', n=25)
sum(holm_ps>.05) #14
BH_ps = p.adjust(pvals, method='BH', n=25)
sum(BH_ps >.05)

##############################################
#Post Double Selection
library(stringr)
ave_dependent = function(name, cols, df=hajj){
  dependent_vars = hajj[, cols]
  if(length(cols) > 1){dependent_var = rowMeans(dependent_vars, na.rm=TRUE)}
  else{dependent_var = dependent_vars}
  new_name=str_c('ave', name, sep='_')
  print(new_name)
  df[,new_name] = dependent_var
  return(df)
  }

hajj = ave_dependent('religious', c('x_s14aq10'))
hajj = ave_dependent('belief_local', c('x_s14aq15d',
                                         'x_s14aq16b',
                                         'x_s14aq16a',
                                         'x_s14aq16d',
                                         'x_s14aq16e',
                                         'x_s14aq16f', 
                                         'x_s14dq6',
                                         'x_s14cq6',
                                         'x_s10dq7',
                                         'x_s14dq7'))
hajj = ave_dependent('belief_orthodox', c('x_s14aq1', 'x_s14aq3','x_s14aq4','x_s14aq5', 
                                            'x_s14aq6',
                                            'x_s14bq2', 
                                            'x_s14aq8',
                                            'x_s14aq9',
                                            'x_s14aq12',
                                            'x_s14aq13'))

hajj = ave_dependent('participate_local', c('x_2_s14aq16b',
                                              'x_2_s14aq16a',
                                              'x_2_s14aq16d',
                                              'x_2_s14aq16f'))

hajj = ave_dependent('view_other_countries', c('x_s15cq5', 
                                                 'x_s15cq6',
                                                 'x_s15cq7',
                                                 'x_s15cq8',
                                                 'x_s15cq10',
                                                 'x_s15cq9' ))

hajj = ave_dependent('view_other_group', c('x_s15b_b6',  
                                             'x_s15b_c6',
                                             'x_s15b_a6'))

hajj = ave_dependent('harmony', c('x_s15cq2',
                                    'x_s15cq3',
                                    'x_s15cq1',
                                    'x_s14aq7'))

hajj = ave_dependent('peaceful', c('x_s10bq4',
                                     'x_s10bq5',
                                     'x_s10cq5', 
                                     'x_s10cq4', 
                                     'x_s10cq3',
                                     'x_s10hq6',
                                     'x_s10hq7',
                                     'x_s10hq8'))

hajj =  ave_dependent('political_islam', c('x_s7q12a',
                                             'x_s7q12f',
                                             'x_s7q12d',
                                             'x_s7q12c',
                                             'x_s7q11c'))

hajj = ave_dependent('view_west', c('x_s10aq3',
                                      'x_s10aq2',
                                      'x_s10aq4',
                                      'x_s10aq5'))

hajj = ave_dependent('view_women', c('x_s10dq15a',
                                       'x_s10dq15b',
                                       'x_s10dq15c',
                                       'x_s10dq18'))

hajj = ave_dependent('women_QoL', c('xda_s10dq19b', 
                                      'xda_s10dq19c', 
                                      'xda_s10dq19d', 
                                      'x_2_s10dq8', 
                                      'xd_s10dq8'))

hajj = ave_dependent('girl_education', c('x_s10eq2',
                                           'x_s10eq5', 
                                           'x_s10eq4', 
                                           'x_s10eq6',  
                                           'x_s10eq3' ))

hajj = ave_dependent('women_in_work', c('x_s1_3q3',
                                          'x_s1_3q4',
                                          'x_s1_3q5r4'))

hajj = ave_dependent('gender_authority', c('x_s10dq15d',
                                             'x_s10dq6',
                                             'x_s10hq3', 
                                             'x_s10hq11',
                                             'x_s10dq1',
                                             'x_s10dq14',
                                             'x_2_s14cq8_9'))

hajj = ave_dependent('K6', c('s3q1a', 
                               's3q1b',
                               's3q1c',
                               's3q1d',
                               's3q1e', 
                               's3q1f'))

hajj = ave_dependent('positive_feeling', c('x_s3q2a',
                                             'x_s3q2b', 
                                             'x_s3q2c',
                                             'x_s3q7',
                                             'x_s3q3'))

hajj = ave_dependent('life_satisfaction', c('x_s3q4', 
                                              'x_s3q5',
                                              'x_s3q6'))

hajj = ave_dependent('physical_health', c('x_s2q1', 
                                            'x_s2q2' ))

hajj = ave_dependent('socio_economic_engagement', c('x_s5_5q1',
                                                      'x_s5_5q2',
                                                      'x_s5_5q3',
                                                      'x_s5_5q4',
                                                      's5_2q1a',
                                                      's5_2q2a',
                                                      's5_2q3a',
                                                      's5_2q1b',
                                                      's5_2q2b',
                                                      's5_2q3b',
                                                      'x_s5_4q1a',
                                                      'x_s5_4q2a',
                                                      'x_s5_4q3a', 
                                                      'x_s12q2',
                                                      'x_s12q1'))

hajj = ave_dependent('engage_politics', c('x_s7q7',
                                            'x_s7q1',
                                            'x_s7q2',
                                            'x_s5_4q5a',
                                            'x_s5_4q4a',
                                            'x_s7q10',
                                            'x_s7q9'))

hajj = ave_dependent('islam_knowledge', c('x_pillars',
                                            'x_s14bq5', 
                                            'x_s14bq7',
                                            'x_s14cq2',
                                            'x_s14cq1', 
                                            'x_s14dq2', 
                                            'x_s14dq1', 
                                            'x_s14dq5', 
                                            'x_s14bq3', 
                                            'x_s14cq3'))

hajj = ave_dependent('diversity_knowledge', c('x_s14dq3', 
                                                'x_2_s14dq6', 
                                                'x_s14dq4'))

hajj = ave_dependent('gender_knowledge', c('x_s14bq4',
                                             'x_s14bq6',
                                             'x_s14cq4',  
                                             'x_2_s14cq6', 
                                             'x_s10dq12',
                                             'xop_s10dq19b',
                                             'xop_s10dq19c',
                                             'xop_s10dq19d'))

hajj = ave_dependent('global_knowledge', c('x_s8q2',
                                             'x_s8q3', 
                                             'x_s8q7',
                                             'x_s8q8',
                                             'x_s8q6',
                                             'x_s8q5'))
#list(religious,belief_local,belief_orthodox, participate_local,view_other_countries,view_other_group,harmony,peaceful,political_islam,view_west,view_women,women_QoL,girl_education, women_in_work, gender_authority, K6, positive_feeling, life_satisfaction, physical_health, socio_economic_engagement, engage_politics, islam_knowledge, diversity_knowledge, gender_knowledge, gender_knowledge)
#hajj$district = as.factor(hajj$district)
library(nlme)
#categorical variables: cat, district, pod, ptygrp
#possible_covars = c('age', 'cat', 'district', 'female', 'literate', 'pod', 'ptygrp')



pod_encoder = onehot(hajj['pod'])
pod = predict(pod_encoder, hajj['pod'])
ptygrp_encoder = onehot(hajj['ptygrp'])
cat_encoder = onehot(hajj['cat'])
ptygrp = predict(ptygrp_encoder, hajj['ptygrp'])
cat = predict(cat_encoder, hajj['cat'])
district_encoder = onehot(hajj['district'])
district = predict(district_encoder, hajj['district'])

format_onehot = function(var){
  cols = colnames(var)
  for(i in seq(length(cols))){
    cols[i] = gsub('=', '',cols[i])
  }
  colnames(var) = cols
  return(var)  
}

district = format_onehot(district)
pod = format_onehot(pod)
ptygrp = format_onehot(ptygrp)
cat = format_onehot(cat)

hajj=merge(hajj, pod, by=0)[-1]
hajj=merge(hajj, cat, by=0)[-1]
hajj=merge(hajj, ptygrp, by=0)[-1]
hajj=merge(hajj, district, by=0)[-1]

dependent_vars = c("ave_belief_local", "ave_belief_orthodox", "ave_diversity_knowledge",
             "ave_engage_politics", "ave_gender_authority", "ave_gender_knowledge", "ave_girl_education",
             "ave_global_knowledge","ave_harmony","ave_islam_knowledge","ave_K6","ave_life_satisfaction", 
             "ave_participate_local", "ave_peaceful" , "ave_physical_health", "ave_political_islam", "ave_positive_feeling",         
             "ave_religious","ave_socio_economic_engagement", "ave_view_other_countries","ave_view_other_group",
             "ave_view_west", "ave_view_women", "ave_women_in_work","ave_women_QoL")

library(glmnet)
library(mice)

x_names = c('age', 'cat1', 'cat2', 'districtATTOCK', 'districtCHAKWAL', 'districtFAISALABAD',
            'districtGUJRAT', 'districtISLAMABAD', 'districtJHELUM', 'districtMULTAN','pod2', 'pod3', 'pod6', 'pod8',
            'ptygrp1', 'ptygrp2', 'ptygrp3', 'ptygrp5', 'ptygrp6', 'female', 'literate', 'urban')
xs = complete(mice(hajj[,x_names])) #turned out to be unnecessary as there was no missing value
hajj[,x_names] = xs

#try getting three way interactions
f = as.formula("success~.*.")
x_int = model.matrix(f, data=hajj[,c(x_names, 'success')])
x_int = x_int[,2:dim(x_int)[2]]
hajj=merge(hajj, x_int[,!colnames(x_int) %in% x_names], by=0)[-1]

success_selection = cv.glmnet(as.matrix(x_int), as.matrix(hajj$success))
cs = coef(success_selection, s='lambda.1se')
success_covars = colnames(x_int)[cs[,1]!=0] #lambda.max = .045

hajj_selection = cv.glmnet(as.matrix(x_int), as.matrix(hajj$hajj2006))
cs = coef(hajj_selection, s='lambda.1se')
hajj_covars = colnames(x_int)[cs[2:length(cs),1]!=0] #lambda.max = .036

xs = as.matrix(xs)
make_AES2_dic = function(dependent_var, all_measures=AES2, model_type=double_selection_model, 
                         num_terms=37){
  result = model_type(dependent_var)
  num_terms = dim(result)[1]
  all_measures[[dependent_var]][['beta']] = result[2]
  all_measures[[dependent_var]][['std_err']] = result[2+num_terms]
  all_measures[[dependent_var]][['t']] = result[2+2*num_terms]
  all_measures[[dependent_var]][['p']] = result[2+3*num_terms]
  return(all_measures)
}

double_selection_model = function(dependent_var){
  ys = hajj[!is.na(hajj[,dependent_var]), dependent_var]
  dependent_selection = cv.glmnet(as.matrix(x_int)[!is.na(hajj[,dependent_var]),], as.matrix(ys))
  cs = coef(dependent_selection, s='lambda.1se')
  dependent_covars = colnames(x_int)[cs[2:length(cs),1]!=0]
  covars_list = unique(c(success_covars, hajj_covars, dependent_covars))
  covars = paste(covars_list, collapse='+')
  hajj2006_right = paste(cbind("hajj2006", covars), collapse='+')
  success_right = paste(cbind("success", covars), collapse='+')
  formula1 = formula(paste(cbind(dependent_var, hajj2006_right), collapse='~'))
  formula2 = formula(paste(cbind('~', success_right), collapse=''))
  
  model=ivreg(formula1, formula2, data=hajj, x=TRUE)
  nonmissing = rep(TRUE, dim(hajj)[1])
  nonmissing[unname(model$na.action)] = FALSE
  model_result = cluster.robust.se(model, clusterid=hajj$clusterid[nonmissing])
  return(model_result)
}

double_tuned_AES = list()
for (var in dependent_vars){
  double_tuned_AES = make_AES2_dic(var, all_measures=double_tuned_AES)
}

double_selection_covars = function(dependent_var){
  ys = hajj[!is.na(hajj[,dependent_var]), dependent_var]
  dependent_selection = cv.glmnet(xs[!is.na(hajj[,dependent_var]),], as.matrix(ys))
  cs = coef(dependent_selection, s='lambda.1se')
  dependent_covars = colnames(x_int)[cs[2:length(cs),1]!=0] #lambda max range from .021 to .059
  covars_list = unique(c(success_covars, hajj_covars, dependent_covars)) #default behavior of lamba.min is  0.0001*lamba.max
  covars = paste(covars_list, collapse='+')
  return(covars)
}

covars = list()
for (var in dependent_vars){
  covars[[var]] = double_selection_covars(var)
}

double_selection_lambdas = function(dependent_var){
  lambdas = sort(seq(.04, .01, length.out=10), decreasing=TRUE)
  lambda_df = data.frame(todrop = rep(NA, 10), row.names = lambdas)
  success_selection = cv.glmnet(as.matrix(x_int), as.matrix(hajj$success), lambda = lambdas) #lambda.max = .045
  hajj_selection = cv.glmnet(as.matrix(x_int), as.matrix(hajj$hajj2006), lambda = lambdas)
  ys = hajj[!is.na(hajj[,dependent_var]), dependent_var]
  dependent_selection = cv.glmnet(xs[!is.na(hajj[,dependent_var]),], as.matrix(ys), lambda=lambdas)
  betas = rep(0, 10)
  std_errs = rep(0, 10)
  ts = rep(0, 10)
  ps = rep(0, 10)
  for( i in seq(10)){
    lambda = lambdas[i]
    cs = coef(success_selection, s=lambda)
    success_covars = colnames(x_int)[cs[,1]!=0]
    cs = coef(hajj_selection, s=lambda)
    hajj_covars = colnames(x_int)[cs[2:length(cs),1]!=0]
    cs = coef(dependent_selection, s=lambda)
    dependent_covars = colnames(x_int)[cs[2:length(cs),1]!=0] #lambda max range from .021 to .059
    covars_list = unique(c(success_covars, hajj_covars, dependent_covars))
    covars = paste(covars_list, collapse='+')
    hajj2006_right = paste(cbind("hajj2006", covars), collapse='+')
    success_right = paste(cbind("success", covars), collapse='+')
    formula1 = formula(paste(cbind(dependent_var, hajj2006_right), collapse='~'))
    formula2 = formula(paste(cbind('~', success_right), collapse=''))
    
    model=ivreg(formula1, formula2, data=hajj, x=TRUE)
    nonmissing = rep(TRUE, dim(hajj)[1])
    nonmissing[unname(model$na.action)] = FALSE
    result = cluster.robust.se(model, clusterid=hajj$clusterid[nonmissing])
    num_terms = dim(result)[1]
    betas[i] = result[2]
    std_errs[i] = result[2+num_terms]
    ts[i] = result[2+2*num_terms]
    ps[i] = result[2+3*num_terms]
  }
  return(data.frame(lambdas, betas, std_errs, ts, ps))
}

orthodox_df = double_selection_lambdas('ave_belief_orthodox')
diversity_df = double_selection_lambdas('ave_diversity_knowledge')
peace_df = double_selection_lambdas('ave_peaceful')
physical_df = double_selection_lambdas('ave_physical_health')
women_df = double_selection_lambdas('ave_view_women')

double_tuned_AES = data.frame(double_tuned_AES)
Holm_double_tuned_ps = p.adjust(double_tuned_AES['p',], method='holm')
dependent_vars[Holm_double_tuned_ps<.05]
BH_double_tuned_ps = p.adjust(double_tuned_AES['p',], method='BH')
dependent_vars[BH_double_tuned_ps<.05]

##############
#The below version was run using lambda.min for cv.glmnet instead of lambda.1se. Not recommended. 
# AES2 = list()
# for (var in dependent_vars){
#   AES2 = make_AES2_dic(var)
# }
# AES2 = data.frame(AES2)
# holm_AES2 = p.adjust(AES2['p',], method='holm', n=25)
# BH_AES2 = p.adjust(AES2['p',], method='BH', n=25)
# print(dependent_vars[BH_AES2<.05])
# print(dependent_vars[holm_AES2<.05])
# print(dependent_vars[AES2['p',]<.05])
# print(dependent_vars[BY_ps<.05])
# print(dependent_vars[BH_ps<.05])
# unique(cbind(dependent_vars[AES2['p',]<.05], dependent_vars[AES2['p',]<.05]))


all_covar_model = function(dependent_var){
  covars = x_names
  hajj2006_right = paste(cbind("hajj2006", covars), collapse='+')
  success_right = paste(cbind("success", covars), collapse='+')
  formula1 = formula(paste(cbind(dependent_var, hajj2006_right), collapse='~'))
  formula2 = formula(paste(cbind('~', success_right), collapse=''))
  
  model=ivreg(formula1, formula2, data=hajj, x=TRUE)
  nonmissing = rep(TRUE, dim(hajj)[1])
  nonmissing[unname(model$na.action)] = FALSE
  model_result = cluster.robust.se(model, clusterid=hajj$clusterid[nonmissing])
  return(model_result)
}
AES_all = list()
for (var in dependent_vars){
  AES_all = make_AES2_dic(var, all_measures=AES_all,model_type=all_covar_model)
}

AES_all = data.frame(AES_all)
holm_AES_all = p.adjust(AES_all['p',], method='holm', n=25)
BH_AES_all = p.adjust(AES_all['p',], method='BH', n=25)
print(dependent_vars[BH_AES_all<.05])
print(dependent_vars[holm_AES_all<.05])
print(dependent_vars[AES2['p',]<.05])
print(dependent_vars[BY_ps<.05])
print(dependent_vars[BH_ps<.05])
unique(cbind(dependent_vars[AES2['p',]<.05], dependent_vars[AES2['p',]<.05]))

example = lm(hajj2006~female, data=hajj)
library('grf')
forest = instrumental_forest(xs, Y=hajj$ave_belief_local, W=hajj$hajj2006, Z=hajj$success)
effect = average_treatment_effect(forest)

options(digits=5)
AES_results = data.frame()
for (col in names(AES_dic)){
  AES_results[col,'beta'] = as.numeric(AES_dic[[col]]['beta'])
  AES_results[col,'std_err'] = as.numeric(AES_dic[[col]]['std_err'])
  AES_results[col,'t'] =  as.numeric(AES_dic[[col]]['t'])
  AES_results[col,'p'] =  as.numeric(AES_dic[[col]]['p'])
}
AES_results[,'FWER corrected p'] = p.adjust(AES_results[,'p'], method='holm', n=25)
AES_results[,'FDR corrected p'] = p.adjust(AES_results[,'p'], method='BH', n=25)

double_results = data.frame()
for (col in names(double_tuned_AES)){
  double_results[col,'beta'] = as.numeric(double_tuned_AES[[col]]['beta'])
  double_results[col,'std_err'] = as.numeric(double_tuned_AES[[col]]['std_err'])
  double_results[col,'t'] =  as.numeric(double_tuned_AES[[col]]['t'])
  double_results[col,'p'] =  as.numeric(double_tuned_AES[[col]]['p'])
}
library(xtable)
double_results[,'FWER corrected p'] = p.adjust(double_results[,'p'], method='holm', n=25)
double_results[,'FDR corrected p'] = p.adjust(double_results[,'p'], method='BH', n=25)
for (col in names(double_tuned_AES)){
double_results[col, 'covariates'] =  covars[[col]]
}
xtmp = xtable(double_results, digits=c(0, 3, 3, 3, 5, 5, 5, 0))

AES_results = AES_results[sort(row.names(AES_results)), ]
result_summary = data.frame(todrop = rep(NA, 25), row.names = row.names(AES_results))
result_summary[,'significance from regular IV regression'] = AES_results$p<.05
result_summary[,'IV regression, FWER corrected'] = AES_results$`FWER corrected p`<.05
result_summary[,'IV regression, FDR corrected'] = AES_results$`FDR corrected p`<.05
result_summary[,'significance with post double selection'] = double_results$p<.05
result_summary[,'post-double-selection, FWER corrected'] = double_results$`FWER corrected p`<.05
result_summary[,'post-double-selection, FDR corrected'] =double_results$`FDR corrected p`<.05
result_summary$todrop = NULL
View(result_summary)
write.csv(result_summary, "~/Desktop/result_summary2.csv")


contrast_methods = function(dependent_name, double_df){
  par(mar=c(6,4,4,2), las=1)
  lambdas = seq(.04, .01, length.out=10)
  low = min(AES_results[dependent_name, 'beta']- AES_results[dependent_name, 'std_err'], min(double_df[,'betas']-double_df[,'std_errs']))
  high = max(AES_results[dependent_name, 'beta']+AES_results[dependent_name, 'std_err'], max(double_df[,'betas']+double_df[,'std_errs']))
  if(low>high){
    x = low
    low=high
    high=x
  }
  plot(rev(lambdas), rev(double_df[,'betas']), type='n',xlab='lambda (L1 penalty parameter) value', ylab='Effect Size', 
       ylim = c(low, high), axes=FALSE)
  lines(rev(lambdas), rev(double_df[,'betas']), lty=1, col='blue')
  lines(rev(lambdas), rev(double_df[,'betas']-double_df[,'std_errs']), lty=2, col='blue')
  lines(rev(lambdas), rev(double_df[,'betas']+double_df[,'std_errs']), lty=2, col='blue')
  lines(rev(lambdas), rep(AES_results[dependent_name, 'beta'], 10), lty=1, col='red')
  lines(rev(lambdas), rep(AES_results[dependent_name, 'beta']- AES_results[dependent_name, 'std_err'], 10), lty=2, col='red')
  lines(rev(lambdas), rep(AES_results[dependent_name, 'beta']+ AES_results[dependent_name, 'std_err'], 10), lty=2, col='red')
  legend(.01, low+.02, c('regular IV regression', 'post-double-selection'), lty=c(1, 1), col=c('red', 'blue'), cex=0.6)
  axis(1, at=rev(lambdas), labels=lapply(rev(lambdas), FUN=(function(x) substr(as.character(x), start=0, stop=4))))
  ys = seq( low, high, length.out=10)
  ys = lapply(ys, FUN=(function(x) substr(as.character(x), start=0, stop=5)))
  axis(2, at=ys, labels=ys)
  title(paste('Effect size comparison between two methods:', dependent_name, collapse=' '))
  #sprintf('The pooled ATEs for each firms are:')
  #print(ATEs)
}
contrast_methods('view_women', women_df)
contrast_methods('belief_orthodox', orthodox_df)
contrast_methods('diversity_knowledge', diversity_df )
contrast_methods('peaceful', peace_df)
contrast_methods('physical_health', physical_df)

