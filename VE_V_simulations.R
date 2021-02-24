require(tidyverse)
require(splitstackshape)

# parameters
parms <- list(trial_length = 300, # end of simulations (so trial actually 200 days long since starting on day 100)
              enrollment_day = 100, # day enroll into trial
              N = 100000, # number of people
              VES = 0.6, # vaccine efficacy against susceptibility
              VED = 0.6, # vaccine efficacy against duration
              VEP = 0.9, # vaccine efficacy against progression to symptoms
              prob_sympt = 0.8, # probability unvaccinated symptomatic
              lat_pd = 3, # latent period (viral load increases after this)
              inf_pd = seq(15,21,1), # for how long viral shedding lasts
              inc_pd = 5, # days symptom onset (incubation period)
              model = "SIS", # options are SIS or SIR
              samp_freq = 1, # how many days looking across
              protocol = "Moderna", # only considering Moderna style (one sample)
              LOD = 0, # no limit of detection for viral load
              X = NA) # relevant for protocols besides Moderna (not done here)

# external force of infection
beta <- 0.001  # 0.001 #0.005
extFOI_flat <- rep(beta,parms[["trial_length"]]) # flat to start

# sampling and VE calculations 
sampling <- function(inf_master,parms){
  
  if(parms[["protocol"]]=="Moderna"){
    samp_freq  <- 1 # looking at every day after enrollment (treating each day independently)
    sample_days <- seq(parms[["enrollment_day"]] + samp_freq,parms[["trial_length"]],samp_freq)
  } 
  
  # VL 
  # Note this is set up to model individual VL trajectory but the actual trajectories are NOT used in the paper as our LOD is 0
  # Leaving in case future iterations want to use but all that matters is tracking which days VL_day > 0 which this does
  inf_master %>%
    ungroup() %>%
    mutate(length = DayRecovered - DayInfectious + 2, # adding one on each end so can have a 0 to start and 0 on day recovered
           VL_rise = round(runif(nrow(inf_master), 1, 4)), # number of days increase
           VL = abs(rnorm(runif(nrow(inf_master)), 8, 1))/(VL_rise+1), # starting VL 
           VL_waning = (VL * 1.95^(VL_rise-1))/(length - 1 - VL_rise)) %>% # how much decrease each day
    expandRows("length") %>%
    group_by(ID,DayExposed) %>%
    mutate(Inf.days = row_number(), # count days since infection
           VL_day = case_when(Inf.days == 1 ~ 0, # make first day zero VL (representing last day of latent period)
                              Inf.days == 2 ~ VL,
                              Inf.days <=(VL_rise+1) ~ VL * 1.95^(Inf.days-2), # increase VL for VL_rise days
                              Inf.days > (VL_rise+1)  ~ VL * 1.95^(VL_rise-1) - (Inf.days-VL_rise-1)*VL_waning), # waning
           VL_day = case_when(VL_day<0 ~ 0, # remove negative values (make zero)
                              TRUE ~ VL_day)) -> VLs
  
  VLs %>% # if rise was greater than duration, last one will not be 0 as it should so manually correct 
            #(other ways to do this by updating VL_rise but just check for now)
    group_by(ID,DayExposed) %>%
    mutate(max_day = max(Inf.days),
           VL_day = case_when(Inf.days == max_day ~ 0,
                              TRUE ~ VL_day)) -> VLs
  
  # sampling
  VLs %>%
    mutate(t = DayInfectious + Inf.days - 2, # convert to actual time
           test_day = case_when(t %in% sample_days ~ 1, # determine if each day would be a test; if sample_freq is 1, 0s are only after study period
                                TRUE ~ 0),
           test_result = case_when(test_day == 1 & VL_day > parms[["LOD"]] ~ 1, # if test day, result based on threshold (here is 0)
                                   TRUE ~ 0))-> VLs_test
    
  # Swab 
  
  VLs_test %>%
    subset(t %in% sample_days) %>%
    group_by(t, VE) %>%
    summarise(npos = sum(test_result)) %>%
    mutate(nneg=parms[["N"]]/2 - npos) -> results # how many people test positive each day grouped by vaccine status
  
  results %>%
    pivot_wider(names_from = VE, values_from = c(npos,nneg)) %>%
    mutate(npos_0 = ifelse(is.na(npos_0),0,npos_0),
           npos_1 = ifelse(is.na(npos_1),0,npos_1),
           nneg_1 = ifelse(is.na(nneg_1),parms[["N"]]/2 ,nneg_1),
           nneg_0 = ifelse(is.na(nneg_0),parms[["N"]]/2 ,nneg_0)) %>%
    subset(!is.na(t)) %>%
    mutate(VEprev = 1 - (npos_1/nneg_1)/(npos_0/nneg_0),
           check = npos_0 + npos_1 + nneg_0 + nneg_1) -> VEprev # calculate VEprev from 1-OR vax/unvax
  
  # Non-symptomatic
  
  VLs_test %>%
    subset(t %in% sample_days & (DaySymptoms>t | is.na(DaySymptoms))) %>% # look only among those that are not yet or never symptomatic who would test positive
    group_by(t, VE) %>%
    summarise(npos = sum(test_result)) -> num_pos_a
  
  results_a <- NULL # need to get number who would not have been tested (had symptoms between enrollment and day)
  for (day in (parms[["enrollment_day"]]+1):(parms[["trial_length"]])){
    VLs_test %>%
      ungroup() %>%
      subset(DaySymptoms > parms[["enrollment_day"]] & DaySymptoms <= day & !is.na(DaySymptoms)) %>%
      dplyr::select(ID,VE) %>%
      unique() %>% # deduplicate
      group_by(VE) %>%
      summarise(n_exclude = n()) %>% # this is the number to exclude from testing
      mutate(t=day) %>%
      bind_rows(results_a) -> results_a# get all who've had symptoms through today
  }
  
  results_a %>%
    full_join(num_pos_a,by=c("t","VE")) %>% # merge back with non-symptomatic testing positive
    mutate(n_exclude = ifelse(is.na(n_exclude),0,n_exclude)) %>%
    mutate(nneg = parms[["N"]]/2 - npos - n_exclude) %>% # negative tests are total - those testing positive - those not tested (symptoms)
    dplyr::select(-n_exclude) %>%
    pivot_wider(names_from = VE, values_from = c(npos,nneg)) %>%
    subset(!is.na(t)) %>%
    mutate(VEprev_a = 1 - (npos_1/nneg_1)/(npos_0/nneg_0),
           check = npos_0 + npos_1 + nneg_0 + nneg_1) -> VEprev_a
  
  # All
  results_all <- NULL # look at combined symptom and tests
  for (day in (parms[["enrollment_day"]]+1):(parms[["trial_length"]])){ # loop through days since enrollment
    VLs_test %>%
      ungroup() %>%
      subset(t==day & test_result == 1) %>%
      dplyr::select(ID,VE) -> tests # get those that test positive this day
    
    VLs_test %>%
      ungroup() %>%
      subset(DaySymptoms > parms[["enrollment_day"]] & DaySymptoms <= day) %>% # get all who've had symptoms through today
      dplyr::select(ID,VE) %>%
      bind_rows(tests) %>% # combine with positive tests
      unique() %>% # deduplicate
      group_by(VE) %>%
      summarise(npos = n()) %>% # summarize
      mutate(t=day) %>%
      mutate(nneg=parms[["N"]]/2 - npos) %>%
      bind_rows(results_all) %>%
      dplyr::select(t,VE,npos,nneg) -> results_all
  }
  
  results_all %>%
    pivot_wider(names_from = VE, values_from = c(npos,nneg)) %>%
    mutate(npos_0 = ifelse(is.na(npos_0),0,npos_0),
           npos_1 = ifelse(is.na(npos_1),0,npos_1),
           nneg_1 = ifelse(is.na(nneg_1),parms[["N"]]/2 ,nneg_1),
           nneg_0 = ifelse(is.na(nneg_0),parms[["N"]]/2 ,nneg_0)) %>%
    subset(!is.na(t)) %>%
    mutate(VEall = 1 - (npos_1/nneg_1)/(npos_0/nneg_0),
           check = npos_0 + npos_1 + nneg_0 + nneg_1) -> VEcombined
  
  if (unique(VEprev$check) != parms[["N"]]){print("Error")}
  if (unique(VEcombined$check) != parms[["N"]]){print("Error")}
  
  list(VEprev,VEcombined,VEprev_a,VLs_test)
}

# loop through parameters
nsims <- 1
protocol <- "Moderna"
extFOI_list <- list(extFOI_flat) 
Infection_master_file <- NULL
results_master_file_VEprev <- NULL
results_master_file_VEprev_a <- NULL
results_master_file_VEall <- NULL
for (model in c("SIS","SIR")){ #,
  for (j in 1:length(extFOI_list)){ 
    for (VES in seq(0,1,0.3)){ #
      for (VED in seq(0,1,0.3)){ #
        for (sympt in c(0.01,0.8)){ #0.01
          
          extFOI <- extFOI_list[[j]]
        
          parms[["model"]] <- model  
          parms[["VES"]] <- VES   
          parms[["VED"]] <- VED  
          parms[["VEP"]] <- 1-0.05/(1-VES) # VESP = 0.95 = 1-(1-VES)*(1-VEP)
          parms[["prob_sympt"]] <- sympt
          VEP1 <- 1-(1-0.95/2)/(1-VES/2) # first dose assumed half efficacy
  
          check <- 1-(1-VES)*(1-parms[["VEP"]]) #0.95
          check2 <- 1-(1-VES/2)*(1-VEP1) #0.475
          
          cat(j,model,VES,VED,check,check2,sympt,"\n")
          
          # matrix for force of infection for everyone (depends on day relative to vaccination and # doses); only half people enrolled
          extFOI_mat <- matrix(extFOI, nrow = (parms[["N"]]), ncol = parms[["trial_length"]],byrow=TRUE)
          extFOI_mat[1:(parms[["N"]]/2),(parms[["enrollment_day"]]+7):(parms[["enrollment_day"]]+34)]<- 
            extFOI_mat[1:(parms[["N"]]/2),(parms[["enrollment_day"]]+7):(parms[["enrollment_day"]]+34)]*(1-parms[["VES"]]/2)
          extFOI_mat[1:(parms[["N"]]/2),(parms[["enrollment_day"]]+35):ncol(extFOI_mat)]<- 
            extFOI_mat[1:(parms[["N"]]/2),(parms[["enrollment_day"]]+35):ncol(extFOI_mat)]*(1-parms[["VES"]])
          
          # conduct bernoulli trials for each person each day to see who/when infected
          inf_prob <- matrix(rbinom(nrow(extFOI_mat)*ncol(extFOI_mat),1,prob=extFOI_mat),nrow = (parms[["N"]]), ncol = parms[["trial_length"]])
           
          dates_exposed <- apply(inf_prob,1,function(x){which(x==1)})
          
          inf_master <- NULL # create list of all infections
          for (i in 1:length(dates_exposed)){
            if(sum(dates_exposed[[i]])>0){
              bind_cols(ID=i,
                        VE = ifelse(i<=(parms[["N"]]/2),1,0),
                        DayExposed = dates_exposed[[i]]) %>%
                bind_rows(inf_master) -> inf_master
            }
          }
          
          inf_master %>%
            subset(DayExposed < (parms[["enrollment_day"]]+7)) -> inf_master_1 # split into 3 relative to vaccine days
          
          inf_master %>%
            subset(DayExposed >= (parms[["enrollment_day"]]+7) & DayExposed < (parms[["enrollment_day"]]+35)) -> inf_master_2
          
          inf_master %>%
            subset(DayExposed >= (parms[["enrollment_day"]]+35)) -> inf_master_3
          
          inf_master_1 %>%
            mutate(DayInfectious = DayExposed + parms[["lat_pd"]],
                   Symptoms = rbinom(nrow(inf_master_1),1,parms[["prob_sympt"]]),
                   DaySymptoms = case_when(Symptoms==1 ~ DayExposed + parms[["inc_pd"]]),
                   DayRecovered = round(DayInfectious + sample(parms[["inf_pd"]],nrow(inf_master_1),replace=TRUE))) -> inf_master_1
          
          inf_master_2 %>%
            mutate(DayInfectious = DayExposed + parms[["lat_pd"]],
                   Symptoms = rbinom(nrow(inf_master_2),1,parms[["prob_sympt"]]*(1-VEP1*VE)), # reduce probability of symptoms for vaccinated
                   DaySymptoms = case_when(Symptoms==1 ~ DayExposed + parms[["inc_pd"]]),
                   DayRecovered = round(DayInfectious + sample(parms[["inf_pd"]],nrow(inf_master_2),replace=TRUE)*(1-parms[["VED"]]/2*VE))) -> inf_master_2
                    # above line reduces duration of infection for vaccinated
          
          inf_master_3 %>%
            mutate(DayInfectious = DayExposed + parms[["lat_pd"]],
                   Symptoms = rbinom(nrow(inf_master_3),1,parms[["prob_sympt"]]*(1-parms[["VEP"]]*VE)),
                   DaySymptoms = case_when(Symptoms==1 ~ DayExposed + parms[["inc_pd"]]),
                   DayRecovered = round(DayInfectious + sample(parms[["inf_pd"]],nrow(inf_master_3),replace=TRUE)*(1-parms[["VED"]]*VE))) %>%
            bind_rows(inf_master_1,inf_master_2) -> inf_master
          
          if (parms[["model"]]=="SIR"){ # only keep first time infected for SIR
            inf_master %>%
              group_by(ID) %>%
              slice_min(DayExposed,n=1) -> inf_master
          } else if (parms[["model"]]=="SIS"){ # can get reinfected but only after recover
            inf_master %>%
              arrange(DayExposed) %>%
              group_by(ID) %>%
              mutate(numrow = row_number()) -> inf_master_rows
            
            for (row in 2:max(inf_master_rows$numrow)){
              inf_master_rows %>%
                mutate(eligible_day = case_when(numrow <= row ~ lag(DayRecovered))) %>%
                subset(DayExposed >= eligible_day | is.na(eligible_day)) -> inf_master_rows
              
            }
            
            inf_master_rows -> inf_master
            
          }
          
          X <- NA # relevant for other types of protocols (not done here)
          sampling_results <- sampling(inf_master,parms)
          
          VEprev <- sampling_results[[1]]
          VEcombined <- sampling_results[[2]]
          VEprev_a <- sampling_results[[3]]
          Tests <- sampling_results[[4]]
          
          VEprev %>%
            bind_cols(ves=VES,
                      ved=VED,
                      sympt = parms[["prob_sympt"]],
                      X = X,
                      model = model,
                      protocol = protocol,
                      foi = j) %>%
            bind_rows(results_master_file_VEprev)-> results_master_file_VEprev
          
          
          VEprev_a %>%
            bind_cols(ves=VES,
                      ved=VED,
                      sympt = parms[["prob_sympt"]],
                      X = X,
                      model = model,
                      protocol = protocol,
                      foi = j) %>%
            bind_rows(results_master_file_VEprev_a)-> results_master_file_VEprev_a
          
          VEcombined %>%
            bind_cols(ves=VES,
                      ved=VED,
                      X = X,
                      model = model,
                      protocol = protocol,
                      sympt = parms[["prob_sympt"]] ,
                      foi = j) %>%
            bind_rows(results_master_file_VEall)-> results_master_file_VEall
          
          inf_master %>%
            bind_cols(ves=VES,
                      ved=VED,
                      X = X,
                      model = model,
                      protocol = protocol,
                      sympt = parms[["prob_sympt"]] ,
                      foi = j) %>%
            bind_rows(Infection_master_file)-> Infection_master_file
        }
      }
    }
  }
}


write.csv(results_master_file_VEprev,paste0(extFOI_flat[1],"_results_master_file_VEprev.csv"))
write.csv(results_master_file_VEprev_a,paste0(extFOI_flat[1],"_results_master_file_VEprev_a.csv"))
write.csv(results_master_file_VEall,paste0(extFOI_flat[1],"_results_master_file_VEall.csv"))
write.csv(Infection_master_file,paste0(extFOI_flat[1],"_Infection_master_file.csv"))




