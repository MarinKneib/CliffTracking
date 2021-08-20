PopDynMod_MainLoop=function(param_list,n_sim,Init_pop,N_yrs){
  
  ## This function is the main loop of the cliff population dynamics model
  # 
  # INPUTS
  #  param_list: list of 46 parameters describing cliff population characteristics
  #  n_sim: number of simulations required
  #  Init_pop: Vector with size of all initial cliffs
  #  N_yrs: number of years over which the simulations are run
  # OUTPUTS
  #  tot_sim_nb: Number of cliffs at each time-step for each run
  #  tot_sim_sz: Total size of cliffs at each time-step for each run
  # 
  # Author: Marin Kneib
  # Work address: Swiss Federal Research Institute WSL
  # Email: marin.kneib@wsl.ch
  
  tot_sim_nb<-c()
  tot_sim_sz<-c()
  
  for (p in 1:n_sim){
    print(p)
    Cliff_pop<-Init_pop
    Cumul_nb<-length(Cliff_pop)
    Cumul_area<-0
    for (i in 1:length(Cliff_pop)){
      Cumul_area<-Cumul_area+Cliff_pop[[i]][1]
    }
    All_cumul_nb<-c(Cumul_nb)
    All_cumul_area<-c(Cumul_area)
    for (n in 1:N_yrs){
      
      # KILL: randomly choosen cliffs based on death rate
      N_kill<-(param_list[1]*length(Cliff_pop)+param_list[2])*(1+rnorm(1,param_list[4],param_list[3]))  # Number of cliffs to kill
      if (N_kill>0){
        for (m in 1:round(N_kill)){
          if (length(Cliff_pop)>0){ # if all the cliffs are dead, not possible to kill more...
            idSize<-rlnorm(1,param_list[5],param_list[6]) # Choose size of cliff to kill based on lognormal distribution
            # Find cliff with size most similar to ideal size
            idx_kill<-1
            Size_diff_kill<-abs(Cliff_pop[[idx_kill]][1]-idSize)
            for (i in 1:length(Cliff_pop)){
              Size_diff<-abs(Cliff_pop[[i]][1]-idSize)   # Size difference of all cliffs with ideal size
              if (Size_diff<Size_diff_kill){
                idx_kill<-i
                Size_diff_kill<-Size_diff
              }
            }
            # delete selected cliff
            Cliff_pop<-Cliff_pop[-idx_kill]
          }
        }   
      }
      
      
      # SPLIT/MERGE: first select which cliffs are to be split/merged and which ones remain. Then apply merging/splitting/area increase.
      
      
      # SPLIT: randomly selected cliffs (N based on split events) are divided into other cliffs
      N_split<-round(rnorm(1,param_list[17],param_list[18]))
      Cliffs2split<-c()
      if (N_split>0){
        for (k in 1:N_split){
          if (length(Cliff_pop>0)){
            idSize<-rlnorm(1,param_list[25],param_list[26]) # Choose size of cliff to split based on lognormal distribution
            # Find cliff with size most similar to ideal size
            idx_split<-1
            Size_diff_split<-abs(Cliff_pop[[idx_split]][1]-idSize)
            for (i in 1:length(Cliff_pop)){
              Size_diff<-abs(Cliff_pop[[i]][1]-idSize)   # Size difference of all cliffs with ideal size
              if ((Size_diff<Size_diff_split) & (sum(Cliffs2split==i) == 0)){     # Make sure that the nearest cliff hasn't already been choosen
                idx_split<-i
                Size_diff_split<-Size_diff
              }
            }
            Cliffs2split<-c(Cliffs2split,idx_split)
          }
        }
      }
      
      # MERGE: randomly selected cliffs (N based on merge events) are combined into one
      N_merge_evts<-length(Cliff_pop)+1
      while (length(Cliffs2split)+2.5*N_merge_evts >= length(Cliff_pop)){
        N_merge_evts<-round(rnorm(1,param_list[9],param_list[10]))
      }
      Cliffs2merge<-c()
      if (N_merge_evts>0){
        N_merge<-0
        while ((N_merge<2*N_merge_evts) | (N_merge>=4*N_merge_evts)){
          N_merge<-round(N_merge_evts*(rnorm(1,param_list[11],param_list[12])))
        }
        for (k in 1:N_merge){
          # Find which cliffs merge (randomly choose from existing cliffs except the ones that split / already merge)
          idx_merge<-0
          while ((idx_merge == 0) | (sum(Cliffs2split==idx_merge) > 0) | (sum(Cliffs2merge==idx_merge) > 0)){
            idx_merge<-sample(1:length(Cliff_pop),1)
          }
          Cliffs2merge<-c(Cliffs2merge,idx_merge)
        }
      }
      
      # MIXED: cliffs that split & merge
      N_mixed_evts<-round(rnorm(1,param_list[27],param_list[28]))
      Cliffs2mix<-c()
      if (N_mixed_evts>0){
        N_mother_mixed<-0
        while ((N_mother_mixed<2*N_mixed_evts) | (N_mother_mixed>=4*N_mixed_evts)){
          N_mother_mixed<-round(N_mixed_evts*rnorm(1,param_list[29],param_list[30]))
        }
        for (k in 1:N_mother_mixed){
          if (length(Cliff_pop)>2){
            idSize<-rlnorm(1,param_list[37],param_list[38]) # Choose size of cliff to mix based on lognormal distribution
            # Find cliff with size most similar to ideal size
            idx_mix<-1
            Size_diff_mix<-abs(Cliff_pop[[idx_mix]][1]-idSize)
            for (i in 1:length(Cliff_pop)){
              Size_diff<-abs(Cliff_pop[[i]][1]-idSize)   # Size difference of all cliffs with ideal size
              if ((Size_diff<Size_diff_mix) & (sum(Cliffs2split==i) == 0) & (sum(Cliffs2merge==i) == 0) & (sum(Cliffs2mix==i) == 0)){     # Make sure that the nearest cliff hasn't already been choosen
                idx_mix<-i
                Size_diff_mix<-Size_diff
              }
            }
            Cliffs2mix<-c(Cliffs2mix,idx_mix)
          }
        }
      }
      
      # AREA CHANGES: for cliffs that do not split/merge
      Cliffs2grow<-c(1:length(Cliff_pop))
      if (length(Cliffs2split)>0){Cliffs2grow<-Cliffs2grow[-Cliffs2split]}
      if (length(Cliffs2merge)>0){Cliffs2grow<-Cliffs2grow[-Cliffs2merge]}
      if (length(Cliffs2mix)>0){Cliffs2grow<-Cliffs2grow[-Cliffs2mix]}
      
      # Change area of selected cliffs
      for (i in Cliffs2grow){
        Init_size<-Cliff_pop[[i]][1]
        Fin_size<-0
        # Calculate area change
        if (param_list[41]*log(Init_size)+param_list[42]>=0){
          while ((Fin_size<param_list[8]) | (Fin_size>param_list[7])){
            frac_grow<-rlnorm(1,param_list[39]*log(Init_size)+param_list[40],param_list[41]*log(Init_size)+param_list[42])
            Fin_size<-Init_size*frac_grow
          }
        }
        if (param_list[41]*log(Init_size)+param_list[42]<0){
          Fin_size<-Init_size
        }
        # Apply area change
        Cliff_pop[[i]][1]<-round(Fin_size)
      }
      
      # SPLIT selected cliffs
      
      # Define number of daughters for each split event
      if (N_split>0){
        Daughters_split<-rep(2,N_split) # Splitting cliffs have at least 2 daughters each
        N_daughters_split<-0
        while ((N_daughters_split<2*N_split) | (N_daughters_split>=4*N_split)){
          N_daughters_split<-round(N_split*rnorm(1,param_list[19],param_list[20]))
        }
        # Add daughters until there are no more 
        if ((N_daughters_split %/% N_split==2 & N_daughters_split/N_split>2) | (N_daughters_split %/% N_split==3 & N_daughters_split/N_split==3)){
          for (j in (1:(N_daughters_split-2*N_split))){
            Daughters_split[j]<-Daughters_split[j]+1
          }
        }
        if (N_daughters_split %/% N_split==3 & N_daughters_split/N_split>3){
          Daughters_split<-rep(3,N_split)
          for (j in (1:(N_daughters_split-3*N_split))){
            Daughters_split[j]<-Daughters_split[j]+1
          }
        }
        idx<-0
        for (i in Cliffs2split){
          idx<-idx+1
          Nb_daughters<-Daughters_split[idx]
          Init_size<-Cliff_pop[[i]][1]
          Fin_size<-0
          # Calculate area change 
          while ((Fin_size<param_list[8]*Nb_daughters) | (Fin_size>param_list[7]*Nb_daughters)){
            if (param_list[23]*log(Init_size)+param_list[24]>=0){
              frac_split<-rlnorm(1,param_list[21]*log(Init_size)+param_list[22],param_list[23]*log(Init_size)+param_list[24])
            }
            if (param_list[23]*log(Init_size)+param_list[24]<0){
              frac_split<-(param_list[8]+1)*Nb_daughters/Init_size
            }
            Fin_size<-Init_size*frac_split
          }
          # Distribute final area to all daughters
          Size_daughters<-rep(0,Nb_daughters)
          for (j in 1:Nb_daughters){
            if (j<Nb_daughters){
              Size_daughters[j]<-runif(1,param_list[8],Fin_size-sum(Size_daughters)-(Nb_daughters-j)*param_list[8])
            }
            if (j==Nb_daughters){
              Size_daughters[j]<-Fin_size-sum(Size_daughters)
            }
          }
          # Modify original cliff size
          Cliff_pop[[i]][1]<-round(Size_daughters[1])
          # Add new cliff of size given by above calculation
          for (j in 2:Nb_daughters){
            Cliff_pop[[length(Cliff_pop)+1]]<-round(Size_daughters[j])
          }
        }
      }
      
      # Cliffs to kill after merge/mixed events
      idx_kill<-c()
      
      # MERGE selected cliffs
      if (N_merge_evts>0){
        Mothers_merge<-rep(2,N_merge_evts)
        # Add mothers until there are no more 
        if ((N_merge %/% N_merge_evts==2 & N_merge/N_merge_evts>2) | (N_merge %/% N_merge_evts==3 & N_merge/N_merge_evts==3)){
          for (j in (1:(N_merge-2*N_merge_evts))){
            Mothers_merge[j]<-Mothers_merge[j]+1
          }
        }
        if (N_merge %/% N_merge_evts==3 & N_merge/N_merge_evts>3){
          Mothers_merge<-rep(3,N_merge_evts)
          for (j in (1:(N_merge-3*N_merge_evts))){
            Mothers_merge[j]<-Mothers_merge[j]+1
          }
        }
        for (i in 1:N_merge_evts){
          # Calculate merging area fraction (normal distribution)
          if (i==1){
            Init_cliffs<-Cliffs2merge[1:Mothers_merge[1]]
          }
          if (i>1){
            Init_cliffs<-Cliffs2merge[(1+sum(Mothers_merge[1:(i-1)])):sum(Mothers_merge[1:i])]
          }
          Init_size<-0
          for (j in Init_cliffs){
            Init_size<-Init_size+Cliff_pop[[j]][1]
          }
          Fin_size<-0
          # Calculate area change 
          while ((Fin_size<param_list[8]) | (Fin_size>param_list[7])){
            if (param_list[15]*log(Init_size)+param_list[16]>=0){
              frac_merge<-rlnorm(1,param_list[13]*log(Init_size)+param_list[14],param_list[15]*log(Init_size)+param_list[16])
            }
            Fin_size<-Init_size*frac_merge
          }
          # Change size of first cliff
          Cliff_pop[[Init_cliffs[1]]][1]<-Fin_size
          # Kill all other
          idx_kill<-c(idx_kill,Init_cliffs[2:length(Init_cliffs)])
        }
      }
      
      # MIXED events
      if (N_mixed_evts>0){
        N_daughter_mixed<-0
        while ((N_daughter_mixed<2*N_mixed_evts) | (N_daughter_mixed>=4*N_mixed_evts)){
          N_daughter_mixed<-round(N_mixed_evts*rnorm(1,param_list[31],param_list[32]))
        }
        Mothers_mixed<-rep(2,N_mixed_evts)
        Daughters_mixed<-rep(2,N_mixed_evts)
        # Add mothers & daughters until there are no more 
        if ((N_mother_mixed %/% N_mixed_evts==2 & N_mother_mixed/N_mixed_evts>2) | (N_mother_mixed %/% N_mixed_evts==3 & N_mother_mixed/N_mixed_evts==3)){
          for (j in (1:(N_mother_mixed-2*N_mixed_evts))){
            Mothers_mixed[j]<-Mothers_mixed[j]+1
          }
        }
        if (N_mother_mixed %/% N_mixed_evts==3 & N_mother_mixed/N_mixed_evts>3){
          Mothers_mixed<-rep(3,N_mixed_evts)
          for (j in (1:(N_mother_mixed-3*N_mixed_evts))){
            Mothers_mixed[j]<-Mothers_mixed[j]+1
          }
        }
        if ((N_daughter_mixed %/% N_mixed_evts==2 & N_daughter_mixed/N_mixed_evts>2) | (N_daughter_mixed %/% N_mixed_evts==3 & N_daughter_mixed/N_mixed_evts==3)){
          for (j in (1:(N_daughter_mixed-2*N_mixed_evts))){
            Daughters_mixed[j]<-Daughters_mixed[j]+1
          }
        }
        if (N_daughter_mixed %/% N_mixed_evts==3 & N_daughter_mixed/N_mixed_evts>3){
          Daughters_mixed<-rep(3,N_mixed_evts)
          for (j in (1:(N_daughter_mixed-3*N_mixed_evts))){
            Daughters_mixed[j]<-Daughters_mixed[j]+1
          }
        }
        for (i in 1:N_mixed_evts){
          # Calculate merging area fraction (normal distribution)
          if (i==1){
            Init_cliffs<-Cliffs2mix[1:Mothers_mixed[1]]
          }
          if (i>1){
            Init_cliffs<-Cliffs2mix[(1+sum(Mothers_mixed[1:(i-1)])):sum(Mothers_mixed[1:i])]
          }
          Init_size<-0
          for (j in Init_cliffs){
            Init_size<-Init_size+Cliff_pop[[j]][1]
          }
          Fin_size<-0
          # Calculate area change 
          while ((Fin_size<param_list[8]*Daughters_mixed[i]) | (Fin_size>param_list[7]*Daughters_mixed[i])){
            if (param_list[35]*log(Init_size)+param_list[36]>=0){
              frac_mix<-rlnorm(1,param_list[33]*log(Init_size)+param_list[34],param_list[35]*log(Init_size)+param_list[36])
            }
            Fin_size<-Init_size*frac_mix
          }
          # Distribute final area to all daughters
          Size_daughters<-rep(0,Daughters_mixed[i])
          for (j in 1:Daughters_mixed[i]){
            if (j<Daughters_mixed[i]){
              Size_daughters[j]<-runif(1,param_list[8],Fin_size-sum(Size_daughters)-(Daughters_mixed[i]-j)*param_list[8])
            }
            if (j==Daughters_mixed[i]){
              Size_daughters[j]<-Fin_size-sum(Size_daughters)
            }
          }
          # Modify cliff population based on number of mothers vs number of daughters
          if (Daughters_mixed[i]<=Mothers_mixed[i]){
            # Change size of first cliffs
            for (j in 1:Daughters_mixed[i]){
              Cliff_pop[Init_cliffs[j]][1]<-round(Size_daughters[j])
            }
            # Kill all other
            if (Daughters_mixed[i]<Mothers_mixed[i]){
              idx_kill<-c(idx_kill,Init_cliffs[(Daughters_mixed[i]+1):length(Init_cliffs)])
            }
          }
          if (Daughters_mixed[i]>Mothers_mixed[i]){
            # Change size of first cliffs
            for (j in 1:Mothers_mixed[i]){
              Cliff_pop[Init_cliffs[j]][1]<-round(Size_daughters[j])
            }
            # Add new ones
            for (j in (Mothers_mixed[i]+1):Daughters_mixed[i]){
              Cliff_pop[[length(Cliff_pop)+1]]<-round(Size_daughters[j])
            }
          }
        }
      }
      
      # Kill cliffs from merge/mixed events
      if (length(idx_kill)>0){
        Cliff_pop<-Cliff_pop[-idx_kill]
      }
      
      # BIRTH: Introduce new cliffs (based on birth rate + stochasticity), with initial means size based on initial size distribution.
      N_new_cliffs<-round(rnorm(1,param_list[43],param_list[44]))  # Number of new cliffs
      for (i in 1:N_new_cliffs){
        Cliff_pop[[length(Cliff_pop)+1]]<-round(rlnorm(1,param_list[45],param_list[46]))
      }
      
      # REMOVE all cliffs smaller than min_size
      idx_rm<-c()
      for (i in 1:length(Cliff_pop)){
        if (Cliff_pop[[i]][1] < param_list[8]){
          idx_rm<-c(idx_rm,i)
        }
      }
      if (length(idx_rm)>0){
        Cliff_pop<-Cliff_pop[-idx_rm]
      }
      
      # Outputs of the loop
      All_cumul_nb<-c(All_cumul_nb,length(Cliff_pop))
      Cumul_area<-0
      for (i in 1:length(Cliff_pop)){
        Cumul_area<-Cumul_area+Cliff_pop[[i]][1]
      }
      All_cumul_area<-c(All_cumul_area,Cumul_area)
    }
    tot_sim_nb<-c(tot_sim_nb,All_cumul_nb)
    tot_sim_sz<-c(tot_sim_sz,All_cumul_area)
  }
  return(list(tot_sim_nb,tot_sim_sz))
}