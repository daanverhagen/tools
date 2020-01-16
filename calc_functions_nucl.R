calc_ss_gap <- function(meth_data, read_stats, s_state, first_last_flag, cl,no_cores) {
  ## calculates mm_gap or uu_gaps, num_a and num_m or num_u in each gap
  ## also calculates the begining and ending mu and um for each ss gap
  ## puts m/u at both ends of each read (using site_stats information)
  
  
  
  first_last_flag<-first_last_flag #decides whether to add groups of the state to count in the beginning and end of the read
  
  meth_data$read_name <- as.character(meth_data$read_name)                                                  #change read names from factors to characters
  meth_data$sample <- as.character(meth_data$sample)                                                        #change sample names from factors to characters
  meth_data <- as.data.frame(meth_data[,c('start', 'end', 'read_name', 'sample', 'hit_state', 'ref_seq')])  #save necessary methylation data
  read_stats <- as.data.frame(read_stats[,c('ref_seq', 'read_name','first_site', 'last_site','sample')])    #save necessary read statistics
  samples <- unique(meth_data$sample)                                                                       #select unique samples
  gap_DF_list <- list()                                                                                     #creates list to collect all data
  
  for(this_sample in samples) {
    
    read_stats_sample <- filter(read_stats, sample == this_sample)                      #table of read statistics for this sample
    meth_data_sample <- filter(meth_data, sample == this_sample)                        #table of methylation data for this sample
    meth_data_sample_read <- split(meth_data_sample, meth_data_sample$read_name)        #methylation data split by read name for this sample
    temp_read_data <- split(read_stats_sample, read_stats_sample$read_name)             #read statistics split by read name for this sample
    read_names <- unique(as.character(meth_data_sample$read_name))                      #list of read names in this sample
    read_chunks <- split(read_names,factor(c(1:no_cores)))                              #latter read names divided into chunks of about equal size for parallelisation
    #num_reads <- length(unique(as.character(meth_data_sample$read_name)))              #number of reads in this sample
    
    
    gap_df_list<-parLapply(cl,read_chunks,function(y){                                  #send one of the read_chunks to each core, collect in list
      data.table::rbindlist(lapply(y,function(i){                                       #iterate over reads in current read_chunk, collect in data::table
        if(nrow(meth_data_sample_read[[i]]) > 1) {
          
          if(first_last_flag){
            # put extra groups (with hit_state = s_state) at both ends of each read or just change the hit_state of the both ending groups to s_state
            
            if(min(meth_data_sample_read[[i]]$start) > temp_read_data[[i]]$first_site) {
              
              # add group at beginning
              
              temp_meth_df_start <- data.frame(start = temp_read_data[[i]]$first_site, end = temp_read_data[[i]]$first_site, read_name = i, sample = this_sample, hit_state = s_state, ref_seq = unique(meth_data_sample_read[[i]]))
              temp_meth_df_start$read_name <- as.character(temp_meth_df_start$read_name)
              temp_meth_df_start$sample <- as.character(temp_meth_df_start$sample)
              meth_data_sample_read[[i]] <- bind_rows(temp_meth_df_start, meth_data_sample_read[[i]])
            }
            
            else if(min(meth_data_sample_read[[i]]$start) == temp_read_data[[i]]$first_site) {
              
              #change first state to hit_state
              
              meth_data_sample_read[[i]]$hit_state[1] = s_state
            }
            
            
            if(max(meth_data_sample_read[[i]]$end) < temp_read_data[[i]]$last_site) {
              
              #add a group at last site
              
              temp_meth_df_end <- data.frame(start = temp_read_data[[i]]$last_site, end = temp_read_data[[i]]$last_site, read_name = i, sample = this_sample, hit_state = s_state, ref_seq = unique(meth_data_sample_read[[i]]))
              temp_meth_df_end$read_name <- as.character(temp_meth_df_end$read_name)
              temp_meth_df_end$sample <- as.character(temp_meth_df_end$sample)
              meth_data_sample_read[[i]] <- bind_rows(meth_data_sample_read[[i]], temp_meth_df_end)
            }
            
            else if(max(meth_data_sample_read[[i]]$end) == temp_read_data[[i]]$last_site) {
              
              #change last state to hit_state
              
              meth_data_sample_read[[i]]$hit_state[nrow(meth_data_sample_read[[i]])] = s_state
            }
          }
          
          
          #calculating all observables
          group_pos <- meth_data_sample_read[[i]]$hit_state == s_state # m/u groups are detected as (TRUE/FALSE)s
          group_pos_temp <- which(group_pos)
          
          if(length(group_pos_temp)>1){
            #if there is more than one group to count
            ss_start = meth_data_sample_read[[i]]$end[group_pos][-length(which(group_pos))]
            ss_end = meth_data_sample_read[[i]]$start[group_pos][-1]
            ss_gap <- ss_end-ss_start-1
            
            
            #count ambiguous and not_s_state for all ss-gaps
            num_ambiguous <- unlist(calc_counts(meth_data_sample_read[[i]][,'hit_state'], group_pos_temp, count_state = 0))
            if (s_state == +1) {
              num_m_or_u <- unlist(calc_counts(meth_data_sample_read[[i]][,'hit_state'], group_pos_temp, count_state = -1))
            }
            else if (s_state == -1) {
              num_m_or_u <- unlist(calc_counts(meth_data_sample_read[[i]][,'hit_state'], group_pos_temp, count_state = +1))
            }
            
            
            #calculate distance to following group
            sx_gap <- meth_data_sample_read[[i]]$start[group_pos_temp+1][-length(group_pos_temp)] - meth_data_sample_read[[i]]$end[group_pos_temp][-length(group_pos_temp)] -1
            sx_gap[num_ambiguous + num_m_or_u == 0] <- NA
            
            
            #calculate distance to previous group
            #sourced out the calculation to cpp due to weird behaviour for negative indices in R
            xs_gap <- create_xs_gap(meth_data_sample_read[[i]]$start,meth_data_sample_read[[i]]$end, group_pos_temp,s_state)
            xs_gap[num_ambiguous + num_m_or_u == 0] <- NA
          }
          
          else {
            #if there is only one group or none, put NA everywhere
            ss_gap<-NA;
            ss_start<-NA;
            ss_end<-NA;
            num_m_or_u<-NA;
            num_ambiguous<-NA;
            sx_gap<-NA;
            xs_gap<-NA;
          }
          
          if (s_state == +1) {
            gap_df_list_i <- data.frame(mm_gap = ss_gap, start = ss_start, end = ss_end, num_u = num_m_or_u, num_a = num_ambiguous, mu_gap = sx_gap, um_gap = xs_gap)
          }
          
          else if (s_state == -1) {
            gap_df_list_i <- data.frame(uu_gap = ss_gap, start = ss_start, end = ss_end, num_m = num_m_or_u, num_a = num_ambiguous, um_gap = sx_gap, mu_gap = xs_gap)
          }
          
          gap_df_list_i$read_name <- meth_data_sample_read[[i]]$read_name[[1]]
          gap_df_list_i$sample <- this_sample
          gap_df_list_i$ref_seq <- (meth_data_sample_read[[i]]$ref_seq)[1]
          
          
          
          gap_df_list_i
        }
        
      })  ) #end of inner lapply
      
    }) #end of outer parLapply
    
    gap_DF_list[[this_sample]] <- data.table::rbindlist(gap_df_list)
    print(this_sample)
  }
  gap_df <- data.table::rbindlist(gap_DF_list)
  
  
  gap_df$read_name <- as.factor(gap_df$read_name)
  gap_df$sample <- factor(gap_df$sample, levels=samples)
  
  
  gap_df
}


calc_ss_gap2 <- function(meth_data, read_stats, s_state) {
  ## calculates mm_gap or uu_gaps, num_a and num_m or num_u in each gap
  ## also calculates the begining and ending mu and um for each ss gap
  ## puts m/u at both ends of each read (using site_stats information)
  
  calc_counts <- function(m_data, group_pos_temp, count_state) {
    count_list <- list()
    for (i in 1:(length(group_pos_temp)-1)) {
      counter <- 0
      for (j in seq(group_pos_temp[i],group_pos_temp[i+1],1)) {
        if(m_data[j,'hit_state'] == count_state) {counter <- counter + 1}
      }
      count_list[i] = counter
    }
    
    count_list
  }
  
  meth_data$read_name <- as.character(meth_data$read_name)
  meth_data$sample <- as.character(meth_data$sample)
  samples <- unique(as.character(meth_data$sample))

  gap_DF_list <- list()
  
  for(this_sample in samples) {
    meth_data_sample <- subset(meth_data, sample == this_sample)
    meth_data_sample_read <- split(meth_data_sample, meth_data_sample$read_name)
    read_names <- unique(as.character(meth_data_sample$read_name))
    num_reads <- length(unique(as.character(meth_data_sample$read_name)))

    read_stats_sample <- subset(read_stats, sample == this_sample)
    
    gap_df_list <- list()
    
    for(i in read_names) {
      if(nrow(meth_data_sample_read[[i]]) > 1) {

        ## put extra groups (with hit_state = s_state) at both ends of each read or just change the hit_state of the both ending groups to s_state
        # temp_read_data <- read_stats_sample[read_stats_sample$read_name == i,]
        temp_read_data <- subset(read_stats_sample, read_name %in% i)
        meth_data_sample_read[[i]] <- meth_data_sample_read[[i]][,c('start', 'end', 'read_name', 'sample', 'hit_state', 'ref_seq')]

        if(min(meth_data_sample_read[[i]]$start) > temp_read_data$first_site) {
          temp_meth_df_start <- data.frame(start = temp_read_data$first_site, end = temp_read_data$first_site, read_name = i, sample = this_sample, hit_state = s_state, ref_seq = unique(meth_data_sample_read[[i]]))
          temp_meth_df_start$read_name <- as.character(temp_meth_df_start$read_name)
          temp_meth_df_start$sample <- as.character(temp_meth_df_start$sample)
          meth_data_sample_read[[i]] <- bind_rows(temp_meth_df_start, meth_data_sample_read[[i]])
        }
        else if(min(meth_data_sample_read[[i]]$start) == temp_read_data$first_site) {
          meth_data_sample_read[[i]]$hit_state[1] = s_state
        }
        if(max(meth_data_sample_read[[i]]$end) < temp_read_data$last_site) {
          temp_meth_df_end <- data.frame(start = temp_read_data$last_site, end = temp_read_data$last_site, read_name = i, sample = this_sample, hit_state = s_state, ref_seq = unique(meth_data_sample_read[[i]]))
          temp_meth_df_end$read_name <- as.character(temp_meth_df_end$read_name)
          temp_meth_df_end$sample <- as.character(temp_meth_df_end$sample)
          meth_data_sample_read[[i]] <- bind_rows(meth_data_sample_read[[i]], temp_meth_df_end)
        }
        else if(max(meth_data_sample_read[[i]]$end) == temp_read_data$last_site) {
          meth_data_sample_read[[i]]$hit_state[nrow(meth_data_sample_read[[i]])] = s_state
        }

        group_pos <- meth_data_sample_read[[i]]$hit_state == s_state # m/u groups are detected as (TRUE/FALSE)s
        group_pos_temp <- which(group_pos)

        ss_gap <- meth_data_sample_read[[i]]$start[group_pos][-1] - meth_data_sample_read[[i]]$end[group_pos][-length(which(group_pos))] - 1
        ss_start = meth_data_sample_read[[i]]$end[group_pos][-length(which(group_pos))]
        ss_end = meth_data_sample_read[[i]]$start[group_pos][-1]

        num_ambiguous <- unlist(calc_counts(meth_data_sample_read[[i]], group_pos_temp, count_state = 0))
        if (s_state == +1) {
          num_m_or_u <- unlist(calc_counts(meth_data_sample_read[[i]], group_pos_temp, count_state = -1))
        }
        else if (s_state == -1) {
          num_m_or_u <- unlist(calc_counts(meth_data_sample_read[[i]], group_pos_temp, count_state = +1))
        }
        
        sx_gap <- meth_data_sample_read[[i]]$start[group_pos_temp+1][-length(group_pos_temp)] - meth_data_sample_read[[i]]$end[group_pos_temp][-length(group_pos_temp)] -1
        sx_gap[num_ambiguous + num_m_or_u == 0] <- NA
        
        xs_gap <- meth_data_sample_read[[i]]$start[group_pos_temp][-1] - meth_data_sample_read[[i]]$end[group_pos_temp-1] -1
        xs_gap[num_ambiguous + num_m_or_u == 0] <- NA
        
        gap_df_list[[i]] <- data.frame(ss_gap = ss_gap, start = ss_start, end = ss_end, num_m_or_u = num_m_or_u, num_a = num_ambiguous, sx_gap = sx_gap, xs_gap = xs_gap)

        gap_df_list[[i]]$read_name <- i
        gap_df_list[[i]]$sample <- this_sample
        gap_df_list[[i]]$ref_seq <- unique(meth_data_sample_read[[i]]$ref_seq)
        if (s_state == +1) {
          names(gap_df_list[[i]])[names(gap_df_list[[i]]) == "num_m_or_u"] <- "num_u"
          names(gap_df_list[[i]])[names(gap_df_list[[i]]) == "ss_gap"] <- "mm_gap"
          names(gap_df_list[[i]])[names(gap_df_list[[i]]) == "sx_gap"] <- "mu_gap"
          names(gap_df_list[[i]])[names(gap_df_list[[i]]) == "xs_gap"] <- "um_gap"
        } 
        else if (s_state == -1) {
          names(gap_df_list[[i]])[names(gap_df_list[[i]]) == "num_m_or_u"] <- "num_m"
          names(gap_df_list[[i]])[names(gap_df_list[[i]]) == "ss_gap"] <- "uu_gap"
          names(gap_df_list[[i]])[names(gap_df_list[[i]]) == "sx_gap"] <- "um_gap"
          names(gap_df_list[[i]])[names(gap_df_list[[i]]) == "xs_gap"] <- "mu_gap"
        }
      }
    }
    gap_DF_list[[this_sample]] <- bind_rows(gap_df_list)
  }
  gap_df <- bind_rows(gap_DF_list)
  
  if(nrow(gap_df) >0 ) {
    gap_df$read_name <- as.factor(gap_df$read_name)
    gap_df$sample <- factor(gap_df$sample, levels=samples)
  }
  
  gap_df
}



calc_mm_gap <- function(meth_data) {
  ## calculates mm_gap, num_a and num_u in each gap
  ## puts m/u at both ends of each read (using only meth_data information)
  
  calc_zeros <- function(m_data, meth_pos_temp) {
    num_zeros_list <- list()

    for (i in 1:(length(meth_pos_temp)-1)) {
      counter_zeros <- 0
      
      for (j in seq(meth_pos_temp[i],meth_pos_temp[i+1],1)) {
        if(m_data[j,'hit_state'] == 0)
          counter_zeros <- counter_zeros + 1
      }
      num_zeros_list[i] = counter_zeros
    }
    
    num_zeros_list
  }
  
  calc_minus_ones <- function(m_data, meth_pos_temp) {
    num_minus_ones_list <- list()

    for (i in 1:(length(meth_pos_temp)-1)) {
      counter_minus_ones <- 0
      
      for (j in seq(meth_pos_temp[i],meth_pos_temp[i+1],1)) {
        if(m_data[j,'hit_state'] == -1)
          counter_minus_ones <- counter_minus_ones + 1
      }
      num_minus_ones_list[i] = counter_minus_ones
    }
    
    num_minus_ones_list
  }
  
  meth_data$read_name <- as.character(meth_data$read_name)
  samples <- unique(as.character(meth_data$sample))
  
  gap_DF_list <- list()
  
  for(this_sample in samples) {
    meth_data_sample <- subset(meth_data, sample == this_sample)
    meth_data_sample_read <- split(meth_data_sample, meth_data_sample$read_name)
    read_names <- unique(as.character(meth_data_sample$read_name))
    num_reads <- length(unique(as.character(meth_data_sample$read_name)))

    gap_df_list <- list()
    
    if(num_reads > 0) {
      for(i in read_names) {
        if(nrow(meth_data_sample_read[[i]]) > 1) {
          meth_pos <- meth_data_sample_read[[i]]$hit_state == 1 #methy groups (TRUE/FALSE)
          
          meth_pos_temp <- which(meth_pos)
          if(meth_pos[1] != TRUE)
            meth_pos_temp <- append(1, meth_pos_temp)
          if(meth_pos[length(meth_pos)] != TRUE)
            meth_pos_temp <- append(meth_pos_temp, length(meth_pos))
  
          if(length(which(meth_pos)) >= 1) {
            mm_gap <- meth_data_sample_read[[i]]$start[meth_pos][-1] - meth_data_sample_read[[i]]$end[meth_pos][-length(which(meth_pos))] - 1
            mm_start = meth_data_sample_read[[i]]$end[meth_pos][-length(which(meth_pos))]
            mm_end = meth_data_sample_read[[i]]$start[meth_pos][-1]
            
            if(meth_pos[1] != TRUE) {
              mm_gap_begining = meth_data_sample_read[[i]]$start[meth_pos][1] - min(meth_data_sample_read[[i]]$start) - 1
              mm_gap <- append(mm_gap_begining, mm_gap)
              mm_start <- append(min(meth_data_sample_read[[i]]$start), mm_start)
              mm_end <- append(meth_data_sample_read[[i]]$start[meth_pos][1], mm_end)
              }
            
            if(meth_pos[length(meth_pos)] != TRUE) {
              mm_gap_ending = max(meth_data_sample_read[[i]]$end) - meth_data_sample_read[[i]]$end[meth_pos][length(which(meth_pos))] -1
              mm_gap <- append(mm_gap, mm_gap_ending)
              mm_start <- append(mm_start, meth_data_sample_read[[i]]$end[meth_pos][length(which(meth_pos))])
              mm_end <- append(mm_end, max(meth_data_sample_read[[i]]$end))
            }
          }
          else if(length(which(meth_pos)) == 0) {
            mm_gap <- max(meth_data_sample_read[[i]]$end) - min(meth_data_sample_read[[i]]$start) - 1
            mm_start = min(meth_data_sample_read[[i]]$start)
            mm_end = max(meth_data_sample_read[[i]]$end)
          }
          
          num_zeros <- unlist(calc_zeros(meth_data_sample_read[[i]], meth_pos_temp))
          num_minus_ones <- unlist(calc_minus_ones(meth_data_sample_read[[i]], meth_pos_temp))
          
          gap_df_list[[i]] <- data.frame(mm_gap = mm_gap, mm_start = mm_start, mm_end = mm_end, num_minus_ones = num_minus_ones, num_zeros = num_zeros)
  
          gap_df_list[[i]]$read_name <- i
          gap_df_list[[i]]$sample <- this_sample
        }
      }
      gap_DF_list[[this_sample]] <- bind_rows(gap_df_list)
    }
  }
  gap_df <- bind_rows(gap_DF_list)
  
  if(nrow(gap_df) > 0 ) {
    gap_df$read_name <- as.factor(gap_df$read_name)
    gap_df$sample <- factor(gap_df$sample, levels=samples)
  }
  
  gap_df
}



calc_island_size <- function(meth_data, read_stats, s_state) {
  ## calculates the size (length) of m-only (or u-ony) islands surrounded by u sites (or m sites)
  
  meth_data$read_name <- as.character(meth_data$read_name)
  meth_data$sample <- as.character(meth_data$sample)
  samples <- unique(as.character(meth_data$sample))
  
  island_DF_list <- list()
  
  for(this_sample in samples) {
    meth_data_sample <- subset(meth_data, sample == this_sample)
    meth_data_sample_read <- split(meth_data_sample, meth_data_sample$read_name)
    read_names <- unique(as.character(meth_data_sample$read_name))
    num_reads <- length(unique(as.character(meth_data_sample$read_name)))
    
    read_stats_sample <- subset(read_stats, sample == this_sample)
    
    island_df_list <- list()
    
    for(i in read_names) {
      a_pos <- meth_data_sample_read[[i]]$hit_state == 0
      if(length(which(a_pos)) > 0)
        meth_data_sample_read[[i]][meth_data_sample_read[[i]]$hit_state == 0,]$hit_state = -1*s_state ## turn the state of all 'a' groups into 'u'
        # meth_data_sample_read[[i]][meth_data_sample_read[[i]]$hit_state == 0,]$hit_state = 1*s_state ## turn the state of all 'a' groups into 'm'

      ## put extra groups (with hit_state = s_state) at both ends of each read or just change the hit_state of the both ending groups to -1*s_state
      temp_read_data <- read_stats_sample[read_stats_sample$read_name == i,]
      meth_data_sample_read[[i]] <- meth_data_sample_read[[i]][,c('start', 'end', 'read_name', 'sample', 'hit_state', 'ref_seq')]

      if(temp_read_data$first_site < min(meth_data_sample_read[[i]]$start) && meth_data_sample_read[[i]]$hit_state[1] == s_state) {
        temp_meth_df_start <- data.frame(start = temp_read_data$first_site, end = temp_read_data$first_site, read_name = i, sample = this_sample, hit_state = -1*s_state, ref_seq = unique(meth_data_sample_read[[i]]))
        temp_meth_df_start$read_name <- as.character(temp_meth_df_start$read_name)
        temp_meth_df_start$sample <- as.character(temp_meth_df_start$sample)
        meth_data_sample_read[[i]] <- bind_rows(temp_meth_df_start, meth_data_sample_read[[i]])
      }
      else if(temp_read_data$first_site == min(meth_data_sample_read[[i]]$start) && meth_data_sample_read[[i]]$hit_state[1] == s_state) {
        meth_data_sample_read[[i]]$hit_state[1] = -1*s_state
      }
      if(max(meth_data_sample_read[[i]]$end) < temp_read_data$last_site && meth_data_sample_read[[i]]$hit_state[nrow(meth_data_sample_read[[i]])] == s_state) {
        temp_meth_df_end <- data.frame(start = temp_read_data$last_site, end = temp_read_data$last_site, read_name = i, sample = this_sample, hit_state = -1*s_state, ref_seq = unique(meth_data_sample_read[[i]]))
        temp_meth_df_end$read_name <- as.character(temp_meth_df_end$read_name)
        temp_meth_df_end$sample <- as.character(temp_meth_df_end$sample)
        meth_data_sample_read[[i]] <- bind_rows(meth_data_sample_read[[i]], temp_meth_df_end)
      }
      else if(max(meth_data_sample_read[[i]]$end) == temp_read_data$last_site && meth_data_sample_read[[i]]$hit_state[nrow(meth_data_sample_read[[i]])] == s_state) {
        meth_data_sample_read[[i]]$hit_state[nrow(meth_data_sample_read[[i]])] = -1*s_state
      }
      
      x_pos <- meth_data_sample_read[[i]]$hit_state == -1*s_state # u groups are detected as (TRUE/FALSE)s
      s_pos <- meth_data_sample_read[[i]]$hit_state == s_state # u groups are detected as (TRUE/FALSE)s

      if(length(which(x_pos)) >= 2 & length(which(s_pos)) > 1) {
        island_start = meth_data_sample_read[[i]]$start[which(x_pos)+1][-length(which(x_pos))]
        island_end = meth_data_sample_read[[i]]$end[which(x_pos)-1]
        island_size <- island_end - island_start

        island_df_list[[i]] <- data.frame(island_size = island_size, start = island_start, end = island_end)
        island_df_list[[i]] <- island_df_list[[i]][island_df_list[[i]]$island_size >= 0,]
        island_df_list[[i]]$read_name <- unique(as.character(i))
        island_df_list[[i]]$sample <- unique(as.character(this_sample))
        island_df_list[[i]]$ref_seq <- unique(meth_data_sample_read[[i]]$ref_seq)
        if(s_state == +1)
          names(island_df_list[[i]])[names(island_df_list[[i]]) == "island_size"] <- "m_island_size"
        else if(s_state == -1)
          names(island_df_list[[i]])[names(island_df_list[[i]]) == "island_size"] <- "u_island_size"
      }
    }
    island_DF_list[[this_sample]] <- bind_rows(island_df_list)
  }
  island_df <- bind_rows(island_DF_list)
  
  if(nrow(island_df) > 0 ) {
    island_df$read_name <- as.factor(island_df$read_name)
    island_df$sample <- factor(island_df$sample, levels=samples)
  }
  
  island_df
}



calc_m_island <- function(umu_gap) {
  ## calculates the size (length) of m-only (or u-ony) islands surrounded by u sites (or m sites)
  
  umu_gap$read_name <- as.character(umu_gap$read_name)
  umu_gap$sample <- as.character(umu_gap$sample)
  samples <- unique(as.character(umu_gap$sample))
  
  island_DF_list <- list()
  
  for(this_sample in samples) {
    umu_gap_sample <- subset(umu_gap, sample == this_sample)
    umu_gap_sample_read <- split(umu_gap_sample, umu_gap_sample$read_name)
    read_names <- unique(as.character(umu_gap_sample$read_name))

    island_df_list <- list()
    
    for(i in read_names) {
      m_island_size <- umu_gap_sample_read[[i]]$uu_gap - umu_gap_sample_read[[i]]$um_gap - umu_gap_sample_read[[i]]$mu_gap
      m_island_start <- umu_gap_sample_read[[i]]$start + umu_gap_sample_read[[i]]$mu_gap
      m_island_end <- umu_gap_sample_read[[i]]$end + umu_gap_sample_read[[i]]$um_gap
      
      island_df_list[[i]] <- data.frame(m_island_size = m_island_size, uu_gap = umu_gap_sample_read[[i]]$uu_gap, start = umu_gap_sample_read[[i]]$start, end = umu_gap_sample_read[[i]]$end, um_gap = umu_gap_sample_read[[i]]$mu_gap, mu_gap = umu_gap_sample_read[[i]]$um_gap)
      island_df_list[[i]]$read_name <- unique(as.character(i))
      island_df_list[[i]]$sample <- unique(as.character(this_sample))
      island_df_list[[i]]$ref_seq <- unique(umu_gap_sample_read[[i]]$ref_seq)
    }
    island_DF_list[[this_sample]] <- bind_rows(island_df_list)
  }
  island_df <- bind_rows(island_DF_list)
  
  if(nrow(island_df) > 0 ) {
    island_df$read_name <- as.factor(island_df$read_name)
    island_df$sample <- factor(island_df$sample, levels=samples)
  }
  
  island_df
}


calc_gap_summary_by_sample <- function(gap_df) {
  
  calc_maxima_under_500 <- function(gaps) {
    paste0(signif(calc_local_density_maxima(gaps, xlim = c(0, 500)), digits = 3), collapse = ", ")
  }
  signif_mean <- function(x) {
    signif(mean(x), digits = 3)
  }
  
  gap_summary_by_sample <- lapply(list(signif_mean), function(f) {aggregate(gap_df[,1], by=list(gap_df$sample), f)[, 2]})
  
  gap_summary_by_sample <- lapply(list(signif_mean, quantile, calc_maxima_under_500, length), function(f) {aggregate(gap_df[,1], by=list(gap_df$sample), f)[, 2]})
  gap_summary_by_sample <- do.call(cbind, gap_summary_by_sample)
  rownames(gap_summary_by_sample) <- unique(gap_df$sample)
  colnames(gap_summary_by_sample)[1:8] <- c("mean", "min", "q25", "median", "q75", "max", "gaps_with_max_dens", "num_of_gaps")

  return(as.data.frame(gap_summary_by_sample))
}



calc_gap_summary_by <- function(gap_df) {

  calc_maxima_under_500 <- function(gaps) {
    paste0(signif(calc_local_density_maxima(gaps, xlim = c(0, 500)), digits = 3), collapse = ", ")
  }
  signif_mean <- function(x) {
    signif(mean(x), digits = 3)
  }
  
  gap_summary_by_col_name <- lapply(list(signif_mean, quantile, calc_maxima_under_500), function(f) {aggregate(gap_df[,1], by=list(gap_df[,4]), f)[, 2]})
  gap_summary_by_col_name <- do.call(cbind, gap_summary_by_col_name)
  rownames(gap_summary_by_col_name) <- unique(gap_df[,4])
  colnames(gap_summary_by_col_name)[1:7] <- c("mean", "min", "q25", "median", "q75", "max", "gaps_with_max_dens")
  
  return(as.data.frame(gap_summary_by_col_name))
}



calc_large_gap_fraction_by_sample <- function(gap_df, length_threshold) {
  
  sum_all_gaps_df = aggregate(gap_df[,1],by=list(gap_df$sample), sum)
  names(sum_all_gaps_df) <- c('sample', 'sum_all_gaps')
  
  gap_df_large <- gap_df[gap_df$gap_length >= length_threshold,]
  sum_large_gaps_df = aggregate(gap_df_large[,1],by=list(gap_df_large$sample), sum)
  names(sum_large_gaps_df) <- c('sample', 'sum_large_gaps')
  
  
  gap_sum_by_sample = merge(sum_all_gaps_df, sum_large_gaps_df, by="sample")
  gap_sum_by_sample$frac_large_gaps <- gap_sum_by_sample$sum_large_gaps/gap_sum_by_sample$sum_all_gaps
  
  return(as.data.frame(gap_sum_by_sample))
}



calc_all_gaps <- function(site_stats) {
  ## calculates the distance between all successive groups

  ## order site stats according to start  
  site_stats_sample <- split(site_stats, site_stats$sample)
  site_stats_sample_ordered <- list()
  for(this_sample in unique(site_stats$sample)) {
    site_stats_ref_seq <- list()
    site_stats_ref_seq_ordered <- list()
    for(this_ref_seq in unique(site_stats_sample[[this_sample]]$ref_seq)){
      site_stats_ref_seq[[this_ref_seq]] <- subset(site_stats_sample[[this_sample]], ref_seq %in% this_ref_seq)
      site_stats_ref_seq_ordered[[this_ref_seq]] <- site_stats_ref_seq[[this_ref_seq]][order(site_stats_ref_seq[[this_ref_seq]]$start),]
    }
    site_stats_sample_ordered[[this_sample]] <- bind_rows(site_stats_ref_seq_ordered)
  }
  
  site_stats <- bind_rows(site_stats_sample_ordered)
  
  ref_seqs <- unique(as.character(site_stats$ref_seq))
  
  gap_df_list <- list()
  
  for(this_ref_seq in ref_seqs) {
    site_stats_ref_seq <- subset(site_stats, ref_seq == this_ref_seq)
    
    if (nrow(site_stats_ref_seq) >= 2) {
      gap_length <- site_stats_ref_seq$start[-1] - site_stats_ref_seq$end[-length(site_stats_ref_seq$end)] - 1
      gap_start = site_stats_ref_seq$end[-length(site_stats_ref_seq$end)]
      gap_end = site_stats_ref_seq$start[-1]
      
      gap_df_list[[this_ref_seq]] <- data.frame(gap_length = gap_length, gap_start = gap_start, gap_end = gap_end)
      gap_df_list[[this_ref_seq]]$ref_seq = this_ref_seq
    }
  }
  all_gaps_df <- bind_rows(gap_df_list)
  
  if(nrow(all_gaps_df) > 0 ) {
    all_gaps_df$ref_seq <- as.factor(all_gaps_df$ref_seq)
  }
  
  all_gaps_df
}



calc_dyad_position <- function(mum_gap_df, nucl_footprint, max_overlap) {
  ## places appropriate number of nucleosomes in m-u-m gaps
  
  mum_gap_df <- subset(mum_gap_df, num_u >= 1)
  
  mum_gap_df$read_name <- as.character(mum_gap_df$read_name)
  samples <- unique(as.character(mum_gap_df$sample))
  
  nucl_DF_list <- list()

  for(this_sample in samples) {
    mum_gap_df_sample <- subset(mum_gap_df, sample == this_sample)
    mum_gap_df_sample_read <- split(mum_gap_df_sample, mum_gap_df_sample$read_name)
    read_names <- unique(as.character(mum_gap_df_sample$read_name))
    num_reads <- length(unique(as.character(mum_gap_df_sample$read_name)))

    nucl_df_list <- list()
    for(i in read_names) {
      num_nucl <- floor((mum_gap_df_sample_read[[i]]$mm_gap + max_overlap)/(nucl_footprint - max_overlap))
      nucl_gap <- num_nucl > 0 # (TRUE/FALSE) #indicates mm_gaps that can accomodate at least one nucl
      dyad_df_list <- list()
      if(length(which(nucl_gap)) > 0) {
        # print(i)
        for (k in 1:length(which(nucl_gap))) {
          overlap <- (mum_gap_df_sample_read[[i]]$mm_gap[nucl_gap][k] - num_nucl[nucl_gap][k]*nucl_footprint)/(num_nucl[nucl_gap][k]+1)
          position <- mum_gap_df_sample_read[[i]]$start[nucl_gap][k] + overlap + 0.5*nucl_footprint + seq(0, (num_nucl[nucl_gap][k]-1), 1)*(overlap + nucl_footprint)
          # position <- mum_gap_df_sample_read[[i]]$start[nucl_gap][k] + mum_gap_df_sample_read[[i]]$mm_gap[nucl_gap][k]*seq(1, 2*num_nucl[nucl_gap][k], 2)/(2*num_nucl[nucl_gap][k])
          dyad_df_list[[k]] <- data.frame(position = floor(position))
          dyad_df_list[[k]]$read_name <- i
          dyad_df_list[[k]]$sample = this_sample
          dyad_df_list[[k]]$ref_seq = unique(as.character(mum_gap_df_sample_read[[i]]$ref_seq))
        }
        nucl_df_list[[i]] <- bind_rows(dyad_df_list)
      }
    }
    nucl_DF_list[[this_sample]] <- bind_rows(nucl_df_list)
  }
  nucl_df <- bind_rows(nucl_DF_list)
  
  if(nrow(nucl_df) > 0) {
    nucl_df$read_name <- as.factor(nucl_df$read_name)
    nucl_df$sample <- as.factor(nucl_df$sample)
    nucl_df$ref_seq <- as.factor(nucl_df$ref_seq)
  }
  
  nucl_df
}



calc_meth_size <- function(meth_data) {
  
  meth_data$read_name <- as.character(meth_data$read_name)
  samples <- unique(as.character(meth_data$sample))
  
  meth_size_DF_list <- list()
  
  for(this_sample in samples) {
    meth_data_sample <- subset(meth_data, sample == this_sample)
    meth_data_sample_read <- split(meth_data_sample, meth_data_sample$read_name)
    read_names <- unique(as.character(meth_data_sample$read_name))
    num_reads <- length(unique(as.character(meth_data_sample$read_name)))
    
    meth_size_df_list <- list()
    
    if(num_reads > 0) {
      for(i in read_names) {
        meth_pos <- meth_data_sample_read[[i]]$hit_state == 1 #methy groups (TRUE/FALSE)
        meth_size <- meth_data_sample_read[[i]]$end[meth_pos] - meth_data_sample_read[[i]]$start[meth_pos]
        m_start = meth_data_sample_read[[i]]$start[meth_pos]
        m_end = meth_data_sample_read[[i]]$end[meth_pos]
        meth_size_df_list[[i]] <- data.frame(read_name = i, sample = this_sample, meth_size = meth_size, m_start = m_start, m_end = m_end)
      }
      meth_size_DF_list[[this_sample]] <- bind_rows(meth_size_df_list)
    }
  }
  meth_size_df <- bind_rows(meth_size_DF_list)

  if(nrow(meth_size_df) > 0) {
    meth_size_df$read_name <- as.factor(meth_size_df$read_name)
    meth_size_df$sample <- factor(meth_size_df$sample, levels=samples)
  }
  
  meth_size_df
}



calc_nucl_spacing <- function(nucl_df)  {
  
  nucl_df$read_name <- as.character(nucl_df$read_name)
  nucl_df$sample <- as.character(nucl_df$sample)
  samples <- unique(as.character(nucl_df$sample))
  
  nucl_DF_read_list <- list()
  
  for (this_sample in samples) {
    nucl_df_sample <- subset(nucl_df, sample == this_sample)
    nucl_df_read <- split(nucl_df_sample, nucl_df_sample$read_name)
    read_names <- unique(as.character(nucl_df_sample$read_name))
    num_reads <- length(unique(as.character(nucl_df_sample$read_name)))
    
    nucl_df_read_list <- list()
    
    for(i in read_names) {
      if (length(nucl_df_read[[i]]$position) >= 2) {
        dyad_dyad_distance <- diff(nucl_df_read[[i]]$position) - 1
        
        nucl_df_read_list[[i]] <- data.frame(dyad_dyad_distance = dyad_dyad_distance)
        nucl_df_read_list[[i]]$read_name = as.character(i)
        nucl_df_read_list[[i]]$sample = as.character(this_sample)
      }
    }
    nucl_DF_read_list[[this_sample]] <- bind_rows(nucl_df_read_list)
  }
  nucl_spacing_df <- bind_rows(nucl_DF_read_list)

  if(nrow(nucl_spacing_df) > 0) {
    nucl_spacing_df$read_name <- as.factor(nucl_spacing_df$read_name)
    nucl_spacing_df$sample <- factor(nucl_spacing_df$sample, levels=samples)
  }
  
  nucl_spacing_df
}



llr2external_pot <- function(meth_data, read_stats, nucl_footprint, region_limit, epsilon) {
  
  read_length <- read_stats$last_site - read_stats$first_site +1
  extern_pot <- rep(0, (read_length + 2*region_limit))

  meth_groups <- meth_data$hit_state == 1 #methy groups (TRUE/FALSE)
  meth_groups_temp <- which(meth_groups) #methy groups

  if(length(meth_groups_temp) >= 1) {
    
    meth_data$start <- meth_data$start + region_limit - read_stats$first_site +1
    meth_data$end <- meth_data$end + region_limit - read_stats$first_site +1

    for (k in meth_groups_temp) {
      ## Isosceles trapezoid potential (flat at the middle and triagular at both sides)
      x <- c((meth_data$start[k] - floor(nucl_footprint/2)):(meth_data$end[k] + floor(nucl_footprint/2)))
      n1 <- c(1:floor(nucl_footprint/2))
      n2 <- c((floor(nucl_footprint/2)+1):(length(x) - floor(nucl_footprint/2)))
      n3 <- c((length(x) - floor(nucl_footprint/2) + 1):length(x))

      extern_pot[x[n1]] <- extern_pot[x[n1]] + sapply(x[n1], function(x_pos) (x_pos - x[1]) * epsilon)
      extern_pot[x[n2]] <- extern_pot[x[n2]] + floor(nucl_footprint/2) * epsilon
      extern_pot[x[n3]] <- extern_pot[x[n3]] + sapply(x[n3], function(x_pos) (x[length(x)] - x_pos) * epsilon)
    }
    meth_data$start <- meth_data$start - region_limit + read_stats$first_site -1
    meth_data$end <- meth_data$end - region_limit + read_stats$first_site -1
    
    ## set the maximum height of the potetial equal to floor(nucl_footprint/2)*epsilon
    extern_pot[extern_pot > floor(nucl_footprint/2)*epsilon] <- floor(nucl_footprint/2)*epsilon
  }
  
  extern_pot
}



calc_phi_interact_pot <- function(epsilon, nucl_footprint) {
  ## nucl-nucl interaction potential ##
  ### measured by summing up all states corresponding to certain binding energy ###
  
  p <- nucl_footprint - 1
  w <- floor(p/(2))
  v <- rep(0,p)
  
  ## Apprximation
  # for (i in 1:p) {
  #     v[i] <- (nucl_footprint-i) * epsilon - log(1 + (nucl_footprint-i) * (1-exp(-epsilon)) )
  # }
  
  ## Exact calculation
  Eexp = matrix(0, ncol = p, nrow = p)
  for(i in 1:(p-1))
    for(j in 1:(p-1))
      Eexp[i,j] <- exp((i+j)*epsilon)
  
  y <- rep(0, p)
  for(n in 1:(w+1))
    for(i in 0:(n-1))
      for(j in 0:(n-1-i))
        y[n] <- y[n] + Eexp[i+1,j+1]
  
  for(n in (w+2):p)
    for(i in 0:w) {
      r = min(w, (n-1-i))
      for(j in 0:r)
        y[n] <- y[n] + Eexp[i+1,j+1]
    }
  
  Omega=0
  for(i in 0:w)
    for(j in 0:w)
      Omega <- Omega + Eexp[i+1,j+1];
  
  for(x in 1:p)
    v[x] <- -log(y[x]) + log(Omega)
  
  v
}



calc_neutralizing_potential <-function(epsilon, nucl_footprint, density_goal) {
  ## calculates an external potential to neutralize/compensate the boundary ##
  
  N = 4000
  a = nucl_footprint
  w = (a-1)/2
  density_goal = 0.006 # in Johannes' paper the density corresponding to the mu=4.012 is 0.006 in the middle half of the system (run with boundary particles)
  region_limit = floor(1.5*a)
  region = c(1:(region_limit))
  R = 6;  # 2^R runs
  dens_diff_to_pot = a

  u_temp = -4.012 * rep(1,N)  # initial external potential
  # u_temp = rep(0,N)  # initial external potential
  phi = c(calc_phi_interact_pot(epsilon, nucl_footprint),rep(0 ,N-2*w))

  u_list <- list()
  u_list[[1]] <- u_temp
  n_list <- list()
  for(i in 1:(2^R)) {
    print(paste0('iteration = ', i))
    u_temp_n_df <- NULL
    u_temp_n_df = minimize_density_fluctuation(u_list[[i]], phi, region, density_goal, dens_diff_to_pot)
    u_list[[i+1]] = u_temp_n_df[,1]
    n_list[[i]] = u_temp_n_df[,2]
  }
  u_list[[2^R+1]] <- NULL

  u_n_list <- list(u_list, n_list)

  u_n_list
}



minimize_density_fluctuation <- function(u, phi, region, density_goal, dens_diff_to_pot) {
  
  N = length(u)
  
  n = pots2dens_only_n(u, phi)
  
  u_new <- 0
  u_new[region] = u[region] + dens_diff_to_pot*(n[region] - density_goal)
  u_new[(max(region)+1):(N-max(region))] = u_new[max(region)]
  u_new[(N-max(region)+1):N] = u[(N-max(region)+1):N] + dens_diff_to_pot*(n[(N-max(region)+1):N] - density_goal)

  u_temp_n_df <- NULL
  u_temp_n_df <- data.frame('u_new' = u_new, 'n' = n)

  u_temp_n_df
}



pots2dens_only_n <-function(external_pot, phi) {
  ## Used notation as in ref [Percus, J. Phys. Condens. Matter 1989] ##
  ## from external (u) and interparticle (phi) potentials -> n1 & n2 ##
  
  read_length <- length(external_pot)
  
  z_matrix <- matrix(0, ncol = read_length, nrow = read_length)
  z_matrix <- diag(exp(-external_pot))
  z_vector <- diag(z_matrix)
  
  w_vec <- exp(-phi)
  w_matrix <- matrix(0, ncol = read_length, nrow = read_length)
  
  for(l in 1:(read_length-1)) {
    m = c((l+1):read_length)
    w_matrix[l,m] <- w_vec[m-l]
  }
  z_w_matrix <- matrix(0, ncol = read_length, nrow = read_length)
  for(l in 1:read_length) { ## faster than normal multiplication (verified)
    z_w_matrix[l,] <- w_matrix[l,]*z_vector[l]
  }
  w_z_matrix <- matrix(0, ncol = read_length, nrow = read_length)
  for(l in 1:read_length) {
    w_z_matrix[,l] <- w_matrix[,l]*z_vector[l]
  }
  
  ## Grand Canonical partition function
  izw_matrix <- diag(read_length) - z_w_matrix
  izw_inverse_matrix <- backsolve(izw_matrix, diag(dim(izw_matrix)[1]))
  izw_inverse_z_matrix <- matrix(0, ncol = read_length, nrow = read_length)
  for(l in 1:read_length) {
    izw_inverse_z_matrix[,l] <- izw_inverse_matrix[,l]*z_vector[l]
  }
  grand_partition_function <- 1 + sum(colSums(izw_inverse_z_matrix))

  ## One-particle density
  iwz_matrix <- diag(read_length) - w_z_matrix
  iwz_inverse_matrix <- backsolve(iwz_matrix, diag(dim(iwz_matrix)[1]))
  colSums_izw_inverse_matrix <- colSums(izw_inverse_matrix)
  rowSums_iwz_inverse_matrix <- rowSums(iwz_inverse_matrix)
  one_particle_density <- rep(0, read_length)
  cS_izw_inv_z_temp <- colSums_izw_inverse_matrix * z_vector
  z_rS_iwz_inv_temp <- z_vector * rowSums_iwz_inverse_matrix
  one_particle_density <- (cS_izw_inv_z_temp * rowSums_iwz_inverse_matrix)/grand_partition_function

  one_particle_density
}



pots2dens <-function(external_pot, phi) {
  ## Used notation as in ref [Percus, J. Phys. Condens. Matter 1989] ##
  ## from external (u) and interparticle (phi) potentials -> n1 & n2 ##
  
  read_length <- length(external_pot)
  
  z_matrix <- matrix(0, ncol = read_length, nrow = read_length)
  z_matrix <- diag(exp(-external_pot))
  z_vector <- diag(z_matrix)
  
  w_vec <- exp(-phi)
  w_matrix <- matrix(0, ncol = read_length, nrow = read_length)
  
  for(l in 1:(read_length-1)) {
    m = c((l+1):read_length)
    w_matrix[l,m] <- w_vec[m-l]
  }
  z_w_matrix <- matrix(0, ncol = read_length, nrow = read_length)
  for(l in 1:read_length) { ## faster than normal multiplication (verified)
    z_w_matrix[l,] <- w_matrix[l,]*z_vector[l]
  }
  w_z_matrix <- matrix(0, ncol = read_length, nrow = read_length)
  for(l in 1:read_length) {
    w_z_matrix[,l] <- w_matrix[,l]*z_vector[l]
  }
  
  ## Grand Canonical partition function
  izw_matrix <- diag(read_length) - z_w_matrix
  izw_inverse_matrix <- backsolve(izw_matrix, diag(dim(izw_matrix)[1]))
  izw_inverse_z_matrix <- matrix(0, ncol = read_length, nrow = read_length)
  for(l in 1:read_length) {
    izw_inverse_z_matrix[,l] <- izw_inverse_matrix[,l]*z_vector[l]
  }
  grand_partition_function <- 1 + sum(colSums(izw_inverse_z_matrix))

  ## One-particle density
  iwz_matrix <- diag(read_length) - w_z_matrix
  iwz_inverse_matrix <- backsolve(iwz_matrix, diag(dim(iwz_matrix)[1]))
  colSums_izw_inverse_matrix <- colSums(izw_inverse_matrix)
  rowSums_iwz_inverse_matrix <- rowSums(iwz_inverse_matrix)
  one_particle_density <- rep(0, read_length)
  cS_izw_inv_z_temp <- colSums_izw_inverse_matrix * z_vector
  z_rS_iwz_inv_temp <- z_vector * rowSums_iwz_inverse_matrix
  one_particle_density <- (cS_izw_inv_z_temp * rowSums_iwz_inverse_matrix)/grand_partition_function

  ## Two-particle density
  two_particle_density <- matrix(0, nrow = read_length, ncol = read_length)
  for(l in 1:read_length) {
    for (m in 1:read_length) {
      two_particle_density[l,m] <- cS_izw_inv_z_temp[l] * w_matrix[l,m] * z_rS_iwz_inv_temp[m]
    }
  }
  
  one_particle_density
}



calc_grand_canon <- function(meth_data, read_stats, neut_extern_pot, interact_pot, nucl_footprint, epsilon, cut_flag) {
  
  meth_data$read_name <- as.character(meth_data$read_name)
  samples <- as.character(unique(factor(meth_data$sample)))

  region_limit = floor(1.5*nucl_footprint)
  
  one_particle_density_DF_list <- list()
  meth_data_new_DF_list <-list()
  
  for(this_sample in samples) {
    print(this_sample)
    meth_data_sample <- subset(meth_data, sample == this_sample)
    meth_data_sample_read <- split(meth_data_sample, meth_data_sample$read_name)
    read_names <- unique(as.character(meth_data_sample$read_name))
    num_reads <- length(unique(as.character(meth_data_sample$read_name)))

    read_stats_sample <- subset(read_stats, sample == this_sample)
    
    one_particle_density_df_list <- list()
    
    for(i in read_names) {
      temp_read_data <- read_stats_sample[read_stats_sample$read_name == i,]
      read_length <- temp_read_data$last_site - temp_read_data$first_site + 1

      external_pot1 <- llr2external_pot(meth_data_sample_read[[i]], temp_read_data, nucl_footprint, region_limit, epsilon)
      external_pot2 <- c(neut_extern_pot, rep(neut_extern_pot[region_limit], read_length), neut_extern_pot[region_limit:1])
      external_pot <- external_pot1 + external_pot2
      
      phi <- c(interact_pot, rep(0, (length(external_pot)-(nucl_footprint-1))))
      
      one_particle_density <- pots2dens_only_n(external_pot, phi)
      
      if(cut_flag) {
        one_particle_density <- one_particle_density[(region_limit+1):(length(one_particle_density)-region_limit)]
        # external_pot <- external_pot[(region_limit+1):(length(external_pot)-region_limit)]
        # one_particle_density_df_list[[i]] <- data.frame(position = c(temp_read_data$first_site:temp_read_data$last_site), one_particle_density = one_particle_density, externl_pot = external_pot)
        one_particle_density_df_list[[i]] <- data.frame(position = c(temp_read_data$first_site:temp_read_data$last_site), one_particle_density = one_particle_density)
      } else {
        # one_particle_density_df_list[[i]] <- data.frame(position = c((temp_read_data$first_site-region_limit):(temp_read_data$last_site+region_limit)), one_particle_density = one_particle_density, externl_pot = external_pot)
        one_particle_density_df_list[[i]] <- data.frame(position = c((temp_read_data$first_site-region_limit):(temp_read_data$last_site+region_limit)), one_particle_density = one_particle_density)
      }
      
      one_particle_density_df_list[[i]]$sample = this_sample
    }
    one_particle_density_DF_list[[this_sample]] <- ldply(one_particle_density_df_list, data.frame)
    colnames(one_particle_density_DF_list[[this_sample]])[1] <- 'read_name'
  } ## for loop over samples
  one_particle_density_df <- bind_rows(one_particle_density_DF_list)
  
  if(nrow(one_particle_density_df) > 0) {
    one_particle_density_df$read_name <- as.factor(one_particle_density_df$read_name)
    one_particle_density_df$sample <- factor(one_particle_density_df$sample, levels=samples)
  }
  
  one_particle_density_df
}



calc_average_gene_1_particle_den <- function(gcanon_den, min_site, max_site) {

  gcanon_den$read_name <- as.character(gcanon_den$read_name)
  warning("Average over all reads (regardless of the sample/pooled_label)\n")
  
  ## Check if the beginning or end of the reads needs to be adjusted (aligned)
  gcanon_den_read <- split(gcanon_den, gcanon_den$read_name)
  read_names <- unique(as.character(gcanon_den$read_name))
  num_reads <- length(read_names)
  
  min_vec <- NULL
  max_vec <- NULL
  for(i in 1:num_reads) {
    min_vec[i] <- min(gcanon_den_read[[i]]$position)
    max_vec[i] <- max(gcanon_den_read[[i]]$position)
  }
  
  if(unique(min_vec) > 1 || unique(max_vec) > 1) {
    ## needs alignement
    density_list <- list()
    for(i in read_names) {
      gcanon_den_read[[i]] <- gcanon_den_read[[i]][order(gcanon_den_read[[i]]$position),]  ## sort position column otherwise this function mixes everything
      gcanon_den_read_temp <- data.frame(read_name = i, position = min_site:max_site, one_particle_density = NaN)

      pos_min = min(gcanon_den_read[[i]]$position)
      pos_max = max(gcanon_den_read[[i]]$position)

      one_particle_den_vec <- rep(NaN, (max_site-min_site+1))
      one_particle_den_vec[(pos_min-min_site+1):(pos_max-min_site+1)] <- gcanon_den_read[[i]]$one_particle_density
      gcanon_den_read_temp$one_particle_density <- one_particle_den_vec

      gcanon_den_read_temp$read_name <- unique(as.character(gcanon_den_read[[i]]$read_name))

      density_list[[i]] <- gcanon_den_read_temp
    }
    gcanon_den_read_new <- bind_rows(density_list)
  }
  else if (unique(min_vec) == 1 || unique(max_vec) == 1) {
    ## No need to alignement
    gcanon_den_read_new = gcanon_den
  }

  average_gene = aggregate(gcanon_den_read_new[, c('one_particle_density')], by = list(gcanon_den_read_new$position), FUN = mean, na.rm = TRUE)
  colnames(average_gene) <- c('position', 'one_particle_density')

  average_gene
}



calc_inter_particle_distance <- function(two_particle_density) {
  
  inter_particle_distance_list <- list()
  read_length <- nrow(two_particle_density)
  for (k in 1:(read_length-1)) {
    inter_particle_distance_list[[k]] <- sapply(c(1:(read_length-k)), function(i) two_particle_density[i,i+k])
  }
  
  mean_inter_particle_distance <- sapply(c(1:(read_length-1)), function(k) mean(inter_particle_distance_list[[k]]))

  mean_inter_particle_distance
}



align_grand_canon_Np1 <- function(density_data, chr_ref_seqs, max_promoter_length, max_ORF_length) {
  
  density_data_ref_seq_list <- list()  # list of density_data with ref_seq as index
  for (ref_seq_label in chr_ref_seqs) {
    density_data_ref_seq_list = c(density_data_ref_seq_list, list(density_data[density_data$ref_seq == ref_seq_label, ]))
  }
  names(density_data_ref_seq_list) <- chr_ref_seqs
  
  gene_data = read.table(file = '../../../external_data/ann_plus1_minus1_filtered.txt', sep = '\t', stringsAsFactors = FALSE, header = TRUE)
  colnames(gene_data)[2] <- c("chr")
  chr_names = data.frame(chr_ref_seqs, as.character(as.roman(1:16)), stringsAsFactors = FALSE)
  colnames(chr_names) <- c('alignment', 'gene_data')
  
  aligned_density_data_list <- list()
  for (i in 1:nrow(gene_data)) {
    chr_temp <- chr_names[chr_names$gene_data == gene_data[i, "chr"], "alignment"]
    
    if (gene_data[i, "strand"] == '+') {
      genome_pos_min <- gene_data[i, "plus1"] - max_promoter_length
      genome_pos_max <- min(max_ORF_length + gene_data[i, "plus1"], gene_data[i, "end"])
      
      new_density_data <- subset(density_data_ref_seq_list[[chr_temp]], position >= genome_pos_min & position <= genome_pos_max)
      new_density_data$position = new_density_data$position - gene_data[i, "plus1"]
      
    } else {
      genome_pos_min <- max(gene_data[i, "start"], gene_data[i, "plus1"] - max_ORF_length)
      genome_pos_max <- gene_data[i, "plus1"] + max_promoter_length
      
      new_density_data <- subset(density_data_ref_seq_list[[chr_temp]], position >= genome_pos_min & position <= genome_pos_max)
      new_density_data$position = gene_data[i, "plus1"] - new_density_data$position

      if(nrow(new_density_data)>0) {
        new_density_data = new_density_data[nrow(new_density_data):1, ]  # revert position order (not needed)
      }
    }
    
    if(nrow(new_density_data)>0) {
      new_density_data$read_name <- paste0(new_density_data$read_name, "_gene_", i)
      new_density_data$aligned_gene_number = i
      new_density_data$aligned_gene_name = gene_data[i, "name"]
      new_density_data$aligned_gene_strand = gene_data[i, "strand"]
      
      aligned_density_data_list[[i]] <- new_density_data
    }
  }
  
  aligned_density_data <- bind_rows(aligned_density_data_list)  # makes one data.frame from the list
  aligned_density_data$read_name <- as.factor(aligned_density_data$read_name)
  
  aligned_density_data
}



filter_aligned_grand_canon <- function(aligned_density_data, pos_min = 0, pos_max = 0) {
  # removes reads that do not cover positions pos_min to pos_max
  
  density_data_starts <- aggregate(aligned_density_data$start, list(as.character(aligned_density_data$read_name)), min)
  colnames(density_data_starts) <- c("read_name", "start")
  density_data_ends <- aggregate(aligned_density_data$end, list(as.character(aligned_density_data$read_name)), max)
  colnames(density_data_ends) <- c("read_name", "end")
  density_data_starts_ends <- merge(density_data_starts, density_data_ends)
  read_names_filtered <- subset(density_data_starts_ends, start < pos_min & end > pos_max)$read_name
  aligned_density_data_filtered <- aligned_density_data[aligned_density_data$read_name %in% read_names_filtered, ]
  
  aligned_density_data_filtered
}



split_reads <- function(meth_data, read_stats, split_length_max) {
  # split reads into smaller fragments
  
  meth_data$read_name <- as.character(meth_data$read_name)
  
  meth_data_sample <- split(meth_data, meth_data$sample)
  samples <- unique(as.character(meth_data$sample))
  
  temp_meth_data_list_sample <- list()
  temp_read_stats_list_sample <- list()
  
  for(j in samples) {
    meth_data_sample_read <- split(meth_data_sample[[j]], meth_data_sample[[j]]$read_name)
    read_names <- unique(as.character(meth_data_sample[[j]]$read_name))
    num_reads <- length(unique(as.character(meth_data_sample[[j]]$read_name)))
    
    read_stats_sample <- subset(read_stats, sample == j)
    
    temp_meth_data_list_read <- list()
    temp_read_stats_list_read <- list()
    
    for(i in read_names) {
      read_stats_sample_read <- subset(read_stats_sample, read_name == i)
      if(nrow(read_stats_sample_read) == 0) {
        stop("Some reads from meth_data are not found in stat_reads file!\n")
      }
      read_length <- read_stats_sample_read$last_site - read_stats_sample_read$first_site +1
      num_fragment <- ceiling(read_length/split_length_max)
      frag_length <- floor(read_length/num_fragment)
      
      ## starting bp of each fragment
      read_start_vec <- read_stats_sample_read$first_site + c(0:(num_fragment-1))*frag_length
      if(num_fragment > 1) {
        in_group_bps <- NULL ## bp that belong to a methy/unmethy/ambig group detected in meth_data
        for(k in 1:nrow(meth_data_sample_read[[i]])) {
          in_group_bps <- c(in_group_bps, c(meth_data_sample_read[[i]]$start[k]:meth_data_sample_read[[i]]$end[k]))
        }
        all_bps <- c(read_stats_sample_read$first_site:read_stats_sample_read$last_site)
        out_group_bps <- setdiff(all_bps, in_group_bps) ## bp that do not belong to a group from meth_data
        for(k in 2:length(read_start_vec)) {
          min_s <- min(abs(read_start_vec[k]-out_group_bps))
          min_e <- min(abs(read_start_vec[k]-1-out_group_bps))
          argmin <- which.min(c(min_s, min_e))
          if(argmin == 1) {
            read_start_vec[k] <- read_start_vec[k] + min_s
          }
          else if(argmin ==2) {
            read_start_vec[k] <- read_start_vec[k] - min_s
          }
        }
      }
      
      ## ending bp of each fragment
      read_end_vec <- c(read_start_vec[-1]-1,read_stats_sample_read$last_site)
      
      ## splitting reads (update read_stats data and meth_data)
      rn <- NULL
      fs <- NULL
      ls <- NULL
      for (k in 1:num_fragment) {
        rn <- c(rn, paste0(read_stats_sample_read$read_name, "_", k))
        fs <- c(fs, read_start_vec[k])
        ls <- c(ls, read_end_vec[k])
      }
      temp_read_stats_df <- data.frame('read_name' = rn, 'first_site' = fs, 'last_site' = ls, 'sample' = j)
      temp_read_stats_df$read_name <- as.character(temp_read_stats_df$read_name)
      temp_read_stats_df$sample <- as.character(temp_read_stats_df$sample)
      # temp_read_stats_df$strand <- read_stats_sample_read$strand
      temp_read_stats_list_read[[i]] <- bind_rows(temp_read_stats_df)
      
      temp_meth_data_df <- meth_data_sample_read[[i]]
      begining_row <- NULL
      ending_row <- NULL
      for (k in 1:num_fragment) {
        begining_row <- c(begining_row, which(meth_data_sample_read[[i]]$start >= (read_start_vec[k]))[1])
      }
      ending_row <- c(begining_row[-1]-1 ,nrow(meth_data_sample_read[[i]]))
      for (k in 1:num_fragment) {
        temp_meth_data_df[begining_row[k]:ending_row[k], ]$read_name <- paste0(meth_data_sample_read[[i]][begining_row[k]:ending_row[k], ]$read_name, "_", k)
      }
      temp_meth_data_list_read[[i]] <- bind_rows(temp_meth_data_df)

    } ## end for loop on i
    
    temp_meth_data_list_sample[[j]] <- bind_rows(temp_meth_data_list_read)
    temp_read_stats_list_sample[[j]] <- bind_rows(temp_read_stats_list_read)
  }
  
  new_data <- list()
  new_data[[1]] <- bind_rows(temp_meth_data_list_sample)
  new_data[[2]] <- bind_rows(temp_read_stats_list_sample)
  
  new_data
}



filter_read_length <- function(meth_data, read_stats, split_length_min) {
  # select reads with read_length >= split_length_min
  
  if(length(setdiff(unique(read_stats$read_name), unique(meth_data$read_name))) >= 1) {
    stop("Different reads in meth_data and stat_reads files!\n")
  }
  
  read_stats_new <- read_stats[read_stats$last_site - read_stats$first_site +1 >= split_length_min, ]
  read_stats_read_name <- unique(as.character(read_stats_new$read_name))
  meth_data_new <- subset(meth_data, read_name %in% read_stats_read_name)
  
  new_data <- list()
  new_data[[1]] <- meth_data_new
  new_data[[2]] <- read_stats_new
  
  new_data
}



cut_reads <- function(meth_data, read_stats, cut_length, pos_min, pos_max, pos_offset, region_limit) {
  # cut reads before pos_min and after pos_max
  
  meth_data$read_name <- as.character(meth_data$read_name)
  
  meth_data_sample <- split(meth_data, meth_data$sample)
  samples <- unique(as.character(meth_data$sample))
  
  temp_meth_data_list_sample <- list()
  temp_read_stats_list_sample <- list()
  
  for(j in samples) {
    meth_data_sample_read <- split(meth_data_sample[[j]], meth_data_sample[[j]]$read_name)
    read_names <- unique(as.character(meth_data_sample[[j]]$read_name))
    num_reads <- length(unique(as.character(meth_data_sample[[j]]$read_name)))
    
    read_stats_sample <- subset(read_stats, sample == j)
    
    temp_meth_data_list_read <- list()
    temp_read_stats_list_read <- list()
    
    for(i in read_names) {
      read_stats_sample_read <- subset(read_stats_sample, read_name == i)
      if(nrow(read_stats_sample_read) == 0) {
        stop("Some reads from meth_data are not found in stat_reads file!\n")
      }
      
      l1 = (pos_min + pos_offset) - read_stats_sample_read$first_site
      l2 = pos_max - pos_min
      l3 = read_stats_sample_read$last_site - (pos_max + pos_offset)

      cut_start = NULL
      cut_end = NULL
      if( (l1+l3) >= (cut_length-l2) ){ ## only select reads that read_length >= cut_length
        if( l1 <= floor((cut_length-l2)/2) ) {
          cut_start = read_stats_sample_read$first_site
          cut_end = read_stats_sample_read$first_site + cut_length
        }
        else if( l3 <= floor((cut_length - l2)/2) ) {
          cut_end = read_stats_sample_read$last_site
          cut_start = read_stats_sample_read$last_site - cut_length
        }
        else {
          cut_start = pos_min + pos_offset - floor((cut_length - l2)/2)
          cut_end = cut_start + cut_length
        }

        ## check if cut positions are inside a group: have to be moved
        in_group_bps <- NULL ## bp that belong to a methy/unmethy/ambig group
        for(k in 1:nrow(meth_data_sample_read[[i]])) {
          in_group_bps <- c(in_group_bps, c(meth_data_sample_read[[i]]$start[k]:meth_data_sample_read[[i]]$end[k]))
        }
        
        if(cut_start %in% in_group_bps) {
          cut_start = meth_data_sample_read[[i]][which(cut_start <= meth_data_sample_read[[i]]$end)[1], 'start']
          begining_row <- which(meth_data_sample_read[[i]]$start >= cut_start)[1]
        }
        else {
          begining_row <- which(meth_data_sample_read[[i]]$start >= cut_start)[1]
        }
        
        if(cut_end %in% in_group_bps) {
          cut_end = meth_data_sample_read[[i]][which(cut_end <= meth_data_sample_read[[i]]$end)[1], 'end']
          ending_row <- which(cut_end <= meth_data_sample_read[[i]]$end)[1]
        }
        else {
          ending_row <- which(cut_end <= meth_data_sample_read[[i]]$end)[1] -1
        }

        ## update read_stats and meth_data
        temp_read_stats_df <- data.frame('read_name' = paste0(read_stats_sample_read$read_name, "_1"), 'first_site' = cut_start, 'last_site' = cut_end, 'sample' = j)
        temp_read_stats_df$read_name <- as.character(temp_read_stats_df$read_name)
        temp_read_stats_df$sample <- as.character(temp_read_stats_df$sample)
        temp_read_stats_list_read[[i]] <- temp_read_stats_df
  
        temp_meth_data_df <- meth_data_sample_read[[i]][begining_row:ending_row, ]
        temp_meth_data_df$read_name <- paste0(unique(meth_data_sample_read[[i]]$read_name), "_1")
        temp_meth_data_list_read[[i]] <- temp_meth_data_df
        
      } ## end if on read length selection
      
    } ## end for loop on i
    
    temp_meth_data_list_sample[[j]] <- bind_rows(temp_meth_data_list_read)
    temp_read_stats_list_sample[[j]] <- bind_rows(temp_read_stats_list_read)
  }
  
  new_data <- list()
  new_data[[1]] <- bind_rows(temp_meth_data_list_sample)
  new_data[[2]] <- bind_rows(temp_read_stats_list_sample)
  
  new_data
}



pool_canon_den <- function(grand_canon_den, sample_pools, pool_labels) {
  # sample_pools: list of vectors of samples, list entries are pooled together
  # pool_labels: vector of pool names
  
  if(missing(sample_pools)) {
    sample_pools <- list(unique(as.character(grand_canon_den$sample)))
    warning("pool_canon_den: pooling all samples.")
  }
  
  if(missing(pool_labels)) {
    pool_labels <- character(0)
    for(i in 1:length(sample_pools)) {
      pool_labels[i] <- paste0("samples ", paste(sample_pools[[i]], collapse = ", "))
    }
  }
  
  canon_den_pooled_list <- vector("list", length = length(sample_pools))
  
  for(i in 1:length(sample_pools)) {
    canon_den_pooled <- subset(grand_canon_den, sample %in% sample_pools[[i]])
    if(nrow(canon_den_pooled) > 0) {
      canon_den_pooled <- canon_den_pooled[,c('read_name', 'position', 'one_particle_density', 'sample')]
      canon_den_pooled$pool_label <- pool_labels[i]

      canon_den_pooled_list[[i]] <- canon_den_pooled
    } else {
      warning(paste0("pool_canon_den: no data to pool for samples ", sample_pools[[i]], "!\n"))
    }
  }
  
  canon_den_pooled <- bind_rows(canon_den_pooled_list)
  
  if(!is.null(canon_den_pooled)) {
    canon_den_pooled$pool_label <- factor(canon_den_pooled$pool_label, pool_labels)
    return(canon_den_pooled)
  } else {
    return(data.frame("pool_label" = NA))
  }
}



align_ss_gap_Np1 <- function(ss_gap_df, chr_ref_seqs, max_promoter_length, max_ORF_length) {
  ## align gaps to the plus_1_nucl
  ## only consider gaps that has start and end between N+1 - promotor_length and ORF_length
  
  ss_gap_df_ref_seq_list <- list()  # list of ss_gap_df with ref_seq as index
  for (ref_seq_label in chr_ref_seqs) {
    ss_gap_df_ref_seq_list = c(ss_gap_df_ref_seq_list, list(ss_gap_df[ss_gap_df$ref_seq == ref_seq_label, ]))
  }
  names(ss_gap_df_ref_seq_list) <- chr_ref_seqs
  
  gene_data = load_scer_gene_data()
  chr_names = data.frame(chr_ref_seqs, paste0("chr", as.character(1:16)), stringsAsFactors = FALSE)
  colnames(chr_names) <- c('alignment', 'gene_data')
  
  aligned_ss_gap_list <- list()
  for (i in 1:nrow(gene_data)) {
    chr_temp <- chr_names[chr_names$gene_data == gene_data[i, "chr"], "alignment"]
    
    if (gene_data[i, "strand"] == "+") {
      genome_pos_min <- gene_data[i, "aligned_pos"] - max_promoter_length
      genome_pos_max <- min(gene_data[i, "aligned_pos"] + max_ORF_length, gene_data[i, "TTS"])

      new_ss_gap_df <- subset(ss_gap_df_ref_seq_list[[chr_temp]], start >= genome_pos_min & end <= genome_pos_max)
      # new_ss_gap_df$start = new_ss_gap_df$start - gene_data[i, "aligned_pos"]
      # new_ss_gap_df$end = new_ss_gap_df$end - gene_data[i, "aligned_pos"]
      
    } else {
      genome_pos_min <- max(gene_data[i, "TTS"], gene_data[i, "aligned_pos"] - max_ORF_length)
      genome_pos_max <- gene_data[i, "aligned_pos"] + max_promoter_length

      new_ss_gap_df <- subset(ss_gap_df_ref_seq_list[[chr_temp]], start >= genome_pos_min & end <= genome_pos_max)
      # new_ss_gap_df_start <- new_ss_gap_df$start  # swapping start and end on '-' strand, so that start <= end again
      # new_ss_gap_df$start = gene_data[i, "aligned_pos"] - new_ss_gap_df$end
      # new_ss_gap_df$end = gene_data[i, "aligned_pos"] - new_ss_gap_df_start
      # 
      # if(nrow(new_ss_gap_df)>0) {
      #   new_ss_gap_df = new_ss_gap_df[nrow(new_ss_gap_df):1, ]  # revert position order (not needed)
      # }
    }
    
    if(nrow(new_ss_gap_df)>0) {
      new_ss_gap_df$read_name <- paste0(new_ss_gap_df$read_name, "_gene_", i)
      new_ss_gap_df$aligned_gene_number = i
      new_ss_gap_df$aligned_gene_name = gene_data[i, "name"]
      new_ss_gap_df$aligned_gene_strand = gene_data[i, "strand"]
      
      aligned_ss_gap_list[[i]] <- new_ss_gap_df
    }
  }
  
  aligned_ss_gap <- bind_rows(aligned_ss_gap_list)  # makes one data.frame from the list (should be faster than do.call(rbind, ...), but is still the slowest step)
  
  aligned_ss_gap
}



load_scer_selected_genes <- function() {
  
  selected_genes <- read.table("../../../external_data/nanopore_individual_genes_to_plot.tab", header = TRUE, stringsAsFactors = FALSE)
  
  return(selected_genes)
}



align_position_df_Np1 <- function(position_df, max_promoter_length, max_ORF_length, ending_point) {
  ## align positions to the plus_1_nucl
  ## consider positions that are between N+1 - promotor_length and gene TTS or N+1 + ORF_length
  
  position_df_ref_seq_list <- list()  # list of position_df with ref_seq as index, much faster than subsetting
  for (ref_seq_label in chr_ref_seqs) {
    position_df_ref_seq_list = c(position_df_ref_seq_list, list(position_df[position_df$ref_seq == ref_seq_label, ]))
  }
  names(position_df_ref_seq_list) <- chr_ref_seqs

  gene_data = load_scer_gene_data()
  chr_names = data.frame(chr_ref_seqs, paste0("chr", as.character(1:16)), stringsAsFactors = FALSE)
  colnames(chr_names) <- c('alignment', 'gene_data')
  
  aligned_position_df_list <- list()

  print('alignement started')
  for (i in 1:nrow(gene_data)) {
    chr_temp <- chr_names[chr_names$gene_data == gene_data[i, "chr"], "alignment"]
    
    if (gene_data[i, "strand"] == "+") {
      
      genome_pos_min <- gene_data[i, "aligned_pos"] - max_promoter_length
      # genome_pos_max <- min(gene_data[i, "aligned_pos"] + max_ORF_length, gene_data[i, "TTS"])
      if (ending_point == "TTS") {
        genome_pos_max <- gene_data[i, "TTS"]
      } else if (ending_point == "max_ORF_length") {
        genome_pos_max <- gene_data[i, "aligned_pos"] + max_ORF_length
      }
      
      new_position_df <- subset(position_df_ref_seq_list[[chr_temp]], position >= genome_pos_min & position <= genome_pos_max)
      new_position_df$position = new_position_df$position - gene_data[i, "aligned_pos"]

    } else {
      
      if (ending_point == "TTS") {
        genome_pos_min <- gene_data[i, "TTS"]
      } else if (ending_point == "max_ORF_length") {
        genome_pos_min <- gene_data[i, "aligned_pos"] - max_ORF_length
      }

      # genome_pos_min <- max(gene_data[i, "TTS"], gene_data[i, "aligned_pos"] - max_ORF_length)
      genome_pos_max <- gene_data[i, "aligned_pos"] + max_promoter_length
      
      new_position_df <- subset(position_df_ref_seq_list[[chr_temp]], position >= genome_pos_min & position <= genome_pos_max)
      new_position_df$position = gene_data[i, "aligned_pos"] - new_position_df$position
      
      if(nrow(new_position_df)>0) {
        new_position_df = new_position_df[nrow(new_position_df):1, ]  # revert position order (not needed)
      }
    }

    if(nrow(new_position_df)>0) {
      new_position_df$read_name <- paste0(new_position_df$read_name, "_gene_", i)
      new_position_df$aligned_gene_number = i
      new_position_df$aligned_gene_name = gene_data[i, "name"]
      new_position_df$aligned_gene_strand = gene_data[i, "strand"]

      aligned_position_df_list[[i]] <- new_position_df
    }
  }
  
  aligned_position_df <- bind_rows(aligned_position_df_list)  # makes one data.frame from the list (should be faster than do.call(rbind, ...), but is still the slowest step)
  rm(position_df_ref_seq_list)
  rm(aligned_position_df_list)

  aligned_position_df
}



filter_aligned_nucl_df <- function(aligned_nucl_df, pos_min = 0, pos_max = 0) {
  # removes reads that do not cover positions pos_min to pos_max
  nucl_df_starts <- aggregate(aligned_nucl_df$position, list(as.character(aligned_nucl_df$read_name)), min)
  colnames(nucl_df_starts) <- c("read_name", "position_min")
  nucl_df_ends <- aggregate(aligned_nucl_df$position, list(as.character(aligned_nucl_df$read_name)), max)
  colnames(nucl_df_ends) <- c("read_name", "position_max")
  nucl_df_starts_ends <- merge(nucl_df_starts, nucl_df_ends)
  read_names_filtered <- subset(nucl_df_starts_ends, position_min < pos_min & position_max > pos_max)$read_name
  aligned_nucl_df_filtered <- aligned_nucl_df[aligned_nucl_df$read_name %in% read_names_filtered, ]
  
  aligned_nucl_df_filtered
}



align_meth_data_Np1_ending_point <- function(meth_data, max_promoter_length, max_ORF_length, ending_point) {
  
  meth_data_ref_seq_list <- list()  # list of meth_data with ref_seq as index
  for (ref_seq_label in chr_ref_seqs) {
    meth_data_ref_seq_list = c(meth_data_ref_seq_list, list(meth_data[meth_data$ref_seq == ref_seq_label, ]))
  }
  names(meth_data_ref_seq_list) <- chr_ref_seqs
  
  gene_data = load_scer_gene_data()
  chr_names = data.frame(chr_ref_seqs, paste0("chr", as.character(1:16)), stringsAsFactors = FALSE)
  colnames(chr_names) <- c('alignment', 'gene_data')
  
  aligned_meth_data_list <- list()
  for (i in 1:nrow(gene_data)) {
    chr_temp <- chr_names[chr_names$gene_data == gene_data[i, "chr"], "alignment"]
    
    if (gene_data[i, "strand"] == "+") {
      genome_pos_min <- gene_data[i, "aligned_pos"] - max_promoter_length
      if (ending_point == "TTS") {
        genome_pos_max <- gene_data[i, "TTS"]
      } else if (ending_point == "max_ORF_length") {
        genome_pos_max <- gene_data[i, "aligned_pos"] + max_ORF_length
      }
      
      new_meth_data <- subset(meth_data_ref_seq_list[[chr_temp]], position >= genome_pos_min & position <= genome_pos_max)
      new_meth_data$position = new_meth_data$position - gene_data[i, "aligned_pos"]
      new_meth_data$start = new_meth_data$start - gene_data[i, "aligned_pos"]
      new_meth_data$end = new_meth_data$end - gene_data[i, "aligned_pos"]
      
    } else {
      if (ending_point == "TTS") {
        genome_pos_min <- gene_data[i, "TTS"]
      } else if (ending_point == "max_ORF_length") {
        genome_pos_min <- gene_data[i, "aligned_pos"] - max_ORF_length
      }
      
      # genome_pos_min <- max(gene_data[i, "TTS"], gene_data[i, "aligned_pos"] - max_ORF_length)
      genome_pos_max <- gene_data[i, "aligned_pos"] + max_promoter_length
      
      new_meth_data <- subset(meth_data_ref_seq_list[[chr_temp]], position >= genome_pos_min & position <= genome_pos_max)
      new_meth_data$position = gene_data[i, "aligned_pos"] - new_meth_data$position
      new_meth_data_start <- new_meth_data$start  # swapping start and end on '-' strand, so that start <= end again
      new_meth_data$start = gene_data[i, "aligned_pos"] - new_meth_data$end
      new_meth_data$end = gene_data[i, "aligned_pos"] - new_meth_data_start
      
      if(nrow(new_meth_data)>0) {
        new_meth_data = new_meth_data[nrow(new_meth_data):1, ]  # revert position order (not needed)
      }
    }
    
    if(nrow(new_meth_data)>0) {
      new_meth_data$read_name <- paste0(new_meth_data$read_name, "_gene_", i)
      new_meth_data$aligned_gene_number = i
      new_meth_data$aligned_gene_name = gene_data[i, "name"]
      new_meth_data$aligned_gene_strand = gene_data[i, "strand"]
      
      aligned_meth_data_list[[i]] <- new_meth_data
    }
  }
  
  aligned_meth_data <- bind_rows(aligned_meth_data_list)  # makes one data.frame from the list (should be faster than do.call(rbind, ...), but is still the slowest step)
  aligned_meth_data$read_name <- as.factor(aligned_meth_data$read_name)
  
  aligned_meth_data
}



calc_occupancy_by_sample <- function(nucl_df_sample, nucl_footprint) {
  ## extends the position to half_nucl_footprint from each side to get the occupancy
  
  half_nucl_footprint = floor(nucl_footprint/2)
  
  nucl_df_sample$read_name <- as.character(nucl_df_sample$read_name)
  nucl_df_sample$sample <- as.character(nucl_df_sample$sample)
  this_sample = unique(nucl_df_sample$sample)
  
  nucl_df_list <- split(nucl_df_sample, nucl_df_sample$read_name)
  read_names <- unique(as.character(nucl_df_sample$read_name))
  
  occ_df_list <- list()

  for(i in read_names) {
    min_footprint = min(nucl_df_list[[i]]$position) - half_nucl_footprint
    max_footprint = max(nucl_df_list[[i]]$position) + half_nucl_footprint

    occ_df_read <- data.frame(position = seq(min_footprint, max_footprint, 1), occu = rep(0, length(max_footprint:min_footprint)))

    for(j in nucl_df_list[[i]]$position) {
      lower_row = j - half_nucl_footprint - min_footprint + 1
      upper_row = j + half_nucl_footprint - min_footprint + 1
      # occ_df_read$occu[lower_row:upper_row] = occ_df_read$occu[lower_row:upper_row] + 1
      occ_df_read$occu[lower_row:upper_row] = 1
    }

    occ_df_list[[i]] <- occ_df_read
    occ_df_list[[i]]$read_name <- unique(as.character(i))
    occ_df_list[[i]]$sample <- unique(as.character(this_sample))
    occ_df_list[[i]]$ref_seq <- unique(nucl_df_list[[i]]$ref_seq)

    nucl_df_col_names <- names(nucl_df_list[[i]])
    if(length(nucl_df_col_names) > 4) {
      for(col_name in nucl_df_col_names[5:length(nucl_df_col_names)]) {
        occ_df_list[[i]]$col_name <- unique(nucl_df_list[[i]][col_name][,1])
        names(occ_df_list[[i]])[names(occ_df_list[[i]]) == "col_name"] <- col_name
      }
    }
  }
  
  occu_df_sample <- bind_rows(occ_df_list)
  rm(nucl_df_list)
  rm(occ_df_list)
  
  if(nrow(occu_df_sample) > 0 ) {
    occu_df_sample$read_name <- as.factor(occu_df_sample$read_name)
    occu_df_sample$sample <- as.factor(this_sample)
  }
  
  occu_df_sample 
}



calc_average_occu_single_gene_by_sample <- function(nucl_df_sample_chr, nucl_footprint, binning_length, max_promoter_length, max_ORF_length, ending_point)  {
  
  occu_df_sample <- calc_occupancy_by_sample(nucl_df_sample_chr, nucl_footprint)
  print('occu calculated')
  
  if(ending_point == "TTS") {
    aligned_genes_occu <- align_position_df_Np1(occu_df_sample, max_promoter_length, max_ORF_length, "TTS")
  } else if(ending_point == "max_ORF_length") {
    aligned_genes_occu <- align_position_df_Np1(occu_df_sample, max_promoter_length, max_ORF_length, "max_ORF_length")
  }
  rm(occu_df_sample)
  print('alignement finished')
  
  aligned_genes_occu = aligned_genes_occu[aligned_genes_occu$position %% binning_length == 0,]
  aligned_genes_occu$aligned_gene_name <- as.character(aligned_genes_occu$aligned_gene_name)
  aligned_genes_occu$aligned_gene_name <- factor(aligned_genes_occu$aligned_gene_name, levels = unique(aligned_genes_occu$aligned_gene_name))
  aligned_occu_single_gene_list <- split(aligned_genes_occu, aligned_genes_occu$aligned_gene_name)
  gene_names <- unique(as.character(aligned_genes_occu$aligned_gene_name))
  rm(aligned_genes_occu)
  
  average_gene_list <- list()
  
  # average_gene_list <- mclapply(aligned_occu_single_gene_list, function(occu_df) { ## calculates mean occu for each gene: cannot allocate memory for large data
  # 
  #   average_occu = aggregate(occu_df[, 'occu'], list(occu_df$position), mean, na.rm = TRUE)
  #   average_occu$occu_sum = aggregate(occu_df[, 'occu'], list(occu_df$position), sum, na.rm = TRUE)[, 2]
  #   colnames(average_occu) = c('position', 'occu_mean', 'occu_sum')
  #   average_occu$sample <- unique(as.character(occu_df[, "sample"]))
  #   average_occu$locus_name <- unique(as.character(occu_df[, "aligned_gene_name"]))
  # 
  #   return(average_occu)
  # }, mc.cores = num_cores)
  
  for(gene_name in gene_names) {
    # print(gene_name)
    average_occu <- NULL
    average_occu = aggregate(aligned_occu_single_gene_list[[gene_name]][, 'occu'], list(aligned_occu_single_gene_list[[gene_name]]$position), mean, na.rm = TRUE)
    average_occu$occu_sum = aggregate(aligned_occu_single_gene_list[[gene_name]][, 'occu'], list(aligned_occu_single_gene_list[[gene_name]]$position), sum, na.rm = TRUE)[, 2]
    colnames(average_occu) = c('position', 'occu_mean', 'occu_sum')
    average_occu$sample <- unique(as.character(aligned_occu_single_gene_list[[gene_name]][, "sample"]))
    average_occu$locus_name <- unique(as.character(aligned_occu_single_gene_list[[gene_name]][, "aligned_gene_name"]))
    
    average_gene_list[[gene_name]] = average_occu
  }
  
  average_occu_single_genes <- bind_rows(average_gene_list)  # makes one data.frame from the list
  rm(average_gene_list)
  
  if(nrow(average_occu_single_genes) > 0 ) {
    average_occu_single_genes$sample <- as.factor(average_occu_single_genes$sample)
    average_occu_single_genes$locus_name <- as.factor(average_occu_single_genes$locus_name)
  }
  
  average_occu_single_genes
}



calc_average_dyad_gene_ind_samples <- function(nucl_df_chr, nucl_footprint, pos_binning, max_promoter_length, max_ORF_length, ending_point)  {

  extend_dyad_df <- function(nucl_df_sample, nucl_footprint) {
    ## extends the position to half_nucl_footprint from each side
    ## score = 1 only at dyad position and 0 every where else
    
    half_nucl_footprint = floor(nucl_footprint/2)
    
    nucl_df_sample$read_name <- as.character(nucl_df_sample$read_name)
    nucl_df_sample$sample <- as.character(nucl_df_sample$sample)
    this_sample = unique(nucl_df_sample$sample)
    
    nucl_df_read <- split(nucl_df_sample, nucl_df_sample$read_name)
    read_names <- unique(as.character(nucl_df_sample$read_name))
    
    extended_dyad_df_list <- list()
    
    for(i in read_names) {
      min_footprint = min(nucl_df_read[[i]]$position) - half_nucl_footprint
      max_footprint = max(nucl_df_read[[i]]$position) + half_nucl_footprint
      
      extended_dyad_df_read <- data.frame(position = seq(min_footprint, max_footprint, 1), score = rep(0, length(max_footprint:min_footprint)))
      
      extended_dyad_df_read$score[nucl_df_read[[i]]$position - min_footprint + 1] = extended_dyad_df_read$score[nucl_df_read[[i]]$position - min_footprint + 1] + 1
      
      extended_dyad_df_list[[i]] <- extended_dyad_df_read
      extended_dyad_df_list[[i]]$read_name <- unique(as.character(i))
      extended_dyad_df_list[[i]]$sample <- this_sample
      extended_dyad_df_list[[i]]$ref_seq <- unique(nucl_df_read[[i]]$ref_seq)
    }
    
    extended_dyad_df_sample <- bind_rows(extended_dyad_df_list)
    
    if(nrow(extended_dyad_df_sample) > 0 ) {
      extended_dyad_df_sample$read_name <- as.factor(extended_dyad_df_sample$read_name)
      extended_dyad_df_sample$sample <- as.factor(this_sample)
    }
    
    extended_dyad_df_sample 
  }
  
  samples <- unique(as.character(nucl_df_chr$sample))
  aligned_genes_list <- vector(mode = "list", length = length(samples))
  average_gene_list <- vector(mode = "list", length = length(samples))
  
  for(i in 1:length(samples)) {
    nucl_df_sample <- subset(nucl_df_chr, sample == samples[i])
    nucl_df_sample$score = 1
    extended_dyad_df <- extend_dyad_df(nucl_df_sample, nucl_footprint)

    if(ending_point == "TTS") {
      aligned_genes_list[[i]] <- align_position_df_Np1(extended_dyad_df, max_promoter_length, max_ORF_length, "TTS")
    } else if(ending_point == "max_ORF_length") {
      aligned_genes_list[[i]] <- align_position_df_Np1(extended_dyad_df, max_promoter_length, max_ORF_length, "max_ORF_length")
    }
    
    aligned_genes_list[[i]]$position_binned = round(aligned_genes_list[[i]]$position/pos_binning,0)*pos_binning
    
    average_gene_list[[i]] = aggregate(aligned_genes_list[[i]][, 'score'], list(aligned_genes_list[[i]]$position_binned), mean, na.rm = TRUE)
    average_gene_list[[i]]$dyad_sum = aggregate(aligned_genes_list[[i]][, 'score'], list(aligned_genes_list[[i]]$position_binned), sum, na.rm = TRUE)[, 2]
    colnames(average_gene_list[[i]]) = c('position', 'dyad_mean', 'dyad_sum')
    average_gene_list[[i]]$sample <- unique(as.character(extended_dyad_df[, "sample"]))

  }
  
  average_gene_ind_samples <- bind_rows(average_gene_list)  # makes one data.frame from the list
  
  if(nrow(average_gene_ind_samples) > 0 ) {
    average_gene_ind_samples$sample <- factor(average_gene_ind_samples$sample, levels=samples)
  }
  
  average_gene_ind_samples
}



select_reads_based_on_gap_size <- function(mm_gap_df, nucl_footprint, max_overlap) {

  mm_gap_df$read_name <- as.character(mm_gap_df$read_name)
  mm_gap_df$sample <- as.character(mm_gap_df$sample)
  samples <- unique(as.character(mm_gap_df$sample))

  new_gap_DF_list <- list()

  for(this_sample in samples) {
    mm_gap_df_sample <- subset(mm_gap_df, sample == this_sample)
    mm_gap_df_read <- split(mm_gap_df_sample, mm_gap_df_sample$read_name)
    read_names <- unique(as.character(mm_gap_df_sample$read_name))
    num_reads <- length(unique(as.character(mm_gap_df_sample$read_name)))

    new_gap_df_list <- list()
    for(i in read_names) {
      suitable_mum_gaps <- mm_gap_df_read[[i]][mm_gap_df_read[[i]]$mm_gap >= (nucl_footprint - 2* max_overlap) & mm_gap_df_read[[i]]$mm_gap <= (2* nucl_footprint - 3* max_overlap) & mm_gap_df_read[[i]]$num_u >= 1,]
      num_suitable_mum_gaps <- nrow(suitable_mum_gaps)
      all_mm_gaps <- mm_gap_df_read[[i]]
      num_all_mm_gaps <- nrow(all_mm_gaps)
      all_mum_gaps <- mm_gap_df_read[[i]][mm_gap_df_read[[i]]$num_u >= 1,]
      num_all_mum_gaps <- nrow(all_mum_gaps)
        
      if(num_all_mum_gaps >= 0.3*num_all_mm_gaps & num_suitable_mum_gaps >= 0.5*num_all_mum_gaps) {
        new_gap_df_list[[i]] <- mm_gap_df_read[[i]][mm_gap_df_read[[i]]$num_u >= 1,]
      }
    }
    new_gap_DF_list[[this_sample]] <- bind_rows(new_gap_df_list)
  }
  
  new_gap_df <- bind_rows(new_gap_DF_list)

  if(nrow(new_gap_df) > 0 ) {
    new_gap_df$sample <- as.factor(new_gap_df$sample)
    new_gap_df$read_name <- as.factor(new_gap_df$read_name)
    new_gap_df$ref_seq <- as.factor(new_gap_df$ref_seq)
  }
  
  new_gap_df
}



calc_occu_df_from_cov_df <- function(cov_df) {
  
  arabic_number_versions <- list("chr1" = "chr01", "chr2" = "chr02", "chr3" = "chr03", "chr4" = "chr04", "chr5" = "chr05", "chr6" = "chr06", "chr7" = "chr07", "chr8" = "chr08", "chr9" = "chr09", 
                                 "chr10" = "chr10", "chr11" = "chr11", "chr12" = "chr12", "chr13" = "chr13", "chr14" = "chr14", "chr15" = "chr15", "chr16" = "chr16")
  
  for(i in 1:16){
    names(cov_df)[i] <- arabic_number_versions[[names(cov_df)[i]]]
  }
  
  occu_df_chr_list <- list()
  for(chr in names(cov_df)) {
    temp_df <- data.frame(ref_seq = chr, position = 1:length(cov_df[[chr]]), occu = cov_df[[chr]], stringsAsFactors = FALSE)
    occu_df_chr_list[[chr]] <- subset(temp_df, position %% 1 == 0)
  }
  occs_df <- bind_rows(occu_df_chr_list)
  return(occs_df)
}



calc_gene_df_of_nucl_data <- function(max_promoter_length, max_ORF_length, ending_point) {
  
  load("/scratch/ge68kuc/Dropbox/Maryam and Michael/abs_occ_paper/cov_list_Chereji_True_Zhang.RData")
  occu_df <- calc_occu_df_from_cov_df(cov_list[["Chereji2018"]])
  
  if(ending_point == "TTS") {
    aligned_genes_occu <- align_position_df_Np1(occu_df, max_promoter_length, max_ORF_length, "TTS")
  } else if(ending_point == "max_ORF_length") {
    aligned_genes_occu <- align_position_df_Np1(occu_df, max_promoter_length, max_ORF_length, "max_ORF_length")
  }
  
  aligned_genes_occu = aligned_genes_occu[aligned_genes_occu$position %% binning_length == 0,]
  aligned_genes_occu$aligned_gene_name <- as.character(aligned_genes_occu$aligned_gene_name)
  aligned_genes_occu$aligned_gene_name <- factor(aligned_genes_occu$aligned_gene_name, levels = unique(aligned_genes_occu$aligned_gene_name))
  aligned_occu_single_gene_list <- split(aligned_genes_occu, aligned_genes_occu$aligned_gene_name)
  
  average_gene_list <- list()
  
  average_gene_list <- mclapply(aligned_occu_single_gene_list, function(occu_df) { ## calculates mean occu for each gene
    
    average_occu = aggregate(occu_df[, 'occu'], list(occu_df$position), mean, na.rm = TRUE)
    colnames(average_occu) = c('position', 'occu_mean')
    average_occu$locus_name <- unique(as.character(occu_df[, "aligned_gene_name"]))
    
    return(average_occu)
  }, mc.cores = num_cores)
  
  average_occu_single_genes <- bind_rows(average_gene_list)  # makes one data.frame from the list
  
  if(nrow(average_occu_single_genes) > 0 ) {
    average_occu_single_genes$locus_name <- as.factor(average_occu_single_genes$locus_name)
  }
  
  return(average_occu_single_genes)
}



# step_function <- function(x, nucl_footprint) {
#   
#   half_nucl_footprint = floor(nucl_footprint/2)
#   h = 0
#   if(x >= -half_nucl_footprint & x <= half_nucl_footprint) {
#     
#     h = 1/(nucl_footprint)
#   }
#   
#   return(h)
# }
# 
# 
# calc_occu_chem_cleav <- function(chem_cleav_dyad_pos, nucl_footprint, binning_length, max_promoter_length, max_ORF_length, ending_point) {
#   
#   half_nucl_footprint = floor(nucl_footprint/2)
#   all_rep = unique(chem_cleav_dyad_pos$rep)
#   chem_cleav_occu <- NULL
#   
#   for(this_rep in all_rep) {
#     
#     chem_cleav_dyad_pos_rep <- subset(chem_cleav_dyad_pos, rep == this_rep)
#     gene_names <- unique(chem_cleav_dyad_pos_rep$locus)
#     
#     chem_cleav_occu_rep <- list()
#     
#     for(gene_name in gene_names) {
#       
#       chem_cleav_dyad_pos_rep_gene <- subset(chem_cleav_dyad_pos_rep, locus == gene_name)
#       
#       chem_cleav_occu_gene <- list()
#       position_vec <- c(-(max_promoter_length-half_nucl_footprint):(max_ORF_length-half_nucl_footprint))
#       occu_vec <- rep(0, length(position_vec))
# 
#       for(i in position_vec) {
#         for(j in -half_nucl_footprint:half_nucl_footprint) {
#           occu_vec[i-min(position_vec)+1] <- occu_vec[i-min(position_vec)+1] + chem_cleav_dyad_pos_rep_gene[chem_cleav_dyad_pos_rep_gene$position == (i+j),]$dyad_density*step_function(j, nucl_footprint)
#         }
#       }
#       chem_cleav_occu_gene[[gene_name]] <- data.frame(position = position_vec, occupancy = occu_vec, locus = gene_name, rep = this_rep)
#     }
#     chem_cleav_occu_rep[[this_rep]] <- bind_rows(chem_cleav_occu_gene)
#   }
# 
#   chem_cleav_occu <- bind_rows(chem_cleav_occu_rep)
#   
#   chem_cleav_occu
# }



align_read_special_gene <- function(read_stats, gene_name) {
  ## consider positions that are between N+1 and TTS
  
  gene_data = load_scer_gene_data()
  chr_names = data.frame(chr_ref_seqs, paste0("chr", as.character(1:16)), stringsAsFactors = FALSE)
  colnames(chr_names) <- c('alignment', 'gene_data')
  
  i = which(gene_data$name == gene_name)
  if(length(i) != 1) {
    warning("None or more than one gene with that name found.")
    read_names <- character(0)
  } else {
    chr_temp <- chr_names[chr_names$gene_data == gene_data[i, "chr"], "alignment"]
    
    if (gene_data[i, "strand"] == "+") {
      
      genome_pos_min <- gene_data[i, "aligned_pos"]
      genome_pos_max <- gene_data[i, "TTS"]
      
      # meth_data_locus <- subset(meth_data, ref_seq == chr_temp & position >= genome_pos_min & position <= genome_pos_max)
      read_stats_locus <- subset(read_stats, ref_seq == chr_temp & last_site >= genome_pos_min & first_site <= genome_pos_max)
      
    } else {
      
      genome_pos_min <- gene_data[i, "TTS"]
      genome_pos_max <- gene_data[i, "aligned_pos"]
      
      # meth_data_locus <- subset(meth_data, ref_seq == chr_temp & position >= genome_pos_min & position <= genome_pos_max)
      read_stats_locus <- subset(read_stats, ref_seq == chr_temp & first_site <= genome_pos_max & last_site >= genome_pos_min)
    }
    
    # meth_read_names <- unique(meth_data_locus$read_name)
    read_names <- unique(read_stats_locus$read_name)
  }
  
  read_names
}


align_reads_one_gene <- function(read_stats, gene_name) {
  ## consider positions that are between N+1 and TTS
  
  gene_data = load_scer_gene_data()
  chr_names = data.frame(chr_ref_seqs, paste0("chr", as.character(1:16)), stringsAsFactors = FALSE)
  colnames(chr_names) <- c('alignment', 'gene_data')
  
  i = which(gene_data$name == gene_name)
  chr_temp <- chr_names[chr_names$gene_data == gene_data[i, "chr"], "alignment"]
  
  if (gene_data[i, "strand"] == "+") {
    
    genome_pos_min <- gene_data[i, "aligned_pos"]
    genome_pos_max <- gene_data[i, "TTS"]
    
    read_stats_locus <- subset(read_stats, ref_seq == chr_temp & last_site >= genome_pos_min & first_site <= genome_pos_max)
  } else {
    
    genome_pos_min <- gene_data[i, "TTS"]
    genome_pos_max <- gene_data[i, "aligned_pos"]
    
    read_stats_locus <- subset(read_stats, ref_seq == chr_temp & first_site <= genome_pos_max & last_site >= genome_pos_min)
  }
  
  read_stats_locus
}