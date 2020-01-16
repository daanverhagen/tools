plot_densities  <- function(df, aesthetic, xlim, adjust, title_addition = "", na.rm = TRUE) {
  p <- ggplot(df, aes_(x=as.name(aesthetic))) +
    geom_density(adjust = adjust, na.rm = na.rm) +
    labs(title=paste0("Density of ", aesthetic, " < ", max(xlim), " ", title_addition))
   
  if(!missing(xlim)) {
    p <- p + coord_cartesian(xlim = xlim) + expand_limits(x = xlim)
    #+ scale_y_log10()
  }

  p
}

plot_histogram  <- function(df, aesthetic, xlim, bins, title_addition = "", na.rm = TRUE) {
  p <- ggplot(df, aes_(x=as.name(aesthetic))) +
    geom_histogram(aes(y=..density..), bins = bins, na.rm = na.rm)
  
  if(!missing(xlim)) {
    p <- p  + labs(title=paste0("Density of ", aesthetic, " < ", max(xlim), " ", title_addition)) + 
      xlim(xlim) + expand_limits(x = xlim)
  }
  
  p
}

plot_2d_densities_1plot <- function(df, aesthetic1, aesthetic2, bins, title_addition = "", na.rm = TRUE) {
  p <- ggplot(df, aes_(x=as.name(aesthetic1), y=as.name(aesthetic2))) +
    geom_bin2d(bins = bins , na.rm = na.rm) + 
    scale_fill_gradientn(trans = "log", colours = my_colours, breaks = 4^seq(0,8)) +
    theme(legend.key.width = unit(0.006, "npc"))

  p
}


plot_individual_reads_with_nucl <- function(meth_data, nucl_df, read_order, reverse_order = FALSE,
                                  pos_min, pos_max, pos_offset = 0, minus_strand = FALSE,
                                  site_radius = 4, nucl_radius = 4, group_footprinting = TRUE, clustering = "off",
                                  filter_length = 10,  # applies to clustering = "bin", "sum" and "sum_and_bin"
                                  my_colours = c("lightskyblue", "grey87", "firebrick", "black", "white")) {
  # read_order contains two columns read names (in the intended order) and order values (used for y-ticks)
  # clustering can be "off", "as_is", "sum", "bin"
  
  if(nrow(nucl_df) == 0) {
    stop("no reads in nucl_df!\n")
  }
  if(any(meth_data$start > meth_data$end)) {
    stop("meth_data$start > meth_data$end for some reads!")
  }
  
  if(clustering != "off" & ( !missing(read_order) || !missing(reverse_order) )) {
    stop("read_order and reverse_order cannot be used together with clustering != off")
  }
  
  if(!missing(pos_min) && !missing(pos_max)) {
    meth_data <- subset(meth_data, end + site_radius < pos_max + pos_offset & start - site_radius > pos_min + pos_offset)  # only take sites within the plotting area
    nucl_df <- subset(nucl_df, position + nucl_radius < pos_max + pos_offset & position - nucl_radius > pos_min + pos_offset)  # only take sites within the plotting area
    
  }
  
  if(!missing(read_order)) {
    if(reverse_order) {
      read_order <- read_order[order(-read_order[, 2]), ]
    } else {
      read_order <- read_order[order(read_order[, 2]), ]
    }
    meth_data <- subset(meth_data, read_name %in% read_order$read_name)
    read_names <- intersect(as.character(read_order$read_name), unique(as.character(meth_data$read_name)))  # intersect keeps the order of the first input set
  } else {
    read_names <-  unique(as.character(meth_data$read_name))
  }
  
  num_reads <- length(read_names)
  if(num_reads == 0) {
    # warning("plot_individual_reads: no read in meth_data! returning NULL...\n")
    # return("test")  # doesn't work for some reason
    stop("plot_individual_reads: no read in meth_data!\n")
  }
  
  if(clustering != "off" && num_reads == 1) {
    clustering <- "off"
    warning("Need more than one read for clustering. Clustering turned off.\n")
  }
  
  if(missing(pos_min)) {
    pos_min <- min(meth_data$start) - site_radius - pos_offset
  }
  if(missing(pos_max)) {
    pos_max <- max(meth_data$end) + site_radius - pos_offset
  }
  
  configurations <- data.frame(position = rep(pos_min:pos_max, num_reads), read = unlist(lapply(1:num_reads, rep, (pos_max-pos_min+1))), hit_state = rep(NaN, num_reads*(pos_max-pos_min+1)))
  hit_states_list <- list()
  hit_states_pure_list <- list()
  
  meth_data$read_name <- as.character(meth_data$read_name)  # split crashes if this is left a factor
  meth_data_read_list <- split(meth_data, meth_data$read_name)  # faster than using subset every time
  nucl_df$read_name <- as.character(nucl_df$read_name)  # split crashes if this is left a factor
  nucl_df_read_list <- split(nucl_df, nucl_df$read_name)  # faster than using subset every time
  
  for(i in 1:num_reads) {
    meth_data_read <- meth_data_read_list[[read_names[i]]]
    hit_states <- rep(NaN, pos_max - pos_min + 1)

    for(j in -site_radius:site_radius) {
      hit_states[floor(meth_data_read$position) + j - pos_min + 1 - pos_offset] <- meth_data_read$hit_state
    }

    if(group_footprinting) {  # enlarge footprint of groups
      meth_data_groups <- subset(meth_data_read, group_size > 1)
      if(nrow(meth_data_groups) > 0) {
        for(j in 1:nrow(meth_data_groups)) {
          footprint <- ( (meth_data_groups[j, "start"]-site_radius) : (meth_data_groups[j, "end"]+site_radius) ) - pos_min + 1 - pos_offset
          if(any(footprint < 0)) browser()
          hit_states[footprint] <- meth_data_groups[j, "hit_state"]
        }
      }
    }
    
    nucl_df_read <- nucl_df_read_list[[read_names[i]]]
    if(length(nucl_df_read$position) > 0) {
      for(j in -nucl_radius:nucl_radius) {
        hit_states[floor(nucl_df_read$position) + j - pos_min + 1 - pos_offset] <- 2
      }
    }
    
    hit_states_list[[i]] <- hit_states
    
    if(clustering == "bin") {
      hit_states_pure <- rep(NaN, pos_max - pos_min + 1)  # without group_footprinting and site radius
      hit_states_pure[floor(meth_data_read$position) - pos_min + 1 - pos_offset] <- meth_data_read$hit_state
      hit_states_pure_list[[i]] <- hit_states_pure
    }
  }
  
  configurations$hit_state <- do.call(c, hit_states_list)
  
  bin_positions <- function(state_matrix, pos_binning = 10) {
    # to be applied to a pure state_matrix (with site_radius = 1, no group footprinting)
    # supposed to decrease the vector size of the hit_state vector by binning for better clustering results
    # when pos_binning <= 10, this binning is just course graining, as neighboring sites are at least 11 bases appart
    
    state_matrix[is.na(state_matrix)] <- 0  # set positions without a site to 0, same as ambiguous sites
    state_matrix_binned <- matrix(0, nrow = nrow(state_matrix), ncol = ceiling(ncol(state_matrix) / pos_binning))
    for(i in 0:(pos_binning-1)) {
      index_vec <- seq(from = 1, to = ncol(state_matrix), by = pos_binning)
      index_vec <- i + index_vec
      if(all(index_vec <= ncol(state_matrix))) {
        state_matrix_binned <- state_matrix_binned + state_matrix[, index_vec]
      } else {
        index_vec <- index_vec[index_vec <= ncol(state_matrix)]
        state_matrix_binned[, -ncol(state_matrix_binned)] <- state_matrix_binned[, -ncol(state_matrix_binned)] + state_matrix[, index_vec]
      }
    }
    state_matrix_binned
  }
  
  sum_positions <- function(state_matrix, filter_length = 10) {
    # to be applied to a state_matrix (with site_radius >= 1, group footprinting on/off)
    # supposed to smooth out the methylation signal and increase the range of sites for better clustering results
    
    state_matrix[is.na(state_matrix)] <- 0  # set positions without a site to 0, same as ambiguous sites
    state_matrix_summed <- matrix(0, nrow = nrow(state_matrix), ncol = ncol(state_matrix)-(filter_length-1))
    for(i in 0:(filter_length-1)) {
      index_vec <- 1:ncol(state_matrix_summed)
      state_matrix_summed <- state_matrix_summed + state_matrix[, i+index_vec]
    }
    state_matrix_summed
  }
  
  if(clustering != "off" && length(read_names) > 1) {
    if(clustering == "as_is") {
      state_matrix <- matrix(data = configurations$hit_state, nrow = num_reads, byrow = TRUE)
      state_matrix_unchanged <- state_matrix
      state_matrix[is.na(state_matrix)] <- 0  # set positions without a site to 0, same as ambiguous sites
      row_order_clustered <- hclust(dist(state_matrix, method = "manhattan"))$order  # for many and long reads dist is really slow
      configurations$hit_state <- as.vector(t(state_matrix_unchanged[row_order_clustered, ]))
    } else if(clustering == "bin") {
      state_matrix <- matrix(data = unlist(hit_states_pure_list), nrow = num_reads, byrow = TRUE)
      state_matrix_unchanged <- matrix(data = configurations$hit_state, nrow = num_reads, byrow = TRUE)
      state_matrix_binned <- bin_positions(state_matrix, pos_binning = filter_length)
      row_order_clustered <- hclust(dist(state_matrix_binned, method = "manhattan"))$order  # binning to speed up dist
      configurations$hit_state <- as.vector(t(state_matrix_unchanged[row_order_clustered, ]))
    } else if(clustering == "sum") {
      state_matrix <- matrix(data = configurations$hit_state, nrow = num_reads, byrow = TRUE)
      state_matrix_unchanged <- state_matrix
      state_matrix_summed <- sum_positions(state_matrix, filter_length = filter_length)
      row_order_clustered <- hclust(dist(state_matrix_summed, method = "manhattan"))$order  # for many and long reads dist is really slow
      configurations$hit_state <- as.vector(t(state_matrix_unchanged[row_order_clustered, ]))
    } else if(clustering == "sum_and_bin") {
      state_matrix <- matrix(data = configurations$hit_state, nrow = num_reads, byrow = TRUE)
      state_matrix_unchanged <- state_matrix
      state_matrix_summed <- bin_positions(sum_positions(state_matrix, filter_length = filter_length), pos_binning = filter_length)
      row_order_clustered <- hclust(dist(state_matrix_summed, method = "manhattan"))$order  # binning the summed state_matrix speeds up dist
      configurations$hit_state <- as.vector(t(state_matrix_unchanged[row_order_clustered, ]))
    }
  }
  
  if(all(!is.na(configurations))) {
    my_colours <- my_colours[c(1, 2, 3)]
    my_labels <- c("unmethylated", "ambiguous / no site", "methylated")
  } else if(any(configurations[!is.na(configurations)] == 0)) {
    my_labels <- c("unmethylated", "ambiguous", "methylated", "nucleosome", "no CpG site")
  } else {
    my_colours <- my_colours[c(1, 3, 4)]
    my_labels <- c("unmethylated", "methylated", "no CpG site")
  }
  
  configurations$hit_state <- factor(configurations$hit_state)
  
  if(minus_strand && pos_offset > 0) {
    configurations$position <- -configurations$position 
  }
  
  p <- ggplot(configurations, aes(x = position, y = read)) + geom_raster(aes(fill=hit_state)) +
    scale_fill_manual(name = "site state", values = my_colours, labels = my_labels) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

  if(!missing(read_order)) {
    my_ybreaks <- seq(1, max(configurations$read), max(1, floor(max(configurations$read) / 6)))
    my_ylabels <- sprintf("%.2f", read_order[my_ybreaks, 2])
    if(any(duplicated(my_ylabels))) {
      my_ylabels <- sprintf("%.3f", read_order[my_ybreaks, 2])
    }
    p <- p + scale_y_continuous(breaks = my_ybreaks, labels = my_ylabels) + ylab(colnames(read_order)[2])
  }

  if(minus_strand && pos_offset == 0) {
    p <- p + scale_x_reverse()
  }
  
  p
}


plot_individual_reads_1_particle_den <- function(density_data, read_order, reverse_order = FALSE,
                                            pos_min, pos_max, pos_offset = 0, minus_strand = FALSE,
                                            clustering = "off",
                                            filter_length = 10  # applies to clustering = "bin", "sum" and "sum_and_bin"
                                            ) {
  # read_order contains two columns read names (in the intended order) and order values (used for y-ticks)
  # clustering can be "off", "as_is", "sum", "bin"

  if(clustering != "off" & ( !missing(read_order) || !missing(reverse_order) )) {
    stop("read_order and reverse_order cannot be used together with clustering != off")
  }

  if(!missing(pos_min) && !missing(pos_max)) {
    density_data <- subset(density_data, bp < pos_max + pos_offset & bp > pos_min + pos_offset)  # only take sites within the plotting area
  }

  if(!missing(read_order)) {
    if(reverse_order) {
      read_order <- read_order[order(-read_order[, 2]), ]
    } else {
      read_order <- read_order[order(read_order[, 2]), ]
    }
    density_data <- subset(density_data, read_name %in% read_order$read_name)
    read_names <- intersect(as.character(read_order$read_name), unique(as.character(density_data$read_name)))  # intersect keeps the order of the first input set
  } else {
    read_names <-  unique(as.character(density_data$read_name))
  }

  if(clustering != "off" && length(read_names) <= 1) {
    warning("Need more than one read for clustering. Clustering turned off.\n")
  }

  if(missing(pos_min)) {
    pos_min <- min(density_data$bp) - pos_offset
  }
  if(missing(pos_max)) {
    pos_max <- max(density_data$bp) - pos_offset
  }

  num_reads <- length(read_names)
  if(num_reads == 0) {
    warning("plot_individual_reads: no reads in density_data! returning NULL...\n")
    return(NULL)
  }

  # bin_positions <- function(state_matrix, pos_binning = 10) {
  #   # to be applied to a pure state_matrix (with site_radius = 1, no group footprinting)
  #   # supposed to decrease the vector size of the hit_state vector by binning for better clustering results
  #   # when pos_binning <= 10, this binning is just course graining, as neighboring sites are at least 11 bases appart
  # 
  #   state_matrix[is.na(state_matrix)] <- 0  # set positions without a site to 0, same as ambiguous sites
  #   state_matrix_binned <- matrix(0, nrow = nrow(state_matrix), ncol = ceiling(ncol(state_matrix) / pos_binning))
  #   for(i in 0:(pos_binning-1)) {
  #     index_vec <- seq(from = 1, to = ncol(state_matrix), by = pos_binning)
  #     index_vec <- i + index_vec
  #     if(all(index_vec <= ncol(state_matrix))) {
  #       state_matrix_binned <- state_matrix_binned + state_matrix[, index_vec]
  #     } else {
  #       index_vec <- index_vec[index_vec <= ncol(state_matrix)]
  #       state_matrix_binned[, -ncol(state_matrix_binned)] <- state_matrix_binned[, -ncol(state_matrix_binned)] + state_matrix[, index_vec]
  #     }
  #   }
  #   state_matrix_binned
  # }
  # 
  # sum_positions <- function(state_matrix, filter_length = 10) {
  #   # to be applied to a state_matrix (with site_radius >= 1, group footprinting on/off)
  #   # supposed to smooth out the methylation signal and increase the range of sites for better clustering results
  # 
  #   state_matrix[is.na(state_matrix)] <- 0  # set positions without a site to 0, same as ambiguous sites
  #   state_matrix_summed <- matrix(0, nrow = nrow(state_matrix), ncol = ncol(state_matrix)-(filter_length-1))
  #   for(i in 0:(filter_length-1)) {
  #     index_vec <- 1:ncol(state_matrix_summed)
  #     state_matrix_summed <- state_matrix_summed + state_matrix[, i+index_vec]
  #   }
  #   state_matrix_summed
  # }
  # 
  # if(clustering != "off" && length(read_names) > 1) {
  #   if(clustering == "as_is") {
  #     state_matrix <- matrix(data = density_data$density, nrow = num_reads, byrow = TRUE)
  #     state_matrix_unchanged <- state_matrix
  #     state_matrix[is.na(state_matrix)] <- 0  # set positions without a site to 0, same as ambiguous sites
  #     row_order_clustered <- hclust(dist(state_matrix, method = "manhattan"))$order  # for many and long reads dist is really slow
  #     density_data$density <- as.vector(t(state_matrix_unchanged[row_order_clustered, ]))
  #   } else if(clustering == "bin") {
  #     state_matrix <- matrix(data = unlist(hit_states_pure_list), nrow = num_reads, byrow = TRUE)
  #     state_matrix_unchanged <- matrix(data = density_data$density, nrow = num_reads, byrow = TRUE)
  #     state_matrix_binned <- bin_positions(state_matrix, pos_binning = filter_length)
  #     row_order_clustered <- hclust(dist(state_matrix_binned, method = "manhattan"))$order  # binning to speed up dist
  #     density_data$density <- as.vector(t(state_matrix_unchanged[row_order_clustered, ]))
  #   } else if(clustering == "sum") {
  #     state_matrix <- matrix(data = density_data$density, nrow = num_reads, byrow = TRUE)
  #     state_matrix_unchanged <- state_matrix
  #     state_matrix_summed <- sum_positions(state_matrix, filter_length = filter_length)
  #     row_order_clustered <- hclust(dist(state_matrix_summed, method = "manhattan"))$order  # for many and long reads dist is really slow
  #     density_data$density <- as.vector(t(state_matrix_unchanged[row_order_clustered, ]))
  #   } else if(clustering == "sum_and_bin") {
  #     state_matrix <- matrix(data = density_data$density, nrow = num_reads, byrow = TRUE)
  #     state_matrix_unchanged <- state_matrix
  #     state_matrix_summed <- bin_positions(sum_positions(state_matrix, filter_length = filter_length), pos_binning = filter_length)
  #     row_order_clustered <- hclust(dist(state_matrix_summed, method = "manhattan"))$order  # binning the summed state_matrix speeds up dist
  #     density_data$density <- as.vector(t(state_matrix_unchanged[row_order_clustered, ]))
  #   }
  # }
  # 
  # density_data$density <- factor(density_data$density)
  # 
  # if(minus_strand && pos_offset > 0) {
  #   density_data$bp <- -density_data$bp
  # }

  p <- ggplot(density_data, aes(x = bp, y = read_name)) + geom_raster(aes(fill=one_particle_density)) +
    # scale_fill_manual(name = "one particle density") #+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  
  if(!missing(read_order)) {
    my_ybreaks <- seq(1, max(density_data$read_name), max(1, floor(max(density_data$read_name) / 6)))
    p <- p + scale_y_continuous(breaks = my_ybreaks) + ylab(colnames(read_order)[2])
  }

  if(minus_strand && pos_offset == 0) {
    p <- p + scale_x_reverse()
  }
  
  p
}



plot_one_read <- function(meth_data, pos_min, pos_max, pos_offset = 0, minus_strand = FALSE,
                          site_radius = 4, 
                          my_colours = c("lightskyblue", "grey87", "firebrick", "white")) {
  
  if(any(meth_data$start > meth_data$end)) {
    stop("meth_data$start > meth_data$end for some reads!")
  }
  
  if(!missing(pos_min) && !missing(pos_max)) {
    meth_data <- subset(meth_data, end + site_radius < pos_max + pos_offset & start - site_radius > pos_min + pos_offset)  # only take sites within the plotting area
  }
  
  read_name <-  unique(as.character(meth_data$read_name))
  
  num_reads <- length(read_name)
  if(num_reads == 0 || num_reads > 1) {
    stop("no read in meth_data!\n")
  }
  if(num_reads > 1) {
    stop("more than one read in meth_data!\n")
  }

  if(missing(pos_min)) {
    pos_min <- min(meth_data$start) - site_radius - pos_offset
  }
  if(missing(pos_max)) {
    pos_max <- max(meth_data$end) + site_radius - pos_offset
  }
  
  print(length(which(meth_data$hit_state == 1)))
  
  configurations <- data.frame(position = (pos_min - site_radius):(pos_max + site_radius), read = rep(1, pos_max-pos_min+1+2*site_radius), hit_state = rep(NaN, pos_max-pos_min+1+2*site_radius))

  hit_states <- rep(NaN, pos_max - pos_min + 1 + 2* site_radius)
  
  for(i in 1:nrow(meth_data)) {
    hit_states[(meth_data$start[i]- site_radius):(meth_data$end[i]+site_radius) - pos_min + 1 - pos_offset] <- meth_data$hit_state[i]
  }
  
  # for(j in -site_radius:site_radius) {
  #   hit_states[floor(meth_data_read$position) + j - pos_min + 1 - pos_offset] <- meth_data_read$hit_state
  #   log_lik_ratios[floor(meth_data_read$position) + j - pos_min + 1 - pos_offset] <- meth_data_read$log_lik_ratio
  # }
  
  configurations$hit_state <- hit_states
  

  if(all(!is.na(configurations))) {
    my_colours <- my_colours[c(1, 2, 3)]
    my_labels <- c("unmethylated", "ambiguous / no site", "methylated")
  } else if(any(configurations[!is.na(configurations)] == 0)) {
    my_labels <- c("unmethylated", "ambiguous", "methylated", "no CpG site")
  } else {
    my_colours <- my_colours[c(1, 3, 4)]
    my_labels <- c("unmethylated", "methylated", "no CpG site")
  }
  
  configurations$hit_state <- factor(configurations$hit_state)

  if(minus_strand && pos_offset > 0) {
    configurations$position <- -configurations$position 
  }
  
  p <- ggplot(configurations, aes(x = position, y = read)) + geom_raster(aes(fill=hit_state)) + 
    scale_fill_manual(name = "site state", values = my_colours, labels = my_labels) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

  if(minus_strand && pos_offset == 0) {
    p <- p + scale_x_reverse()
  }

  p
}



plot_chem_cleav_data <- function(chem_cleav_data, locus, gene, promoter_length, gene_length, colname, title_addition, plot_string) {
  ## plots dyad_density/occupancy data for 3 replicates and the avg values in a grid
  p <- ggplot(chem_cleav_data, aes(x=position)) + 
       facet_grid(rows = vars(rep)) + 
       geom_line(aes_string(y=colname), colour="blue") +
       ylab(title_addition) + xlab("Genome position (TSS = 0)") +
       ggtitle(paste0("Chem cleavage data [Chereji et al. 2018] for locus ", locus, ' (gene ', gene, ')')) +
       theme(strip.text.y = element_text(angle = 0)) +
       scale_x_continuous(breaks = seq(-promoter_length, gene_length + 100, 200))
  ggsave(paste0(plot_folder, "chem_cleav_", colname, "_", locus, "_", gene, plot_string), plot = p, width = 20, height = 14, units = "cm")

  p
}



plot_chem_cleav_data_overlap <- function(chem_cleav_data, locus, gene, promoter_length, gene_length, colname, title_addition, alpha_value, plot_string) {
  ## plots the overlap of dyad_density/occupancy data for 3 replicates and the avg values
  p <- ggplot(chem_cleav_data, aes(x = position)) + 
       geom_line(aes(y = chem_cleav_data[,2], col = rep), alpha = alpha_value) + 
       ylab(title_addition) + xlab("Genome position (TSS = 0)") +
       ggtitle(paste0("Chem cleavage data [Chereji et al. 2018] for locus", locus, ' (gene ', gene, ')')) +
       labs(color='Chemi-cleav data') +
       scale_x_continuous(breaks = seq(-promoter_length, gene_length + 100, 200))
  ggsave(paste0(plot_folder, "chem_cleav_", colname, "2_", locus, "_", gene, plot_string), plot = p, width = 20, height = 14, units = "cm")
  
  p
}


plot_large_gap_fraction <- function(gap_df) {
  
  calc_large_gap_fraction_by_ref_seq <- function(gap_df, length_threshold) {
    
    sum_all_gaps_df = aggregate(gap_df[,1],by=list(gap_df$ref_seq), sum)
    names(sum_all_gaps_df) <- c('ref_seq', 'sum_all_gaps')
    
    gap_df_large <- gap_df[gap_df$gap_length >= length_threshold,]
    sum_large_gaps_df = aggregate(gap_df_large[,1],by=list(gap_df_large$ref_seq), sum)
    names(sum_large_gaps_df) <- c('ref_seq', 'sum_large_gaps')
    
    
    gap_sum = merge(sum_all_gaps_df, sum_large_gaps_df, by="ref_seq")
    gap_sum$frac_large_gaps <- gap_sum$sum_large_gaps/gap_sum$sum_all_gaps
    
    return(as.data.frame(gap_sum))
  }
  
  gap_summry_list <- list()
  for (large_length in seq(120, 180, 10)) {
    gap_summry_list[[large_length]] <- calc_large_gap_fraction_by_ref_seq(gap_df, large_length)
    gap_summry_list[[large_length]]$large_length <- large_length
  }
  
  gap_summary_df <- bind_rows(gap_summry_list)
  
  p <- ggplot(gap_summary_df, aes(x = large_length, y = frac_large_gaps)) + geom_line() + geom_point() + facet_wrap(~ref_seq) +
    ylab("Genome fraction") + xlab("Large gap threshold") +
    ggtitle("Genome fraction covered by large gaps")
  
  p
}

