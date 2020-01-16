plot_hist_with_origin <- function(df, bins=21, xlim, filename_addition="", na.rm=TRUE, title, y_axis_density=TRUE) {  # df: first column for histogram, second column is origin
  if(nrow(df)==0) return()
  cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  df_name <- strsplit(as.character(substitute(df)), "\\[")[[2]]
  p1 <- ggplot(df, aes_(x=as.name(names(df)[1]))) + ggtitle(df_name)
  if(y_axis_density) p1 <- p1 + geom_histogram(aes(y=..density..), bins=bins, na.rm = na.rm)
  else p1 <- p1 + geom_histogram(bins=bins, na.rm = na.rm)
  p2 <- ggplot(df, aes_(x=as.name(names(df)[1]), fill=as.name(names(df)[2]))) +
    geom_histogram(bins=bins, position = "fill", na.rm = na.rm) + scale_fill_manual(values = rep(cbPalette, 4)) + ylab("origin of counts")
  if(ncol(df) > 2 && "sample" %in% colnames(df)) {
    p1 <- p1 + facet_wrap(~sample)
    p2 <- p2 + facet_wrap(~sample)
  }
  if(ncol(df) > 2 && "pool_label" %in% colnames(df)) {
    p1 <- p1 + facet_wrap(~pool_label)
    p2 <- p2 + facet_wrap(~pool_label)
  }
  if(missing("xlim")) {
    xlim <- round(quantile(df[, 1], c(0.01, 0.99), na.rm = na.rm))
  }
  p1 <- p1 + xlim(xlim = xlim)
  p2 <- p2 + xlim(xlim = xlim)
  if(!missing("title")) {
    p1 <- p1 + ggtitle(title)
  }
  if(length(unique(df[, 2])) == 1) {
    p12 <- p1 + xlab(paste0(as.name(names(df)[1]), " (only ", unique(df[, 2]), ")"))
  } else {
    p1 <- p1 + ggplot_larger_plot
    p2 <- p2 + ggplot_larger_plot
    p12 <- plot_grid(p1, p2 + theme(legend.position = "none"), ncol=1, align="v") + theme(aspect.ratio = 1)
    p12 <- plot_grid(p12, get_legend(p2), rel_widths = c(3, 1))
  }
  if(!missing("title")) {
    p12 <- p12 + ggtitle(title)
  }
  filename <- paste0(df_name, "_hist_origin_", names(df)[1], "_", names(df)[2], filename_addition)
  p12
}


# multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
#   # Multiple plot function (from the 'Cookbook for R')
#   #
#   # ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
#   # - cols:   Number of columns in layout
#   # - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#   #
#   # If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
#   # then plot 1 will go in the upper left, 2 will go in the upper right, and
#   # 3 will go all the way across the bottom.
#   #
#   
#   require(grid)
#   
#   # Make a list from the ... arguments and plotlist
#   plots <- c(list(...), plotlist)
#   
#   numPlots = length(plots)
#   
#   # If layout is NULL, then use 'cols' to determine layout
#   if (is.null(layout)) {
#     # Make the panel
#     # ncol: Number of columns of plots
#     # nrow: Number of rows needed, calculated from # of cols
#     layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
#                      ncol = cols, nrow = ceiling(numPlots/cols))
#   }
#   
#   if (numPlots==1) {
#     print(plots[[1]])
#   } else {
#     # Set up the page
#     grid.newpage()
#     pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
#     
#     # Make each plot, in the correct location
#     for (i in 1:numPlots) {
#       # Get the i,j matrix positions of the regions that contain this subplot
#       matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
#       
#       print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
#                                       layout.pos.col = matchidx$col))
#     }
#   }
# }


plot_2d_densities <- function(df, bins = 30, filename_addition = "", na.rm = TRUE, title = "") {
  if(nrow(df)==0) return()
  n <- ncol(df)
  plot_list <- list()
  list_ind <- 0
  plot_cols <- max(2,ncol(df)-2)
  plot_rows <- ceiling(n*(n-1)/2/plot_cols)
  barheights <- c(NaN, 0.07, 0.04, 0.025, 0.02)
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      list_ind <- list_ind + 1
      plot_list[[list_ind]] <- ggplot(df, aes_(x=as.name(names(df)[i]), y=as.name(names(df)[j]))) +
        geom_bin2d(bins = bins, na.rm = na.rm) + ggplot_larger_plot +
        scale_fill_gradientn(trans = "log", colours = my_colours, breaks = 4^seq(0,8)) +
        theme(legend.key.width = unit(0.006, "npc"))
    }
  }
  
  p <- plot_grid(plotlist = plot_list, ncol = plot_cols, align = "v")
  
  if(title != "") {
    p_title <- ggdraw() + draw_label(title, fontface = 'bold', size = 10)
    p <- plot_grid(p_title, p, ncol = 1, rel_heights = c(0.1, 1))
  }
  
  df_name <- strsplit(as.character(substitute(df)), "\\[")[[2]]
  filename <- paste0(df_name, "_2d_densities", filename_addition)
  ggsave(paste0(plot_folder, filename, plot_string), plot=p, width=1.5*1.5*plot_cols, height=1.5*plot_rows, units="in")

  p
}


# plot_2d_densities_multiplot <- function(df, bins = 30, filename) {
#   n <- ncol(df)
#   plot_list <- list()
#   list_ind <- 0
#   plot_cols <- max(2,ncol(df)-2)
#   plot_rows <- ceiling(n*(n-1)/2/plot_cols)
#   barheights <- c(NaN, 0.07, 0.04, 0.025, 0.02)
#   for (i in 1:(n-1)) {
#     for (j in (i+1):n) {
#       list_ind <- list_ind + 1
#       plot_list[[list_ind]] <- ggplot(df, aes_(x=as.name(names(df)[i]), y=as.name(names(df)[j]))) + geom_bin2d(bins = bins) +
#         scale_fill_gradientn(trans = "log", colours = my_colours, breaks = 4^seq(0,8)) +
#         background_grid(major = "xy", minor = "none") +
#         theme(axis.text = element_text(size = 6)) +  # changes axis labels
#         theme(axis.title = element_text(size = 6)) +  # change axis titles
#         theme(text = element_text(size = 6)) +  # this will change all text size (except geom_text)
#         theme(legend.key.width = unit(0.01, "npc"), legend.key.height = unit(barheights[plot_rows], "npc"))
#     }
#   }
#   multiplot(plotlist = plot_list, cols = plot_cols)
#   if (!missing(filename)) {
#     dev.copy(pdf, width=6, height=5, paste0(plot_folder, filename, ".pdf"))
#     invisible(dev.off())
#   }
# }


# plot_2d_densities_facet_wrap <- function(df) {
#   n <- ncol(df)
#   df_plot <- data.frame(x = numeric(0), y = numeric(0), x_ind = integer(0), y_ind = integer(0), plot_name = character(0))
#   for (i in 1:(n-1)) {
#     for (j in (i+1):n) {
#       plot_name <- paste0("y = ", colnames(df)[j], " vs x = ", colnames(df)[i])
#       df_plot <- rbind(df_plot, data.frame(df[, i], df[, j], i, j, plot_name))
#     }
#   }
#   colnames(df_plot) <- c("x", "y", "x_ind", "y_ind", "plot_name")
#   p <- ggplot(df_plot, aes(x=x, y=y)) + 
#     geom_bin2d(bins = 15) +
#     scale_fill_gradientn(trans = "log", colours = my_colours, breaks = c(1,4,16,64)) +
#     background_grid(major = "xy", minor = "none") +
#     facet_wrap(~plot_name, ncol = 2, scales = "free") +
#     theme(axis.text = element_text(size = 6)) +  # changes axis labels
#     theme(axis.title = element_text(size = 6)) +  # change axis titles
#     theme(text = element_text(size = 6)) +  # this will change all text size (except geom_text)
#     theme(legend.key.width = unit(0.02, "npc"))
#   p
# }


plot_ordered_configs <- function(meth_data_matrix, order_vector, row_order, title = 'Individual reads', ylabel = 'reads', position_offset = 0, legend_title = "state") {
  if (!missing('order_vector') && missing('row_order')) row_order <- order(order_vector)
  else if (missing('order_vector') && missing('row_order')) row_order <- 1:nrow(meth_data_matrix)
  else if (!missing('order_vector') && !missing('row_order')) stop('Can only give order_vector or row_order!')
  
  df_meth_matrix <- melt(as.matrix(t(meth_data_matrix[row_order, ])))
  colnames(df_meth_matrix) <- c('position','read','hit_state')
  df_meth_matrix$position <- df_meth_matrix$position * pos_binning + position_offset
  df_meth_matrix$hit_state <- factor(df_meth_matrix$hit_state)
  
  if(any(meth_data_matrix[!is.na(meth_data_matrix)] == 0)) {
    my_colours <- c("blue", "yellow", "red", "grey")
    my_labels <- c("unmethylated", "ambiguous", "methylated", "no CpG site")
  } else {
    my_colours <- c("blue", "red", "grey")
    my_labels <- c("unmethylated", "methylated", "no CpG site")
  }
  
  if (missing('order_vector')) {
    p <- ggplot(df_meth_matrix, aes(x = position, y = read)) + 
      geom_raster(aes(fill=hit_state)) + scale_fill_manual(name = legend_title, values = my_colours, labels = my_labels) +
      ggtitle(title) +
      ylab(ylabel)
  } else {
    my_ylabels <- seq(0,1,0.1)
    my_ylabels <- my_ylabels[my_ylabels>min(order_vector) & my_ylabels<max(order_vector)]
    my_ybreaks <- c()
    for(p in my_ylabels) {
      my_ybreaks <- c(my_ybreaks, min(which(sort(order_vector)>=p)))
    }
    p <- ggplot(df_meth_matrix, aes(x = position, y = read)) + 
      geom_raster(aes(fill=hit_state)) + scale_fill_manual(name = legend_title, values = my_colours, labels = my_labels) +
      ggtitle(title) +
      ylab(ylabel) + scale_y_continuous(breaks = my_ybreaks, labels = my_ylabels)
  }
  p
}


plot_individual_reads <- function(meth_data, read_order, reverse_order = FALSE,
                                  pos_min, pos_max, pos_offset = 0, minus_strand = FALSE,
                                  site_radius = 4, group_footprinting = TRUE, clustering = "off",
                                  filter_length = 10,  # applies to clustering = "bin", "sum" and "sum_and_bin"
                                  my_colours = c("lightskyblue", "grey87", "firebrick", "white"),
                                  continuous_scale = FALSE) {
  # read_order contains two columns read names (in the intended order) and order values (used for y-ticks)
  # clustering can be "off", "as_is", "sum", "bin"
  # in case (continuous_scale = False) uses hit_state and in case (continuous_scale = TRUE) uses log_lik_ratio
  
  if(any(meth_data$start > meth_data$end)) {
    stop("meth_data$start > meth_data$end for some reads!")
  }
  
  if(clustering != "off" & ( !missing(read_order) || !missing(reverse_order) )) {
    stop("read_order and reverse_order cannot be used together with clustering != off")
  }
  
  if(!missing(pos_min) && !missing(pos_max)) {
    meth_data <- subset(meth_data, end + site_radius < pos_max + pos_offset & start - site_radius > pos_min + pos_offset)  # only take sites within the plotting area
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
  
  configurations <- data.frame(position = rep(pos_min:pos_max, num_reads), read = unlist(lapply(1:num_reads, rep, (pos_max-pos_min+1))), hit_state = rep(NaN, num_reads*(pos_max-pos_min+1)), log_lik_ratio = rep(NaN, num_reads*(pos_max-pos_min+1)))
  # configurations$position <- configurations$position + genes_to_plot$TSS  
  hit_states_list <- list()
  hit_states_pure_list <- list()
  log_lik_ratios_list <- list()

  meth_data$read_name <- as.character(meth_data$read_name)  # split crashes if this is left a factor
  meth_data_read_list <- split(meth_data, meth_data$read_name)  # faster than using subset every time
  
  for(i in 1:num_reads) {
    meth_data_read <- meth_data_read_list[[read_names[i]]]
    hit_states <- rep(NaN, pos_max - pos_min + 1)
    log_lik_ratios <- rep(NaN, pos_max - pos_min + 1)
    for(j in -site_radius:site_radius) {
      hit_states[floor(meth_data_read$position) + j - pos_min + 1 - pos_offset] <- meth_data_read$hit_state
      log_lik_ratios[floor(meth_data_read$position) + j - pos_min + 1 - pos_offset] <- meth_data_read$log_lik_ratio
    }
    
    if(group_footprinting) {  # enlarge footprint of groups
      meth_data_groups <- subset(meth_data_read, group_size > 1)
      if(nrow(meth_data_groups) > 0) {
        for(j in 1:nrow(meth_data_groups)) {
          footprint <- ( (meth_data_groups[j, "start"]-site_radius) : (meth_data_groups[j, "end"]+site_radius) ) - pos_min + 1 - pos_offset
          if(any(footprint < 0)) browser()
          hit_states[footprint] <- meth_data_groups[j, "hit_state"]
          log_lik_ratios[footprint] <- meth_data_groups[j, "log_lik_ratio"]
        }
      }
    }
    ## force the log_lik_ratios to be between -5 and 5 (for a better plot)
    log_lik_ratios <- ifelse(log_lik_ratios >= 5, 5, log_lik_ratios)
    log_lik_ratios = ifelse(log_lik_ratios <= -5, -5, log_lik_ratios)
      
    hit_states_list[[i]] <- hit_states
    log_lik_ratios_list[[i]] <- log_lik_ratios

    if(clustering == "bin") {
      hit_states_pure <- rep(NaN, pos_max - pos_min + 1)  # without group_footprinting and site radius
      hit_states_pure[floor(meth_data_read$position) - pos_min + 1 - pos_offset] <- meth_data_read$hit_state
      hit_states_pure_list[[i]] <- hit_states_pure
    }
  }
  
  configurations$hit_state <- do.call(c, hit_states_list)
  configurations$log_lik_ratio <- do.call(c, log_lik_ratios_list)
  
  if (continuous_scale == FALSE) {
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
      my_labels <- c("unmethylated", "ambiguous", "methylated", "no CpG site")
    } else {
      my_colours <- my_colours[c(1, 3, 4)]
      my_labels <- c("unmethylated", "methylated", "no CpG site")
    }
    
    configurations$hit_state <- factor(configurations$hit_state)
  }
  
  if(minus_strand && pos_offset > 0) {
    configurations$position <- -configurations$position 
  }
    
  if (continuous_scale == FALSE) {
    p <- ggplot(configurations, aes(x = position, y = read)) + geom_raster(aes(fill=hit_state)) + 
      scale_fill_manual(name = "site state", values = my_colours, labels = my_labels) +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  } else if (continuous_scale == TRUE) {
    p <- ggplot(configurations, aes(x = position, y = read)) + geom_raster(aes(fill=log_lik_ratio)) + 
      scale_fill_viridis(na.value="white") +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  }
  
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

plot_individual_reads_daan <- function(FirstSiteMin, FirstSiteMax, ReadLength, SampleName = " ", RefSeq = " ", Barcode = " ", GeneLength, PromotorLength) {
  ####filter data for length####
  #Gene to be plotted has to be in genes_to_plot
  #1: select region of interest in meth_data
  read_stats_chr <- subset(read_stats_combined, ref_seq %in% chr_ref_seqs)  # only chromosomes
  #1: select region of interest in meth_data with min and max values for starting position
  read_stats_chr <- subset(read_stats_chr, first_site>=FirstSiteMin)
  read_stats_chr <- subset(read_stats_chr, first_site<=FirstSiteMax)
  #2: select how long you want the read to be maximally (no sense to make longer than region of interest)
  read_stats_chr <- subset(read_stats_chr, dist_first_to_last_site >= ReadLength)
  #3: select which sample and on which chromosome it is
  read_stats_chr <- subset(read_stats_chr, sample == SampleName)
  read_stats_chr <- subset(read_stats_chr, ref_seq == RefSeq)
  #4: make new vector with the names of the reads
  long_reads <- read_stats_chr$read_name
  #5.0: load meth_data_chr. can be used to subset for barcode (sample) and chromosome (ref_seq)
  meth_data_chr <-  load_meth_data(barcode_subset = c("barcode03"), ref_seq_subset = "")
  #5: use these read_names to subset in meth_data_chr
  meth_data_chr <- meth_data_chr[meth_data_chr$read_name %in%  long_reads,]
  nrow(meth_data_chr)
  
  for(i in 1:nrow(genes_to_plot)) {
    locus <- genes_to_plot$locus_name[i]
    cat("### ", locus, "aka", genes_to_plot$gene[i],"  \n")
    
    chr_roman_to_arabic <- list("chrI" = "chr01", "chrII" = "chr02", "chrIII" = "chr03", "chrIV" = "chr04", "chrV" = "chr05", "chrVI" = "chr06", "chrVII" = "chr07", "chrVIII" = "chr08", "chrIX" = "chr09", "chrX" = "chr10", "chrXI" = "chr11", "chrXII" = "chr12", "chrXIII" = "chr13", "chrXIV" = "chr14", "chrXV" = "chr15", "chrXVI" = "chr16")
    # gene length ultimately determines the length of the x-axis. reads will be cut off at max length of gene if required
    gene_length <- GeneLength
    # promotor length is how much the x-axis goes into the negative
    promoter_length <- PromotorLength
    site_stats_locus_list <- list()
    TSS <- genes_to_plot$TSS[i]
    chr <- chr_roman_to_arabic[[genes_to_plot$chr[i]]]
    
    if(genes_to_plot$strand[i] == 1) {
      site_stats_locus <- subset(site_stats_chr, ref_seq==chr & position >= TSS - promoter_length & position <= TSS + gene_length)
      site_stats_locus$pos_TSS <- site_stats_locus$position - TSS
      
      site_stats_locus_pooled <- pool_site_stats(site_stats_locus, sample_pools_chr, pool_labels_chr)
      site_stats_locus_pooled$pos_TSS <- site_stats_locus_pooled$position - TSS
      
      read_stats_locus <- subset(read_stats_chr, ref_seq==chr & last_site >= TSS - promoter_length & first_site <= TSS + gene_length)
      minus_strand <- FALSE
      pos_min <- -promoter_length  # wrt TSS offset
      pos_max <- gene_length
    } else if(genes_to_plot$strand[i] == -1) {
      site_stats_locus <- subset(site_stats_chr, ref_seq==chr & position <= TSS + promoter_length & position >= TSS - gene_length)
      site_stats_locus$pos_TSS <- TSS - site_stats_locus$position
      
      site_stats_locus_pooled <- pool_site_stats(site_stats_locus, sample_pools_chr, pool_labels_chr)
      site_stats_locus_pooled$pos_TSS <- TSS - site_stats_locus_pooled$position
      
      read_stats_locus <- subset(read_stats_chr, ref_seq==chr & first_site <= TSS + promoter_length & last_site >= TSS - gene_length)
      minus_strand <- TRUE
      pos_min <- -gene_length  # wrt TSS offset
      pos_max <- promoter_length
    } else {
      stop("Invalid strand")
    }
  }
  
  if(nrow(meth_data_chr) > 0) {
    pp(plot_individual_reads(meth_data_chr, clustering = "sum_and_bin", pos_min = pos_min, pos_max = pos_max,
                             pos_offset = TSS, minus_strand =  minus_strand, group_footprinting = TRUE, filter_length = 150) +
         ggtitle(paste0("Individual reads: ", SampleName)))
    #changed clustering from "sum" to "sum_and_bin"
    # rm(meth_data_temp)
  } else {
    cat("\nNo individual reads to plot.\n\n")
  }
}


plot_log_lik_ratio_densities <- function(meth_data, title_addition="", plot_group_sizes = FALSE) {
  if(!plot_group_sizes) {
    p <- ggplot(meth_data, aes(x=log_lik_ratio))
  } else {
    meth_data <- subset(meth_data, group_size <= 3)
    meth_data$group_size <- factor(meth_data$group_size)
    p <- ggplot(meth_data, aes(x=log_lik_ratio, colour=group_size)) + ggplot_colours
  }
  p <- p + geom_density(n = 2^12) + xlim(c(-100, 100)) + coord_cartesian(xlim = c(-20, 20)) + 
    geom_vline(xintercept = c(-log_lik_ratio_threshold, log_lik_ratio_threshold), colour = "red") + 
    ggtitle(paste0("Log-likelihood-ratio density ", title_addition)) + facet_wrap(~sample)
  
  p
}


plot_coverage_ecdf <- function(coverage, title_addition ="", my_xlim, my_xticks) {
  if(nrow(coverage)==0) return()
  p <- ggplot(coverage, aes(x = read_count, colour = sample)) + stat_ecdf() + ylab("empirical cdf") + ggtitle(paste0("Coverage ", title_addition)) + ggplot_colours
  if(missing(my_xlim)) {
    my_xlim <- c(0, quantile(coverage$read_count, 0.99))
  } 
  p <- p + coord_cartesian(xlim = my_xlim)
  if(!missing(my_xticks)) {
    p <- p + scale_x_continuous(breaks = my_xticks)
  }
  p
}


plot_channel_stats <- function(channel_stats) {
  p_list <- list(qplot(channel_stats$log10_aligned_bases,  main='Log10 of aligned bases') + xlab("log10 aligned bases") + ggplot_larger_plot, 
                 qplot(channel_stats$reads,  main='Number of reads') + xlab("reads") + ggplot_larger_plot, 
                 qplot(channel_stats$site_fraction_ambiguous,  main='Average channel ambiguity') + xlab("ambiguity") + ggplot_larger_plot,
                 qplot(channel_stats$site_call_fraction_methylated, main='Average channel methylation') + xlab("methylation") + ggplot_larger_plot,
                 qplot(channel_stats$alignment_clipping_loss/channel_stats$bases_incl_hardclipping, main='Alignment clipping loss of total bases') + xlab("alignment clipping loss") + ggplot_larger_plot)
  plot_grid(plotlist = p_list)
}


pp <- function(p) {  # used to plot ggplots inside loops
  if(!is.null(p)) {
    cat("\n")
    print(p)
    cat("  \n")  # needed for Rmarkdown to render properly
  }
}


knit_from_tools <- function() {
  dir <- getwd()
  rmarkdown::render(input="../../tools/nanopore_analysis.Rmd", output_file=paste0(gsub('^.*nanopore/\\s*|\\s*/analysis.*$', '', dir), ".html"), intermediates_dir = dir, output_dir=dir, knit_root_dir=dir)
}


my_knit <- function() {
  dir <- getwd()
  rmarkdown::render(input="nanopore_analysis.Rmd", output_file=paste0(gsub('^.*nanopore/\\s*|\\s*/analysis.*$', '', dir), ".html"))
}
