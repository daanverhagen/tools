combine_cols_for_plot <- function(df, new_col_names, col_types) {  # restructure df to plot multiple cols of different types marked by for example different color
  
  # example: df with cols x_1, x_2 to be restructured to col x and col type with entries "type_1" for x_1 values and "type_2" for x_2 values. 
  # Then use new_col_names = c("x", "type") and col_types = c("x_1" = "type_1", "x_2" = "type_2").
  
  if(length(new_col_names) != 2) stop("Need two new column names!")
  
  df_list <- list()
  for(col in names(col_types)) {
    df_new_cols <- data.frame(df[, col], col_types[col], row.names = NULL, stringsAsFactors = FALSE)
    colnames(df_new_cols) <- new_col_names
    df_list[[col]] <- cbind(df, df_new_cols)
  }
  
  df_new <- bind_rows(df_list)
  df_new[, names(col_types)] <- NULL
  
  return(df_new)
}


helper_plot <- function(... , col=1, alpha.f=0.2, pch=19) {
  if(class(col) == "character") {
    plot(... ,  pch=pch, col=adjustcolor(col, alpha.f=alpha.f))
  } else {
    plot(... ,  pch=pch, col=adjustcolor(palette()[col], alpha.f=alpha.f))
  }
}


compare_two_occs_dfs_old <- function(sample_1, df_1, sample_2, df_2, method = "occ_mean_in_radius", min_coverage = 20, plot_folder = NA, 
                                     mean_occ_goal = 0.75, plot_types = c("as_is", "normalized", "rescaled"), no_title = FALSE, my_alpha = NA, 
                                     use_pdf = FALSE, use_ggplot = TRUE, full_info = !no_title) {
  
  # sample_1 should be the sample with less sites, usually the RE sample if RE is compared to BS
  # sites of sample_2 within given radius of a fixed site of sample_1 are averaged
  
  colnames(df_1)[grepl("mean_occ", colnames(df_1))] <- "occ"
  colnames(df_2)[grepl("mean_occ", colnames(df_2))] <- "occ"
  
  if("counts" %in% colnames(df_1)) {
    cat("Removed", round(100*sum(df_1$counts < min_coverage)/nrow(df_1)), "% low coverage sites in sample 1.\n")
    df_1 <- subset(df_1, counts >= min_coverage)
  } else {
    cat("No counts column in df_1, min_coverage not checked.\n")
  }
  
  if(grepl("^occ_mean_in_radius_[0-9]*$", method)) {
    radius = as.numeric(sub("[a-z_]*", "", method))
    
    if("counts" %in% colnames(df_2)) {
      cat("Removed", round(100*sum(df_2$counts < min_coverage)/nrow(df_2)), "% low coverage sites in sample 2.\n")
      df_2 <- subset(df_2, counts >= min_coverage)
    } else {
      cat("No counts column in df_2, min_coverage not checked.\n")
    }
    
    if(radius == 0) {
      
      chr_pos_1 <- with(df_1, paste0(chr, "_", pos))
      chr_pos_2 <- with(df_2, paste0(chr, "_", pos))
      
      chr_pos_1_good <- with(df_1, chr_pos_1[!is.na(occ)])
      chr_pos_2_good <- with(df_2, chr_pos_2[!is.na(occ)])
      
      chr_pos_valid <- intersect(chr_pos_1_good, chr_pos_2_good)
      occs_1 <- with(df_1, occ[chr_pos_1 %in% chr_pos_valid])
      occs_2 <- with(df_2, occ[chr_pos_2 %in% chr_pos_valid])
      
    } else {  # take sites of df_1 and calc the occ mean of df_2 near these sites
      
      pos_valid_list <- list()
      occs_1_list <- list()
      occs_2_radius_mean_list <- list()
      for(this_chr in unique(df_1$chr)) {
        
        df_1_chr <- subset(df_1, chr == this_chr & !is.na(occ))
        df_2_chr <- subset(df_2, chr == this_chr & !is.na(occ))
        
        pos_1_good <- df_1_chr$pos
        pos_2_good <- df_2_chr$pos
        
        occs_1 <- df_1_chr$occ
        
        occs_2_vec <- rep(NA, max(c(pos_1_good + radius, pos_2_good)))
        occs_2_vec[pos_2_good] <- df_2_chr$occ
        
        occs_2_radius_mean <- sapply(pos_1_good, function(x) mean(occs_2_vec[x+(-radius:radius)], na.rm = TRUE))
        
        pos_valid <- pos_1_good[is.finite(occs_2_radius_mean)]
        occs_2_radius_mean <- occs_2_radius_mean[is.finite(occs_2_radius_mean)]
        
        pos_valid_list[[this_chr]] <- pos_valid
        occs_1_list[[this_chr]] <- subset(df_1_chr, pos %in% pos_valid)$occ
        occs_2_radius_mean_list[[this_chr]] <- occs_2_radius_mean
      }
      occs_1 <- unlist(occs_1_list)
      occs_2 <- unlist(occs_2_radius_mean_list)
    }
  } else if(grepl("^count_mean_in_radius_[0-9]*$", method)) {
    radius = as.numeric(sub("[a-z_]*", "", method))
    if(radius == 0) stop("Use occ_mean_in_radius_0 instead.")
    
    if(!("counts" %in% colnames(df_2))) {
      stop("df_2 needs to have an occ and counts column")
    
    } else {  # take sites of df_1 and calc the count mean of df_2 near these sites
      pos_valid_list <- list()
      occs_1_list <- list()
      occs_2_radius_mean_list <- list()
      counts_2_radius_sum_list <- list()
      for(this_chr in unique(df_1$chr)) {
        
        df_1_chr <- subset(df_1, chr == this_chr & !is.na(occ))
        df_2_chr <- subset(df_2, chr == this_chr & !is.na(occ))  # this coverage will be checked after summing counts within radius
        
        pos_1_good <- df_1_chr$pos
        pos_2_good <- df_2_chr$pos
        
        occs_1 <- df_1_chr$occ
        
        counts_2_vec <- rep(NA, max(c(pos_1_good + radius, pos_2_good)))
        counts_2_vec[pos_2_good] <- df_2_chr$counts
        
        hits_2_vec <- rep(NA, max(c(pos_1_good + radius, pos_2_good)))
        hits_2 <- df_2_chr$counts * df_2_chr$occ
        if(any(abs(hits_2 - round(hits_2)) > 1e-6)) stop("occs_2 and counts_2 do not multiply to integer!")
        hits_2_vec[pos_2_good] <- round(hits_2)
        
        counts_2_radius_sum <- sapply(pos_1_good, function(x) sum(counts_2_vec[x+(-radius:radius)], na.rm = TRUE))
        
        occs_2_radius_mean <- sapply(pos_1_good, function(x) sum(hits_2_vec[x+(-radius:radius)], na.rm = TRUE) / sum(counts_2_vec[x+(-radius:radius)], na.rm = TRUE))
        
        #cat("Removed", round(100*sum(counts_2_radius_sum < min_coverage)/length(counts_2_radius_sum)), "% low coverage sites in sample 2.\n")
        
        #occs_2_radius_mean[counts_2_radius_sum < min_coverage] <- NA
        
        good_indeces <- counts_2_radius_sum > 0 & is.finite(occs_2_radius_mean)
        
        pos_valid <- pos_1_good[good_indeces]
        occs_2_radius_mean <- occs_2_radius_mean[good_indeces]
        counts_2_radius_sum <- counts_2_radius_sum[good_indeces]
        
        pos_valid_list[[this_chr]] <- pos_valid
        occs_1_list[[this_chr]] <- subset(df_1_chr, pos %in% pos_valid)$occ
        occs_2_radius_mean_list[[this_chr]] <- occs_2_radius_mean
        counts_2_radius_sum_list[[this_chr]] <- counts_2_radius_sum
      }
      occs_1 <- unlist(occs_1_list)
      occs_2 <- unlist(occs_2_radius_mean_list)
      counts_2 <- unlist(counts_2_radius_sum_list)
      
      good_indeces <- counts_2 >= min_coverage
      cat("Removed", round(100*sum(counts_2 < min_coverage)/length(counts_2)), "% low coverage sites in sample 2.\n")
      occs_1 <- occs_1[good_indeces]
      occs_2 <- occs_2[good_indeces]
      counts_2 <- counts_2[good_indeces]
    }
  } else if(grepl("^count_mean_in_radius_[0-9]*_meth$", method)) {  # this method is symmetric in swapping sample_1 with sample_2 when comparing GpC with CpG
    radius = as.numeric(sub("[a-z_]*", "", sub("_meth", "", method)))
    stop("Not implemented properly, yet!")
    if(radius == 0) stop("Use occ_mean_in_radius_0 instead.")
    
    if(!("counts" %in% colnames(df_2))) {
      stop("df_2 needs to have an occ and counts column")
      
    } else {  # take sites of df_1 and calc the count mean of df_2 near these sites
      pos_valid_list <- list()
      occs_1_list <- list()
      occs_2_radius_mean_list <- list()
      counts_2_radius_sum_list <- list()
      for(this_chr in unique(df_1$chr)) {
        
        df_1_chr <- subset(df_1, chr == this_chr & !is.na(occ))
        df_2_chr <- subset(df_2, chr == this_chr & !is.na(occ))  # this coverage will be checked after summing counts within radius
        
        pos_1_good <- df_1_chr$pos
        pos_2_good <- df_2_chr$pos
        
        occs_1 <- df_1_chr$occ
        
        counts_2_vec <- rep(NA, max(c(pos_1_good + radius, pos_2_good)))
        counts_2_vec[pos_2_good] <- df_2_chr$counts
        
        hits_2_vec <- rep(NA, max(c(pos_1_good + radius, pos_2_good)))
        hits_2 <- df_2_chr$counts * df_2_chr$occ
        if(any(abs(hits_2 - round(hits_2)) > 1e-6)) stop("occs_2 and counts_2 do not multiply to integer!")
        hits_2_vec[pos_2_good] <- round(hits_2)
        
        counts_2_radius_sum <- sapply(pos_1_good, function(x) sum(counts_2_vec[x+(-radius:radius)], na.rm = TRUE))
        
        occs_2_radius_mean <- sapply(pos_1_good, function(x) sum(hits_2_vec[x+(-radius:radius)], na.rm = TRUE) / sum(counts_2_vec[x+(-radius:radius)], na.rm = TRUE))
        
        #cat("Removed", round(100*sum(counts_2_radius_sum < min_coverage)/length(counts_2_radius_sum)), "% low coverage sites in sample 2.\n")
        
        #occs_2_radius_mean[counts_2_radius_sum < min_coverage] <- NA
        
        good_indeces <- counts_2_radius_sum > 0 & is.finite(occs_2_radius_mean)
        
        pos_valid <- pos_1_good[good_indeces]
        occs_2_radius_mean <- occs_2_radius_mean[good_indeces]
        counts_2_radius_sum <- counts_2_radius_sum[good_indeces]
        
        pos_valid_list[[this_chr]] <- pos_valid
        occs_1_list[[this_chr]] <- subset(df_1_chr, pos %in% pos_valid)$occ
        occs_2_radius_mean_list[[this_chr]] <- occs_2_radius_mean
        counts_2_radius_sum_list[[this_chr]] <- counts_2_radius_sum
      }
      occs_1 <- unlist(occs_1_list)
      occs_2 <- unlist(occs_2_radius_mean_list)
      counts_2 <- unlist(counts_2_radius_sum_list)
      
      good_indeces <- counts_2 >= min_coverage
      cat("Removed", round(100*sum(counts_2 < min_coverage)/length(counts_2)), "% low coverage sites in sample 2.\n")
      occs_1 <- occs_1[good_indeces]
      occs_2 <- occs_2[good_indeces]
      counts_2 <- counts_2[good_indeces]
    }
  } else {
    stop("method not implemented.")
  }
  
  if(mean_occ_goal == "sample_1") {
    gamma_value_1 <- 1
    mean_occ_goal <- mean(occs_1)
  } else if(is.numeric(mean_occ_goal)) {
    gamma_value_1 <- optimize_gamma(occs_1, mean_occ_goal)
  } else {
    stop("mean_occ_goal not valid.")
  }
  gamma_value_2 <- optimize_gamma(occs_2, mean_occ_goal)
  
  occs_1_rescaled <- rescale_in_zero_to_one(occs_1, gamma_value_1)
  occs_2_rescaled <- rescale_in_zero_to_one(occs_2, gamma_value_2)
  
  delta = mean(abs(occs_1 - occs_2))
  delta_2 = mean((occs_1 - occs_2)^2)
  
  beta = mean(abs(occs_1/mean(occs_1) - occs_2/mean(occs_2)))
  beta_2 = mean((occs_1/mean(occs_1) - occs_2/mean(occs_2))^2)
  
  alpha = mean(abs(occs_1_rescaled - occs_2_rescaled))
  alpha_2 = mean((occs_1_rescaled - occs_2_rescaled)^2)
  
  r_value = cor(occs_1, occs_2)
  r_value_rescaled = cor(occs_1_rescaled, occs_2_rescaled)
  
  d_value = delta + abs(mean(occs_1) - mean(occs_2))
  
  num_comparisons <- length(occs_1)
  
  if(!is.na(plot_folder)) {
    
    if(is.na(my_alpha)) my_alpha <- min(0.25, 2500/num_comparisons)
    
    if(!use_ggplot) {
      
      if(use_pdf) {
        pdf(paste0(plot_folder, "scatterplot_", method, "_" , sample_2, "_vs_", sample_1, "_min_cov_", min_coverage, ".pdf"), width=length(plot_types)*2.75, height=2.75, pointsize=6)
      } else {
        png(paste0(plot_folder, "scatterplot_", method, "_" , sample_2, "_vs_", sample_1, "_min_cov_", min_coverage, ".png"), width=length(plot_types)*1200, height=1200, res=220)
      }
      par(mfrow=c(1, length(plot_types)))
      
      if("as_is" %in% plot_types) {
        occ_limits <- c(0, 1)
        this_title <- ifelse(no_title, "", paste0(sample_2, " (", signif(mean(occs_2), 3), ") vs\n", sample_1, " (", signif(mean(occs_1), 3), ")"))
        helper_plot(occs_1, occs_2, col="black", alpha.f=my_alpha, pch=20, xlim=occ_limits, ylim=occ_limits, main=this_title,
                    xlab=paste0(sample_1, " site occupancy"), ylab=paste0(sample_2, " site occupancy"))
        lines(x=occ_limits, y=occ_limits, col="blue", lty=2, lwd=1)
        
        mtext(side=3, text=paste0("         radius = ", radius, ", ", num_comparisons, " sites"))
        mtext(side=3, adj=0, text=paste0("<|x-y|> = ", sprintf("%.3f", round(delta, 3))))
        mtext(side=3, adj=1, text=paste0("r = ", sprintf("%.3f", round(r_value, 3))))
      }
      if("normalized" %in% plot_types) {
        occ_limits <- c(0, max(1/mean(occs_1), 1/mean(occs_2)))
        helper_plot(occs_1/mean(occs_1), occs_2/mean(occs_2), col="black", alpha.f=my_alpha, pch=20, xlim=occ_limits, ylim=occ_limits, main=paste0(sample_2, " vs ", sample_1),
                    xlab=paste0(sample_1, " normalized site occupancy"), ylab=paste0(sample_2, " normalized site occupancy"))
        lines(x=occ_limits, y=occ_limits, col="blue", lty=2, lwd=1)
        
        mtext(side=3, adj=0, text=paste0("mean_occ_goal * <|x/mean(x)-y/mean(y)|> = ", sprintf("%.3f", round(mean_occ_goal * beta, 3))))
        mtext(side=3, adj=1, text=paste0("r = ", sprintf("%.3f", round(r_value, 3))))
      }
      if("rescaled" %in% plot_types) {
        occ_limits <- c(0, 1)
        helper_plot(occs_1_rescaled, occs_2_rescaled, col="black", alpha.f=my_alpha, pch=20, xlim=occ_limits, ylim=occ_limits, main=paste0(sample_2, " vs ", sample_1),
                    xlab=paste0(sample_1, " rescaled site occupancy"), ylab=paste0(sample_2, " rescaled site occupancy"))
        lines(x=occ_limits, y=occ_limits, col="blue", lty=2, lwd=1)
        
        mtext(side=3, adj=0, text=paste0("<|x'-y'|> = ", sprintf("%.3f", round(alpha, 3))))
        mtext(side=3, adj=1, text=paste0("r_rescaled = ", sprintf("%.3f", round(r_value_rescaled, 3))))
      }
      
      dev.off()
    
    } else {
    
      my_alpha <- max(0.01, my_alpha)
      
      this_title <- paste0(sample_2, " vs ", sample_1)
      
      occs_df <- data.frame(occs_1 = occs_1, occs_2 = occs_2, occs_1_rescaled = occs_1_rescaled, occs_2_rescaled = occs_2_rescaled)
      
      geom_text_size <- theme_paper[[1]]$axis.text$size
      
      occ_limits <- c(0, 1)
      p <- ggplot(occs_df, aes(x=occs_1, y=occs_2)) + 
        geom_abline(slope = 1, intercept = 0, linetype = "dashed",  color = "blue", size = 0.4) +
        geom_point(alpha = my_alpha, size = 0.4, shape = 16, stroke = 0) + 
        scale_x_continuous(labels = percent_no_sign, limits = (occ_limits)) + scale_y_continuous(labels = percent_no_sign, limits = (occ_limits)) +
        xlab(paste0(sample_1, " absolute occupancy / %")) + ylab(paste0(sample_2, " abs. occupancy / %")) + theme_paper
        # draw_label(label = paste0("<|x-y|> = ", format(delta, digits=3), 
        #                           ",  <x> = ", format(mean(occs_df$occs_1), digits=3), 
        #                           ",  <y> = ", format(mean(occs_df$occs_2), digits=3)), x=0.01, y=1.035*occ_limits[2], 
        #            size=geom_text_size, hjust=0) +
        # draw_label(label = paste0("radius = ", radius, ",  ", num_comparisons, " sites,  r = ", format(r_value, digits=3)), 
        #            x=0.99, y=-0.015*occ_limits[2], size=geom_text_size, hjust=1)

      plot_infos_1 <- paste0("<|x-y|> = ", sprintf("%.3f", round(delta, 3)), 
                             ",  <x> = ", sprintf("%.3f", round(mean(occs_df$occs_1), 3)), 
                             ",  <y> = ", sprintf("%.3f", round(mean(occs_df$occs_2), 3)), 
                             ",  d = ", sprintf("%.3f", round(d_value, 3)),
                             ",\nradius = ", radius, ",  ", num_comparisons, " sites,  r = ", sprintf("%.3f", round(r_value, 3)))
      
      plot_infos_2 <- paste0(num_comparisons, " sites,  <|x-y|> = ", sprintf("%.3f", round(delta, 3)), 
                             ",  <x> = ", sprintf("%.3f", round(mean(occs_df$occs_1), 3)), 
                             ",  <y> = ", sprintf("%.3f", round(mean(occs_df$occs_2), 3)))

      if(full_info) {
        plot_infos <- plot_infos_1
      } else {
        plot_infos <- plot_infos_2
      }
      
      if(!no_title) {
        p <- p + labs(title = this_title, subtitle = plot_infos)
      } else {
        p <- p + labs(subtitle = plot_infos)
      }
      
      if(use_pdf) {
        ggsave(plot = p, filename = paste0(plot_folder, "scatterplot_", method, "_" , sample_2, "_vs_", sample_1, "_min_cov_", min_coverage, "_old.pdf"),
               width = 6, height = 6, units = "cm")
      } else {
        ggsave(plot = p, filename = paste0(plot_folder, "scatterplot_", method, "_" , sample_2, "_vs_", sample_1, "_min_cov_", min_coverage, "_old.png"),
               width = 6, height = 6, units = "cm")
      }
    
    }
    
  }
  
  comparison_table <- data.frame(sample_1=sample_1, sample_2=sample_2, mean_occ_1=mean(occs_1), mean_occ_2=mean(occs_2), radius = radius, min_coverage = min_coverage,
                                 delta=delta, delta_2=delta_2, beta = beta, beta_2 = beta_2, alpha = alpha, alpha_2 = alpha_2, 
                                 r_value=r_value, r_value_rescaled=r_value_rescaled, num_comparisons=num_comparisons, d_value = d_value,
                                 stringsAsFactors = FALSE)
  return(comparison_table)
}


compare_two_occs_dfs <- function(method, sample_1, df_1, sample_2, df_2, min_coverage = NA, plot_folder = NA, 
                                 no_title = FALSE, my_alpha = NA, use_pdf = FALSE, full_info = !no_title) {
  
  # Sample properties:
  
  # RE sites are at least 200bp apart from each other. For this comparison we define a RE site area as the 10bp where the cut happens between 5bp and 6bp.
  
  # NP sites (called by nanopolish) are at least 11bp apart, as methylation sites that are closer are grouped together to one site. 
  # For nanopolish the sites on the + and - strand in the "CG" motif are always grouped together. 
  # NP dfs should have a start and end column such that end - start is at least 10bp
  
  # BS sites are arbitrary close to each other ("CG" already has two neighbouring sites with distance 1, the C on the + strand and the C on the - strand)
  
  # methods:
  
  # RE_vs_BS_5: BS sites within a fixed RE site are averaged (i.e. 10bp bins at RE cut positions)
  # RE_vs_RE_5: RE_1 and RE_2 sites are compared if their area overlap is at least 5 (i.e. 10bp bins at RE_1 cut positions)
  # BS_vs_BS_0: Compare identical sites 
  # BS_vs_BS_5: Compare average in bins of 10bp (genomic bins)
  
  # RE_vs_NP_5: RE and NP sites if the RE site touches the NP site
  # NP_vs_BS_5: Average all BS sites within the NP site (i.e. bins at NP sites with boundarys start and end)
  # NP_vs_NP_0: Compare identical sites (since so far only CpG NP data is available)
  
  # All sites need to have minimal coverage of min_coverage if count data is provided
  
  colnames(df_1)[grepl("mean_occ", colnames(df_1))] <- "occ"  # RE maps have a mean_occ column instead of an occ column
  colnames(df_2)[grepl("mean_occ", colnames(df_2))] <- "occ"
  
  colnames(df_1)[grepl("mean_eff_coverage", colnames(df_1))] <- "counts"  # RE maps have a mean_eff_coverage column instead of an counts column
  colnames(df_2)[grepl("mean_eff_coverage", colnames(df_2))] <- "counts"
  
  if(is.finite(min_coverage)) {
    if("counts" %in% colnames(df_1)) {
      cat("Removed", round(100*sum(df_1$counts < min_coverage)/nrow(df_1)), "% low coverage sites in sample 1.\n")
      df_1 <- subset(df_1, counts >= min_coverage)
    } else {
      cat("No counts column in df_1, min_coverage not checked.\n")
    }
    if("counts" %in% colnames(df_2)) {
      cat("Removed", round(100*sum(df_2$counts < min_coverage)/nrow(df_2)), "% low coverage sites in sample 2.\n")
      df_2 <- subset(df_2, counts >= min_coverage)
    } else {
      cat("No counts column in df_2, min_coverage not checked.\n")
    }
  }
  
  radius <- as.numeric(sub("[A-Za-z_]*", "", method))
  method_str <- sub("_[0-9]*$", "", method)
  
  if(grepl("(RE_vs_NP|NP_vs_BS)", method) & radius != 5) {  # NP methods have fixed radius given by the NP calling resolution
    stop("NP method with radius != 5 not implemented")
  }
  
  if(grepl("^RE", method) | method == "BS_vs_BS_0" | method == "NP_vs_BS_5" | method == "NP_vs_NP_0") {
    if(method_str == "RE_vs_BS") {
      bin <- -radius:(radius-1)
    } else if(method_str == "RE_vs_RE" | method == "BS_vs_BS_0" | method == "NP_vs_NP_0") {
      bin <- -radius:radius
    } else if (method == "NP_vs_BS_5") {
    } else {
      print(method)
      stop("Invalid method!")
    }
    
    pos_valid_list <- list()
    occs_1_list <- list()
    occs_2_bin_mean_list <- list()
    for(this_chr in unique(df_1$chr)) {
      
      df_1_chr <- subset(df_1, chr == this_chr & !is.na(occ))
      df_2_chr <- subset(df_2, chr == this_chr & !is.na(occ))
      
      pos_1_good <- df_1_chr$pos
      pos_2_good <- df_2_chr$pos
      
      occs_1 <- df_1_chr$occ
      
      if(method != "NP_vs_BS_5") {
        occs_2_vec <- rep(NA, max(c(pos_1_good + radius, pos_2_good)))
        occs_2_vec[pos_2_good] <- df_2_chr$occ
        
        occs_2_bin_mean <- sapply(pos_1_good, function(x) mean(occs_2_vec[x+bin], na.rm = TRUE))
      } else {
        occs_2_vec <- rep(NA, max(c(df_1_chr$end, pos_2_good)))
        occs_2_vec[pos_2_good] <- df_2_chr$occ
        
        occs_2_bin_mean <- mapply(function(s, e) mean(occs_2_vec[s:e], na.rm = TRUE), df_1_chr$start, df_1_chr$end)
      }
      
      pos_valid <- pos_1_good[is.finite(occs_2_bin_mean)]
      occs_2_bin_mean <- occs_2_bin_mean[is.finite(occs_2_bin_mean)]
      
      pos_valid_list[[this_chr]] <- pos_valid
      occs_1_list[[this_chr]] <- subset(df_1_chr, pos %in% pos_valid)$occ
      occs_2_bin_mean_list[[this_chr]] <- occs_2_bin_mean
    }
    occs_1 <- unlist(occs_1_list)
    occs_2 <- unlist(occs_2_bin_mean_list)
    
  } else if(grepl("^BS_vs_BS", method) & radius > 0) {  # bin BS sites of all samples (e.g radius=5: pos_binned = 10 <==> pos = 5,6,...,13,14)
    
    pos_binning <- 2*radius
    
    calc_binned_occs_df <- function(df) {
      df$pos_binned <- round(df$pos/pos_binning) * pos_binning
      df$chr_pos <- with(df, paste0(chr, "_", pos_binned))
      
      df_binned <- aggregate(df$occ, list(df$chr_pos), mean, na.rm = TRUE)
      colnames(df_binned) <- c("chr_pos", "occ")
      df_binned <- df_binned[!is.na(df_binned$occ), ]
      df_binned <- df_binned[order(df_binned$chr_pos), ]
      return(df_binned)
    }
    
    df_1 <- calc_binned_occs_df(df_1)
    df_2 <- calc_binned_occs_df(df_2)
    
    valid_chr_pos <- intersect(df_1$chr_pos, df_2$chr_pos)
    
    occs_1 <- with(df_1, occ[chr_pos %in% valid_chr_pos])
    occs_2 <- with(df_2, occ[chr_pos %in% valid_chr_pos])
    
  } else {
    print(method)
    stop("Method not valid or to be implemented!")
  }
  
  delta = mean(abs(occs_1 - occs_2))
  delta_2 = mean((occs_1 - occs_2)^2)
  
  r_value = cor(occs_1, occs_2)
  d_value = delta + abs(mean(occs_1) - mean(occs_2))
  
  num_comparisons <- length(occs_1)
  
  # plotting
  
  if(!is.na(plot_folder)) {
    
    if(is.na(my_alpha)) my_alpha <- min(0.3, 2500/num_comparisons)
    
    my_alpha <- max(0.01, my_alpha)
    
    this_title <- paste0(sample_1, " vs ", sample_2)
    
    occs_df <- data.frame(occs_1 = occs_1, occs_2 = occs_2)
    
    geom_text_size <- theme_paper[[1]]$axis.text$size
    
    occ_limits <- c(0, 1)
    p <- ggplot(occs_df, aes(x=occs_1, y=occs_2)) + 
      geom_abline(slope = 1, intercept = 0, linetype = "dashed",  color = "blue", size = 0.4) +
      geom_point(alpha = my_alpha, size = 0.4, shape = 16, stroke = 0) + 
      scale_x_continuous(labels = percent_no_sign, limits = (occ_limits)) + scale_y_continuous(labels = percent_no_sign, limits = (occ_limits)) +
      xlab(paste0("abs. occupancy ", sample_1, " / %")) + ylab(paste0("abs. occupancy ", sample_2, " / %")) + theme_paper
    
    plot_infos_1 <- paste0("<|x-y|> = ", sprintf("%.3f", round(delta, 3)), 
                           ",  <x> = ", sprintf("%.3f", round(mean(occs_df$occs_1), 3)), 
                           ",  <y> = ", sprintf("%.3f", round(mean(occs_df$occs_2), 3)), 
                           ",  d = ", sprintf("%.3f", round(d_value, 3)),
                           ",\nradius = ", radius, ",  ", num_comparisons, " sites,  r = ", sprintf("%.3f", round(r_value, 3)))
    
    plot_infos_2 <- paste0(num_comparisons, " sites,  <|x-y|> = ", sprintf("%.3f", round(delta, 3)), 
                           ",  <x> = ", sprintf("%.3f", round(mean(occs_df$occs_1), 3)), 
                           ",  <y> = ", sprintf("%.3f", round(mean(occs_df$occs_2), 3)))
    
    if(full_info) {
      plot_infos <- plot_infos_1
    } else {
      plot_infos <- plot_infos_2
    }
    
    if(!no_title) {
      p <- p + labs(title = this_title, subtitle = plot_infos)
    } else {
      p <- p + labs(subtitle = plot_infos)
    }
    
    if(use_pdf) {
      ggsave(plot = p, filename = paste0(plot_folder, "scatterplot_", method, "_" , sample_1, "_vs_", sample_2, "_min_cov_", min_coverage, ".pdf"),
             width = 6, height = 6, units = "cm")
    } else {
      ggsave(plot = p, filename = paste0(plot_folder, "scatterplot_", method, "_" , sample_1, "_vs_", sample_2, "_min_cov_", min_coverage, ".png"),
             width = 6, height = 6, units = "cm", type = "cairo")
    }
  }
  
  comparison_table <- data.frame(method=method, sample_1=sample_1, sample_2=sample_2, mean_occ_1=mean(occs_1), mean_occ_2=mean(occs_2), radius = radius, min_coverage = min_coverage,
                                 delta=delta, delta_2=delta_2, r_value=r_value, d_value = d_value, num_comparisons=num_comparisons,
                                 stringsAsFactors = FALSE)
  return(comparison_table)
}


add_ <- function(char_vec, sep = "_") {
  return(ifelse(char_vec == "", "", paste0(sep, char_vec)))
}


rescale_in_zero_to_one <- function(x, gamma) return(x / (x + (1-x) * gamma))


optimize_gamma <- function(vec, mean_goal) {  # entries in vec in [0, 1], i.e. occ or acc values; gamma for occ is 1/gamma for acc
  result <- optimize(function(gamma) sum(( mean(rescale_in_zero_to_one(vec, gamma)) - mean_goal ) ^ 2), lower = 0.01, upper = 100, maximum = FALSE)
  if(result$objective < 0.000001) {
    return(result$minimum) 
  } else {
    stop("Optimization not successful.")
  }
}


calc_Np1_positions <- function(dayd_df, max_tss_Np1_dist = 220) {
  start_path <- "../../"
  
  gene_data <- read.delim(paste0(start_path, "external_data/+1coordiantesETC_tirosh_32U.tab"))  # best results
  test <- gene_data[gene_data$name %in% gene_data$name[duplicated(gene_data$name)], ]
  
  gene_data$chr <- paste0("chr", as.numeric(as.roman(sub("chr", "", gene_data$chr))))
  gene_data_chr_list <- list()
  for(this_chr in unique(gene_data$chr)) {
    dyad_chr_df <- subset(dyad_df, chr == this_chr)
    for(s in c("+", "-")) {
      this_gene_data <- subset(gene_data, chr == this_chr & strand == s)
      this_gene_data$tss_Np1_dist <- sapply(this_gene_data$tss, function(tss_pos) calc_tss_Np1_distance(dyad_chr_df$aligned_pos, tss_pos, s))
      this_gene_data$Np1_pos <- with(this_gene_data, tss + ifelse(s=="+", 1, -1) * tss_Np1_dist)
      gene_data_chr_list[[paste0(this_chr, "_", s)]] <- this_gene_data
    }
  }
  gene_data <- bind_rows(gene_data_chr_list)
  gene_data <- gene_data[gene_data$tss_Np1_dist <= max_tss_Np1_dist, ]
  gene_data$aligned_pos <- gene_data$Np1_pos
  
  gene_data <- gene_data[!is.na(gene_data$aligned_pos), c("ID", "chr", "start", "end", "strand", "name", "aligned_pos")]
  
  return(gene_data)
}

calc_tss_Np1_distance <- function(vec, pos, s) {
  if(s == "+") {
    d <- vec - pos
  } else if(s == "-") {
    d <- pos - vec
  } else {
    stop("strand not valid!")
  }
  return(min(d[d>=0]))
}

load_scer_gene_data <- function(data_set = "BY0_CL_pm9") {
  start_path <- "../../"
  
  if(grepl("/nanopore/", getwd())) {  # for Elisa's NP runs
    start_path <- "../../../"
  }
  
  if(data_set == "transcripts_steinmetz") {  # TSS alignment
    gene_data <- read.delim(paste0(start_path, "external_data/occ_paper/transcripts_R64-1-1_steinmetz.txt"))
    gene_data$chr <- paste0("chr", as.numeric(as.roman(sub("chr", "", gene_data$chr))))
    gene_data$aligned_pos <- with(gene_data, ifelse(strand == "+", start, end))
  
  } else if(grepl("NPD", getwd())) {  # Daan's NP runs, needs NPD in working directory name for this to be TRUE
    gene_data <- read.delim(paste0("../../external_data/Yeast_annotation_Genome_Chereji_2018_Daan2.tsv"))
    colnames(gene_data) <- c("name", "chr", "start", "end", "strand", "TSS", "TTS", "plus_1_nucl", "minus_1_nucl",	"NDR_center",	"NDR_width")
    gene_data$aligned_pos <- gene_data$plus_1_nucl
    gene_data$ID <- gene_data$name  # some downstream functions use the ID instead of name
    gene_data$strand <- with(gene_data, ifelse(strand == 1, "+", "-"))
    gene_data$chr <- paste0("chr", as.numeric(as.roman(sub("chr", "", gene_data$chr))))
  
  } else {
    #gene_data <- read.delim(paste0(start_path, "external_data/ann_plus1_filtered.txt"))
    #gene_data$aligned_pos <- gene_data$plus1
    
    #gene_data <- read.delim(paste0(start_path, "external_data/+1coordiantesETC_DANPOS_32U.tab"))
    gene_data <- read.delim(paste0(start_path, "external_data/+1coordiantesETC_tirosh_32U.tab"))  # best results
    gene_data$BY0_all <- gene_data$plus1NK
    gene_data$plus1NK <- NULL
    
    gene_data$chr <- paste0("chr", as.numeric(as.roman(sub("chr", "", gene_data$chr))))
    
    Np1_files <- grep(".rda", dir(paste0(start_path, "external_data/occ_paper/called_dyads_and_Np1")), value = TRUE)
    for(f in Np1_files) {
      load(paste0(start_path, "external_data/occ_paper/called_dyads_and_Np1/", f))
      col <- sub(".rda", "", f)
      df <- centers[, c("ID", "plus1")]
      colnames(df)[2] <- col
      gene_data <- merge(gene_data, df[, c("ID", col)], all.x = TRUE, by = "ID")
    }
    
    # ggplot(gene_data, aes(x=BY0_CL_pm14-BY0_CL_pm9)) + geom_histogram(binwidth = 5) + xlim(c(-200, 200))
    # ggplot(gene_data, aes(x=BY0_CL_all-BY0_CL_pm9)) + geom_histogram(binwidth = 5) + xlim(c(-200, 200))
    # ggplot(gene_data, aes(x=BY0_CL_all-BY0_CL_pm14)) + geom_histogram(binwidth = 5) + xlim(c(-200, 200))
    # 
    # ggplot(gene_data, aes(x=BY0_all-BY0_CL_all)) + geom_histogram(binwidth = 5) + xlim(c(-200, 200))
    # ggplot(gene_data, aes(x=BY0_all-BY0_CL_pm9)) + geom_histogram(binwidth = 5) + xlim(c(-200, 200))
    # 
    # ggplot(gene_data, aes(x=Zhang-BY0_CL_all)) + geom_histogram(binwidth = 5) + xlim(c(-200, 200))
    # ggplot(gene_data, aes(x=Zhang-BY0_all)) + geom_histogram(binwidth = 5) + xlim(c(-200, 200))
    
    gene_data$aligned_pos <- gene_data[, data_set]
    gene_data$TSS <- gene_data$tss
    gene_data$tss <- NULL
    gene_data$TTS <- gene_data$tts
    gene_data$tts <- NULL
  }
  
  if(any(duplicated(gene_data$ID))) stop("Found duplicated gene IDs!")
  
  gene_data$name <- as.character(gene_data$name)
  
  if("TSS" %in% colnames(gene_data) & "TTS" %in% colnames(gene_data)) {
    gene_data <- gene_data[!is.na(gene_data$aligned_pos), c("ID", "chr", "start", "end", "strand", "name", "aligned_pos", "TSS", "TTS")]
  } else {
    gene_data <- gene_data[!is.na(gene_data$aligned_pos), c("ID", "chr", "start", "end", "strand", "name", "aligned_pos")]
  }
  gene_data <- gene_data[order(gene_data$chr, gene_data$start), ]
  rownames(gene_data) <- NULL
  
  return(gene_data)
}


load_dyad_positions <- function(data_set) {
  if(data_set == "BY0_CL_pm9") {
    dyad_df <- read.delim("../../external_data/occ_paper/called_dyads_and_Np1/BY0_CL_pm9.csv", sep = ";", header = TRUE)
    dyad_df$chr <- sub("chr0", "chr", dyad_df$chr)
    dyad_df <- subset(dyad_df, grepl("chr[0-9]*$", chr))
    dyad_df <- dyad_df[, c("chr", "smt_pos")]
    colnames(dyad_df) <- c("chr", "aligned_pos")
    dyad_df$name <- paste0("nucl", 1:nrow(dyad_df))
  } else if(data_set == "chereji_email") {
    dyad_df <- read.delim("../../external_data/occ_paper/called_dyads_and_Np1/chereji_email_my_header.bed", sep = "\t", header = TRUE)
    dyad_df$start <- dyad_df$start_bed + 1
    dyad_df$aligned_pos <- with(dyad_df, floor((start + end) / 2))
    dyad_df$chr <- paste0("chr", as.numeric(as.roman(sub("chr", "",  dyad_df$chr))))
    dyad_df <- dyad_df[, c("chr", "aligned_pos", "name")]
    dyad_df$name <- as.character(dyad_df$name)
  } else if(data_set == "zhang") {
    dyad_df <- read.delim("../../external_data/occ_paper/called_dyads_and_Np1/Zhang.csv", sep = "\t", header = TRUE)
    dyad_df$chr <- paste0("chr", as.numeric(as.roman(sub("chr", "",  dyad_df$chr))))
    dyad_df <- subset(dyad_df, grepl("chr[0-9]*$", chr))
    dyad_df <- dyad_df[, c("chr", "smt_pos")]
    colnames(dyad_df) <- c("chr", "aligned_pos")
    dyad_df$name <- paste0("nucl", 1:nrow(dyad_df))
  }
  dyad_df <- dyad_df[order(dyad_df$chr, dyad_df$aligned_pos), ]
  dyad_df$ID <- dyad_df$name
  dyad_df$strand <- "+"
  
  unique(dyad_df$chr)
  sapply(unique(dyad_df$chr), function(this_chr) max(subset(dyad_df, chr == this_chr)[, "aligned_pos"]))
  any(duplicated(dyad_df$name))
  
  return(dyad_df)
}

calc_aligned_positions <- function(occs_df, max_promoter_length = 1000, max_ORF_length = 1000, cut_at_gene_ends = FALSE, data_set = "BY0_CL_pm9", 
                                   dyads_instead_of_genes = FALSE) {
  
  occs_df <- occs_df[!is.na(occs_df$occ), grepl("^(chr|pos|occ|sd_occ|counts)$", colnames(occs_df))]  
  
  chrs <- sort(unique(as.character(occs_df$chr)))
  
  expected_labels <- sort(c(paste0("chr", 1:9), paste0("chr", 10:16)))  # NO LEADING ZERO
  if(length(chrs)==0 || !all(diag(sapply(expected_labels, grepl, chrs)))) {
    stop("Check chr labels in input occs_df! Not all 16 chromosomes were found!\n")
  }
  
  occs_df_chr_list <- split(occs_df, list(occs_df$chr))
  if(nrow(occs_df) > 1500000) {
    method <- "access_extended_df"
    occs_df_chr_list <- lapply(occs_df_chr_list, function(occs_df) merge(occs_df, data.frame(pos = 1:max(occs_df$pos)), by = "pos", all = TRUE))
  } else {
    method <- "subset_normal_df"
  }
  occs_df_chr_list <- lapply(occs_df_chr_list, function(occs_df) occs_df[order(occs_df$pos), ])
  rm(occs_df)
  
  if(!dyads_instead_of_genes) {
    gene_data <- load_scer_gene_data(data_set)
  } else {
    gene_data <- load_dyad_positions(data_set)
  }

  aligned_genes_list <- list()
  for (i in 1:nrow(gene_data)) {
    chr_temp <- gene_data$chr[i]
    aligned_pos <- gene_data[i, "aligned_pos"]
    if (gene_data[i, "strand"] == '+') {
      genome_pos_min <- aligned_pos - max_promoter_length
      if(cut_at_gene_ends) {
        genome_pos_max <- min(max_ORF_length + aligned_pos, gene_data[i, "end"])
      } else {
        genome_pos_max <- max_ORF_length + aligned_pos
      }
      if(method == "subset_normal_df") {
        new_gene <- subset(occs_df_chr_list[[chr_temp]], pos >= genome_pos_min & pos <= genome_pos_max)
      } else if(method == "access_extended_df") {
        new_gene <- occs_df_chr_list[[chr_temp]][genome_pos_min:genome_pos_max, ]
        new_gene <- new_gene[!is.na(new_gene$chr), ]
      } else {
        stop("invalid method!")
      }
      new_gene$pos = new_gene$pos - aligned_pos
    } else {
      if(cut_at_gene_ends) {
        genome_pos_min <- max(gene_data[i, "start"], aligned_pos - max_ORF_length)
      } else {
        genome_pos_min <- aligned_pos - max_ORF_length
      }
      genome_pos_max <- aligned_pos + max_promoter_length
      if(method == "subset_normal_df") {
        new_gene <- subset(occs_df_chr_list[[chr_temp]], pos >= genome_pos_min & pos <= genome_pos_max)
      } else if(method == "access_extended_df") {
        new_gene <- occs_df_chr_list[[chr_temp]][genome_pos_min:genome_pos_max, ]
        new_gene <- new_gene[!is.na(new_gene$chr), ]
      } else {
        stop("invalid method!")
      }
      new_gene$pos = aligned_pos - new_gene$pos
      new_gene = new_gene[nrow(new_gene):1,]  # revert pos order (not needed)
    }
    if(nrow(new_gene)>0) {
      new_gene$gene_number = i
      new_gene$gene_name = gene_data[i, "name"]
      new_gene$gene_id = gene_data[i, "ID"]
      new_gene$chr <- NULL
      aligned_genes_list[[i]] <- new_gene
    }
  }
  aligned_genes <- bind_rows(aligned_genes_list)  # makes one data.frame from the list
  
  aligned_genes
}


calc_average_gene_occ_df_from_genes_df <- function(genes_df) {
  calc_average_gene_occ_df(genes_df = genes_df)
}

calc_average_gene_occ_df <- function(occs_df, genes_df, pos_binning = 10, filter_size = 1, max_promoter_length = 1000, max_ORF_length = 1000, 
                                     cut_at_gene_ends = FALSE, normalize_cov_for_each_gene = FALSE, data_set = "BY0_CL_pm9", dyads_instead_of_genes = FALSE) {
  if(missing(genes_df) & !missing(occs_df)) {
    genes_df_given <- FALSE
    if(!("occ" %in% colnames(occs_df)) & ("cov" %in% colnames(occs_df))) {
      coverage_given <- TRUE
      occs_df$occ <- occs_df$cov  # treat coverage as occupancy
      occs_df$cov <- NULL
    } else {
      coverage_given <- FALSE
    }
    genes_df <- calc_aligned_positions(occs_df, max_promoter_length, max_ORF_length, cut_at_gene_ends, data_set, dyads_instead_of_genes = dyads_instead_of_genes)
   
  } else {
    genes_df_given <- TRUE
    if(!("occ" %in% colnames(genes_df)) & ("cov" %in% colnames(genes_df))) {
      coverage_given <- TRUE
      genes_df$occ <- genes_df$cov  # treat coverage as occupancy
      genes_df$cov <- NULL
    } else {
      coverage_given <- FALSE
    }
  }
  
  if(coverage_given) {
    if(normalize_cov_for_each_gene) {
      genes_df_list <- split(genes_df, list(genes_df$gene_id))
      genes_df_list <- lapply(genes_df_list, function(df) {
        df$occ <- df$occ / mean(df$occ)
        return(df)
      })
      genes_df <- bind_rows(genes_df_list)
    } else {  # normalize with the global mean
      genes_df$occ <- genes_df$occ / mean(genes_df$occ, na.rm = TRUE)
    }
  }
  
  genes_df$pos_binned = round(genes_df$pos/pos_binning, 0) * pos_binning
  av_gene_df = aggregate(genes_df[, c('occ')], list(genes_df$pos_binned), mean, na.rm = TRUE)
  colnames(av_gene_df) <- c("pos", "occ")
  
  av_gene_df_2 = aggregate(genes_df[, c('occ')], list(genes_df$pos_binned), length)
  colnames(av_gene_df_2) <- c("pos", "num_sites_in_bin")
  
  av_gene_df <- merge(av_gene_df, av_gene_df_2)
  
  if("counts" %in% colnames(genes_df)) {
    av_gene_df_3 = aggregate(genes_df[, c('counts')], list(genes_df$pos_binned), mean, na.rm = TRUE)
    colnames(av_gene_df_3) <- c("pos", "mean_counts_per_site")
    
    av_gene_df <- merge(av_gene_df, av_gene_df_3)
  }
  
  if("sd_occ" %in% colnames(genes_df)) {
    av_gene_df_4 = aggregate(genes_df[, c('sd_occ')], list(genes_df$pos_binned), mean, na.rm = TRUE)
    colnames(av_gene_df_4) <- c("pos", "mean_sd_occ")
    
    av_gene_df <- merge(av_gene_df, av_gene_df_4)
  }
  
  av_gene_df$occ_smoothed <- NA
  av_gene_df$num_sites_smoothed <- NA
  n <- nrow(av_gene_df)
  av_gene_df$occ_smoothed[ceiling(filter_size/2):(n-floor(filter_size/2))] <- rollmean(av_gene_df$occ, filter_size)
  av_gene_df$num_sites_smoothed[ceiling(filter_size/2):(n-floor(filter_size/2))] <- rollmean(av_gene_df$num_sites_in_bin, filter_size)
  
  if(!coverage_given & !dyads_instead_of_genes & data_set != "transcripts_steinmetz") {
  
    lag_min_max <- round(c(100/pos_binning, 250/pos_binning))
    
    autocor <- acf(av_gene_df$occ[av_gene_df$pos >= 0 & av_gene_df$pos <= min(1000, max_ORF_length)], na.action=na.pass, lag.max = lag_min_max[2], plot = FALSE)[["acf"]][lag_min_max[1]:lag_min_max[2]]
    autocor_max <- max(autocor)
    autocor_lag <- lag_min_max[1]-1 + mean(which(autocor == max(autocor)))
    
    # max min differences
    #period <- lag_min_max[1]-1 + which(autocor == max(autocor))[1]
    period <- round(180 / pos_binning)
    if(is.finite(period)) {
      maxima <- rep(NA, 3)
      minima <- rep(NA, 3)
      for(i in 1:3) {
        max_area <- which(av_gene_df$pos == 0) + period * (i-1) + seq(-round(period/2), round(period/2), 1)
        min_area <- which(av_gene_df$pos == 0) + period * (i-1) + seq(0, period, 1)
        maxima[i] <- max(av_gene_df$occ[max_area])
        minima[i] <- min(av_gene_df$occ[min_area])
      }
      mean_max_min_diff <- mean(maxima - minima)
    } else {
      mean_max_min_diff <- NA
    }
  } else {
    autocor_lag <- NA
    autocor_max <- NA
    mean_max_min_diff <- NA
  }
  if(coverage_given) {  # bin coverage to reduce memory needed
    temp_df <- aggregate(genes_df[, c('occ')], list(paste0(genes_df$gene_id, "__", genes_df$gene_name, "__", genes_df$pos_binned)), mean, na.rm = TRUE)
    names(temp_df) <- c("info", "cov")
    temp_info <- str_split(temp_df$info, pattern = "__")
    temp_df$gene_id <- sapply(temp_info, "[", 1)
    temp_df$gene_name <- sapply(temp_info, "[", 2)
    temp_df$pos <- as.integer(sapply(temp_info, "[", 3))
    temp_df$info <- NULL
    genes_df <- temp_df
    
    av_gene_df$cov <- av_gene_df$occ
    av_gene_df$occ <- NULL
  } else {
    genes_df$pos_binned <- NULL
  }
  
  if(!genes_df_given) {
    genes_df$gene_id <- as.character(genes_df$gene_id)
    return(list("av_gene_df" = av_gene_df, "genes_df" = genes_df, "info_df" = data.frame(pos_binning, autocor_lag, autocor_max, mean_max_min_diff)))
  } else {
    return(list("av_gene_df" = av_gene_df, "info_df" = data.frame(pos_binning, autocor_lag, autocor_max, mean_max_min_diff)))
  }
}


add_sample_names_to_aligned_dfs_list <- function(aligned_dfs_list, name_add = "") {
  
  saved_names <- names(aligned_dfs_list)
  aligned_dfs_list <- lapply(names(aligned_dfs_list), function(name) {
    aligned_dfs_list[[name]][[1]]$sample <- paste0(name_add, name)
    aligned_dfs_list[[name]][[2]]$sample <- paste0(name_add, name)
    aligned_dfs_list[[name]][[3]] <- cbind(data.frame(sample = paste0(name_add, name)), aligned_dfs_list[[name]][[3]])
    return(aligned_dfs_list[[name]])
  })
  names(aligned_dfs_list) <- saved_names
  return(aligned_dfs_list)
}


plot_sorted_genes <- function(genes_df, sorted_genes_df, list_name, map_name, plot_folder, pos_binnings = c(5, 10, 20), file_type = "png", 
                              window = NA, plot_quintiles_heatmap = TRUE) {
  
  occ_given <- "occ" %in% colnames(genes_df)
  if(!occ_given) {
    genes_df$occ <- genes_df$cov
  }
  
  if(missing(sorted_genes_df)) {
    if(all(!is.na(window))) {
      genes_df_to_sort <- genes_df[genes_df$pos %in% window, ]
    } else {
      genes_df_to_sort <- genes_df
    }
    sorted_genes_df <- aggregate(genes_df_to_sort$occ, list(genes_df_to_sort$gene_id), mean, na.rm = TRUE)
    colnames(sorted_genes_df) <- c("gene_id", "mean_signal")
    sorted_genes_df <- sorted_genes_df[order(sorted_genes_df$mean_signal), ]
    sorted_genes_df$sorted_number <- 1:nrow(sorted_genes_df)
    sorted_genes_df$sample <- map_name
    if(missing(list_name)) list_name <- paste0("self_sorted_from_", min(window), "_to_", max(window))
    
    min_cov <- 1.5*quantile(sorted_genes_df$mean_signal, 0.01, na.rm = TRUE)
    max_cov <- 1.5*quantile(sorted_genes_df$mean_signal, 0.99, na.rm = TRUE)
    
    if(!is.na(plot_folder)) {
      p <- ggplot(sorted_genes_df, aes(x=mean_signal)) + geom_histogram(bins = 40) + theme_paper +
        xlab(paste0(ifelse(occ_given, "mean abs. gene occupancy / %", "mean norm. gene coverage"), 
                    ifelse(all(!is.na(window)), paste0(" (", min(window), " to ", max(window), ")"), ""))) +
        ggtitle(map_name)
      if(occ_given) {
        p <- p + scale_x_continuous(labels = percent_no_sign, limit = c(0, 1), breaks = seq(0, 1, 0.25))
      } else {
        p <- p + scale_x_continuous(limit = c(ifelse(min_cov < 0, min_cov, 0), max_cov))
      }
      if(file_type == "png") {
        ggsave(plot = p, paste0(plot_folder, "signal_histogram_", list_name, "_using_", map_name, ".", file_type), width = 7, height = 5, units = "cm", type = "cairo")
      } else {
        ggsave(plot = p, paste0(plot_folder, "signal_histogram_", list_name, "_using_", map_name, ".", file_type), width = 7, height = 5, units = "cm")
      }
    }
  }
  
  if(plot_quintiles_heatmap & !is.na(plot_folder)) {
  
    genes_df <- merge(genes_df, sorted_genes_df)
  
    for(pos_binning in pos_binnings) {
      genes_df$pos_binned <- round(genes_df$pos/pos_binning) * pos_binning
      binned_df <- unique(genes_df[, c("gene_id", "sorted_number", "pos_binned")])
      binned_df$ID <- with(binned_df, paste0(gene_id, "_", pos_binned))
      temp_df <- aggregate(genes_df$occ, list(with(genes_df, paste0(gene_id, "_", pos_binned))), mean)
      colnames(temp_df) <- c("ID", "occ")
      binned_df <- merge(binned_df, temp_df, by = "ID")
      rm(temp_df)
      max_cov <- quantile(binned_df$occ, probs = 0.95, na.rm = TRUE)
  
      ggplot(binned_df, aes(x=pos_binned, y=sorted_number, fill=occ)) + geom_raster() +
        scale_fill_viridis(direction = -1, limits = c(0, ifelse(occ_given, 1, max_cov))) +
        theme_paper + labs(fill=ifelse(occ_given, "abs. occupancy", "norm. coverage"))
      if(file_type == "png") {
        ggsave(paste0(plot_folder, "heatmap_", list_name, "_using_", map_name, "_pos_binning_", pos_binning, ".", file_type), width = 15, height = 25, units = "cm", type = "cairo")
      } else {
        ggsave(paste0(plot_folder, "heatmap_", list_name, "_using_", map_name, "_pos_binning_", pos_binning, ".", file_type), width = 15, height = 25, units = "cm")
      }
    }
  
    genes_df$quintile <- with(genes_df, ceiling(5 * sorted_number / max(sorted_number)))
    genes_df_list <- split(genes_df, list(genes_df$quintile))
    temp_list <- lapply(genes_df_list, calc_average_gene_occ_df_from_genes_df)
    quintiles_df <- bind_rows(lapply(names(temp_list), function(name) {
      df <- temp_list[[name]][[1]]
      df$quintile <- name
      return(df)
    }))
    rm(temp_list)
    max_cov <- max(quintiles_df$occ, na.rm = TRUE)
  
    p <- ggplot(quintiles_df, aes(x=pos, y=occ, color=quintile)) +
      geom_line(size = 0.35, alpha = 0.8) + scale_color_viridis(direction = -1, discrete = TRUE, end = 0.9) + theme_paper +
      ylab(ifelse(occ_given, "gene averaged abs. occupancy / %", "gene averaged norm. coverage")) + xlab("distance to in vivo +1 nucleosome / bp")
    if(occ_given) {
      p <- p + scale_y_continuous(labels = percent_no_sign, limit = c(0, 1), breaks = seq(0, 1, 0.25))
    } else {
      p <- p + scale_y_continuous(limit = c(0, max_cov), breaks = seq(0, max_cov, 0.5))
    }
    if(file_type == "png") {
      ggsave(plot = p, paste0(plot_folder, "averaged_quintiles_", list_name, "_using_", map_name, ".", file_type), width = 7, height = 5, units = "cm", type = "cairo")
    } else {
      ggsave(plot = p, paste0(plot_folder, "averaged_quintiles_", list_name, "_using_", map_name, ".", file_type), width = 7, height = 5, units = "cm")
    }
  }
  
  return(sorted_genes_df)
}

