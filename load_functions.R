load_read_infos <- function(fastq_read_info_file, bam_read_info_file, barcode_file) {
  fastq_read_infos <- read.table(fastq_read_info_file, sep = '\t', header = TRUE, stringsAsFactors = FALSE)
  bam_read_infos <- read.table(bam_read_info_file, sep = '\t', header = TRUE, stringsAsFactors = FALSE)
  read_infos <- merge(bam_read_infos, fastq_read_infos)

  read_infos$start_time <- strptime(read_infos$start_time, format='%Y-%m-%dT%H:%M:%SZ', tz='UTC')
  read_infos$start_time <- as.numeric(difftime(read_infos$start_time, min(read_infos$start_time), units = "hours"))
  
  read_infos$ref_seq_detailed <- read_infos$ref_seq
  read_infos$ref_seq <- sub(".*\\|", "", read_infos$ref_seq_detailed)
  read_infos$ref_seq_detailed <- as.factor(read_infos$ref_seq_detailed)
  read_infos$ref_seq <- as.factor(read_infos$ref_seq)
  
  if(all(is.na(read_infos$barcode))) {
    warning("All barcode information is NA. Assuming no barcoding and setting the barcode to 00 (just to have it defined)")
    cat("\n")
    read_infos$barcode <- "barcode00"
  }
  
  read_infos$barcode <- as.factor(sub("barcode", "", read_infos$barcode))
  read_infos$strand <- factor(read_infos$strand, levels=c('+', '-'))
  
  read_infos$sample <- NA
  read_infos$methylation_time <- NA
  
  read_infos$log10_aligned_reference_length <- log10(read_infos$aligned_reference_length)
  
  if (!missing(barcode_file)) {
    barcode_infos <- read.table(barcode_file, sep = '\t', header = TRUE, stringsAsFactors = FALSE)
    barcodes <- sub("barcode", "", barcode_infos$barcode)

    for(barcode in barcodes) {
      reads_temp <- read_infos$barcode == barcode
      read_infos[reads_temp, "sample"] <- barcode_infos[barcodes == barcode, "sample"]
      if("methylation_time" %in% colnames(barcode_infos)) {
        read_infos[reads_temp, "methylation_time"] <- barcode_infos[barcodes == barcode, "methylation_time"]
      }
    }
    if("methylation_time" %in% colnames(barcode_infos)) {
      read_infos$methylation_time <- factor(read_infos$methylation_time, unique(barcode_infos$methylation_time))
    }
    read_infos$sample <- factor(read_infos$sample, barcode_infos$sample)
  }
  
  read_infos <- subset(read_infos, mapping_quality >= 30)
  
  na_samples <- is.na(read_infos$sample)
  if(any(na_samples)) {
    warning(paste0("Ignoring ", length(na_samples[na_samples]), " reads with unspecified barcode/sample.\n"))
  }
  
  read_infos$read_identifier <- with(read_infos, paste0(read_name, "_", ref_seq_detailed, "_", reference_start, "_", aligned_reference_length))
                             
  read_infos <- read_infos[!na_samples, ]
  read_infos
}


load_read_stats <-function(fastq_read_info_file, bam_read_info_file, read_stats_file, barcode_file) {
  read_infos <- load_read_infos(fastq_read_info_file, bam_read_info_file, barcode_file)
  read_stats <- read.table(read_stats_file, sep = '\t', header = TRUE, stringsAsFactors = FALSE)

  # Shifting start and end positions by +1 to give correct methylation cytosine position
  read_stats$first_site <- read_stats$first_site + 1
  read_stats$last_site <- read_stats$last_site + 1
  
  read_stats$dist_first_to_last_site <- read_stats$last_site - read_stats$first_site
  
  read_stats$site_fraction_methylated <- read_stats$sites_methylated / read_stats$sites_total
  read_stats$site_call_fraction_methylated <- read_stats$sites_methylated / (read_stats$sites_total - read_stats$sites_ambiguous)
  read_stats$site_fraction_ambiguous <- read_stats$sites_ambiguous / read_stats$sites_total
  
  read_stats$group_fraction_methylated <- read_stats$groups_methylated / read_stats$groups_total
  read_stats$group_call_fraction_methylated <- read_stats$groups_methylated / (read_stats$groups_total - read_stats$groups_ambiguous)
  read_stats$group_fraction_ambiguous <- read_stats$groups_ambiguous / read_stats$groups_total
  
  read_stats <- merge(read_stats[, -2], read_infos)  # exclude chromosome column (ref_seq already in read_infos)
  read_stats$barcode <- as.character(read_stats$barcode)
  
  read_stats
}


load_barcoded_read_stats <-function(fastq_read_info_file, bam_read_info_file, read_stats_file, barcode_file) {
  barcode_infos <- read.table(barcode_file, sep = '\t', header = TRUE, stringsAsFactors = FALSE)
  barcode_folders <- barcode_infos$barcode
  
  read_stats_list  <- vector("list", length = length(barcode_folders))
  for(i in 1:length(barcode_folders)) {
    barcode_folder <- barcode_folders[i]
    fastq_read_info_file_temp = paste0("../", barcode_folder, "/", fastq_read_info_file)
    bam_read_info_file_temp = paste0("../", barcode_folder, "/", bam_read_info_file)
    read_stats_file_temp =  paste0("../", barcode_folder, "/", read_stats_file)
    read_stats_list[[i]] <- load_read_stats(fastq_read_info_file = fastq_read_info_file_temp, bam_read_info_file = bam_read_info_file_temp, read_stats_file = read_stats_file_temp, barcode_file = barcode_file)
    read_stats_list[[i]]$ref_seq_detailed <- as.character(read_stats_list[[i]]$ref_seq_detailed)
    read_stats_list[[i]]$ref_seq <- as.character(read_stats_list[[i]]$ref_seq)
  }
  read_stats <- bind_rows(read_stats_list)
  read_infos$ref_seq_detailed <- as.factor(read_infos$ref_seq_detailed)
  read_infos$ref_seq <- as.factor(read_infos$ref_seq)

  read_stats
}


load_site_stats <- function(site_stats_file) {
  site_stats <- read.table(file = site_stats_file, sep = '\t', header = TRUE, stringsAsFactors = FALSE);

  names(site_stats)[1] <- "ref_seq_detailed"
  site_stats$ref_seq <- sub(".*\\|", "", site_stats$ref_seq_detailed)
  site_stats$ref_seq_detailed <- as.factor(site_stats$ref_seq_detailed)
  site_stats$ref_seq <- as.factor(site_stats$ref_seq)
  
  site_stats$hits_called <- site_stats$hits_total - site_stats$hits_ambiguous
  site_stats$hit_fraction_methylated <- site_stats$hits_methylated / site_stats$hits_total
  site_stats$hit_call_fraction_methylated <- site_stats$hits_methylated / site_stats$hits_called
  site_stats$hit_fraction_ambiguous <- site_stats$hits_ambiguous / site_stats$hits_total
  
  # Shifting start and end positions by +1 to give correct methylation cytosine position
  site_stats$start <- site_stats$start + 1
  site_stats$end <- site_stats$end + 1
  
  site_stats$position <- (site_stats$start + site_stats$end) / 2
  
  site_stats
}


load_barcoded_site_stats <- function(site_stats_file, barcode_file, project_folder = "../") {
  barcode_infos <- read.table(paste0(project_folder, barcode_file), sep = '\t', header = TRUE, stringsAsFactors = FALSE)
  barcode_folders <- barcode_infos$barcode

  site_stats_list <- vector(mode = "list", length = length(barcode_folders))
  for(i in 1:length(barcode_folders)) {
    barcode_folder <- barcode_folders[i]
    site_stats <- load_site_stats(site_stats_file = paste0(project_folder, barcode_folder, "/", site_stats_file))
    site_stats$ref_seq_detailed <- as.character(site_stats$ref_seq_detailed)
    site_stats$ref_seq <- as.character(site_stats$ref_seq)
    site_stats$barcode <- sub("barcode", "", barcode_folder)
    site_stats$sample <- barcode_infos$sample[i]
    if("methylation_time" %in% colnames(barcode_infos)) {
      site_stats$methylation_time <- barcode_infos$methylation_time[i]
    } else {
      site_stats$methylation_time <- NA
    }
    site_stats_list[[i]] <- site_stats
  }
  
  site_stats_combined <- bind_rows(site_stats_list) 
  site_stats_combined$barcode <- factor(site_stats_combined$barcode, sub("barcode", "", barcode_folders))
  site_stats_combined$sample <- factor(site_stats_combined$sample, barcode_infos$sample)
  if("methylation_time" %in% colnames(barcode_infos)) {
    site_stats_combined$methylation_time  <-  factor(site_stats_combined$methylation_time, unique(barcode_infos$methylation_time))
  }
  site_stats_combined$ref_seq_detailed <- as.factor(site_stats_combined$ref_seq_detailed)
  site_stats_combined$ref_seq <- as.factor(site_stats_combined$ref_seq)
  
  site_stats_combined
}


calc_and_save_meth_data <- function(meth_file, barcode_file, log_lik_ratio_threshold = 2.5, remove_columns = c("log_lik_unmethylated", "num_calling_strands")) {
  barcode_infos <- read.table(barcode_file, sep = '\t', header = TRUE, stringsAsFactors = FALSE)
  barcode_folders <- barcode_infos$barcode

  warning("Shifting start and end positions by +1 to give correct methylation cytosine position.\n")
  
  for(i in 1:length(barcode_folders)) {
    barcode_folder <- barcode_folders[i]
    meth_data <- read.table(paste0("../", barcode_folder, "/", meth_file), sep = '\t', header = TRUE, stringsAsFactors = FALSE)
    
    names(meth_data)[grepl("chromosome", names(meth_data))] <- "ref_seq"
    names(meth_data)[grepl("(num_motifs|num_cpgs)", names(meth_data))] <- "group_size"  # num_cpgs was renamed to num_motifs in the output of Nanopolish
    if(!all(c("ref_seq", "group_size") %in% colnames(meth_data))) {
      stop("Column renaming not successful!")
    }
    
    meth_data[, remove_columns] <- NULL
    
    meth_data$ref_seq_detailed <- meth_data$ref_seq
    meth_data$ref_seq <- sub(".*\\|", "", meth_data$ref_seq_detailed)
    meth_data$ref_seq_detailed <- NULL  # don't track ref_seq_detailed in meth_data anymore
    
    if(nrow(meth_data) > 0) {
      meth_data$sample <- barcode_infos$sample[i]
      
      meth_data$start <- as.integer(meth_data$start + 1)  # as.integer saves memory
      meth_data$end <- as.integer(meth_data$end + 1)
      meth_data$position <- (meth_data$start + meth_data$end)/2
      meth_data$hit_state <- as.integer(1*(meth_data$log_lik_ratio >= log_lik_ratio_threshold) - 1*(meth_data$log_lik_ratio <= -log_lik_ratio_threshold))  # 1: methylated, -1: not methylated, 0: ambiguous

      dir.create(paste0("../", barcode_folder, "/methylation/meth_data"), showWarnings = FALSE)
      lapply(split(meth_data, meth_data$ref_seq), function(meth_data) save(meth_data, file = paste0("../", barcode_folder, "/methylation/meth_data/", meth_data$ref_seq[1], ".RData")))
    } else {
      warning(paste0("No meth_data saved for ", barcode_folder, "\n"))
    }
  }
}


load_meth_data <- function(barcode_subset, ref_seq_subset, read_name_subset, remove_columns = c("log_lik_methylated", "sequence")) {
  
  # Include not needed columns in remove_columns to decrease the memory footprint of loaded meth_data!
  
  meth_data_list <- list()
  
  for(barcode in barcode_subset) {
    for(ref_seq in ref_seq_subset) {
      meth_data_file <- paste0("../", barcode, "/methylation/meth_data/", ref_seq, ".RData")
      if(file.exists(meth_data_file)) {
        load(file = meth_data_file)
        meth_data[, remove_columns] <- NULL
        if(!missing(read_name_subset)) {
          meth_data <- subset(meth_data, read_name %in% read_name_subset)
        }
        meth_data_list[[paste0(barcode, "_", ref_seq)]] = meth_data
      } else {
        warning(paste0("meth_data for ", barcode, " comined with ", ref_seq, " not found!\n"))
      }
    }
  }
  
  meth_data <- bind_rows(meth_data_list)
  
  # using factors makes the memory footprint smaller (tested with object_size()) but when saved to disk with save() the file becomes larger
  meth_data$sample <- factor(meth_data$sample, barcode_infos$sample)
  meth_data$ref_seq <- as.factor(meth_data$ref_seq)
  meth_data$read_name <- as.factor(meth_data$read_name)
  
  return(meth_data)
}


load_meth_data_chr <- function(read_name_subset, remove_columns = c("log_lik_methylated", "sequence")) {
  return(load_meth_data(barcode_subset = barcode_infos$barcode, ref_seq_subset = chr_ref_seqs, read_name_subset = read_name_subset, remove_columns = remove_columns))
}
