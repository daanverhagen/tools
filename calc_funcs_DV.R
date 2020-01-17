# Daan Verhagen
# 20191023
# calculation functions

#### site_stats calculations for specific loci

### NOTE: This is an adaptation to the original calc_site_stats script only to be used in combination with aligned reads of a specific loci (ie GRF, N1, ORI)

calc_site_stats_loci <- function(read_names, meth_data) {
  meth_data <- subset(meth_data, read_name %in% read_names)
  meth_data$position <- (meth_data$start + meth_data$end)/2
  meth_data <- meth_data[, c("ref_seq", "start", "end", "group_size", "sample", "hit_state", "position")]
  # meth_data$site_ID <- paste(meth_data$sample, meth_data$ref_seq, meth_data$start, sep = "_")  # results in only hits with the same site_ID will be aggregated
  meth_data$hits_total <- 1
  meth_data$hits_methylated <- as.integer(meth_data$hit_state == 1)
  meth_data$hits_ambiguous <- as.integer(meth_data$hit_state == 0)
  
  site_stats <- aggregate(meth_data[, c("hits_total", "hits_methylated", "hits_ambiguous")], list(meth_data$position), sum) # aggregates the data based on position
  names(site_stats)[1] <- "position"
  
  site_stats$hits_called <- site_stats$hits_total - site_stats$hits_ambiguous
  site_stats$hit_fraction_methylated <- site_stats$hits_methylated / site_stats$hits_total
  site_stats$hit_call_fraction_methylated <- site_stats$hits_methylated / site_stats$hits_called
  site_stats$hit_fraction_ambiguous <- site_stats$hits_ambiguous / site_stats$hits_total
  
  # site_stats$position <- (site_stats$start + site_stats$end) / 2
  
  site_stats_loci <<- site_stats # writes the data to the global environment 
}

site_stats_loci_plot <- function(meth_data, binsize = 10, sample) {
  # make vector of read_names to be used later when calculating site_stats
  aligned_read_names <- meth_data$read_name
  # Add sample information to data.frame
  meth_data$sample <- sample
  # requires custom script in calc_functions_DV.R
  calc_site_stats_loci(read_names = aligned_read_names, meth_data = meth_data)
  
  
  # bin the reads for a smoother line
  binsize = binsize
  site_stats_loci$start_binned <- round(site_stats_loci$position/binsize, 0)*binsize
  print(paste0("Percentage of sites removed due to NaN: ", 100*round(sum(site_stats_loci$hits_methylated == 0 & site_stats_loci$hits_called == 0)/nrow(site_stats_loci), 4),"%")) # percentage of reads that will be removed
  site_stats_loci <-  site_stats_loci[!site_stats_loci$hit_call_fraction_methylated == "NaN",] # remove sites with only ambiguous hits
  ag_site_stats_loci <<- aggregate(site_stats_loci, by = list(site_stats_loci$start_binned), FUN = mean, na.action = na.pass)
  paste0("average occupancy around GRFs: ", round(1-mean(site_stats_loci$hit_call_fraction_methylated), 2))
  print("DONE. Use ag_site_stats_loci for plotting")
}  



# Isolate and calculate a region of the rDNA locus

rDNA_meth_calc <-  function(barcode, ref_seq, rDNA_start, rDNA_end, meth_percentage = 0.7) {
  # to calculate methylation of a certain region you need meth_data_chr
  meth_data_chr <-  load_meth_data(barcode_subset = paste0(barcode), ref_seq_subset = paste0(ref_seq))
  # subset sites below rDNA_start and above rDNA_end
  meth_data_chr_rDNA <- meth_data_chr[!(meth_data_chr$start <= rDNA_start) & !(meth_data_chr$end >= rDNA_end),]
  # find unique reads in dataframe
  rDNA_reads <- unique(meth_data_chr_rDNA$read_name)
  
  # add read_length of read_name to each row but makes it significantly slower
  # final_list <- data.frame()
  # for(readname in rDNA_reads) {
  #   temp_meth_data_chr_rDNA <- subset(meth_data_chr_rDNA, read_name == readname)
  #   temp_meth_data_chr_rDNA$read_length <- max(temp_meth_data_chr_rDNA$position) - min(temp_meth_data_chr_rDNA$position)
  #   final_list <- rbind(final_list, temp_meth_data_chr_rDNA)
  # }
  # meth_data_chr_rDNA <- final_list
  
  meth_data_chr_rDNA$position <- (meth_data_chr_rDNA$start + meth_data_chr_rDNA$end)/2
  meth_data_chr_rDNA <- meth_data_chr_rDNA[, c("ref_seq", "start", "end", "group_size", "sample", "hit_state", "position", "read_name")]
  meth_data_chr_rDNA$hits_total <- 1
  meth_data_chr_rDNA$hits_methylated <- as.integer(meth_data_chr_rDNA$hit_state == 1)
  meth_data_chr_rDNA$hits_ambiguous <- as.integer(meth_data_chr_rDNA$hit_state == 0)
  
  site_stats_rDNA <- aggregate(meth_data_chr_rDNA[, c("hits_methylated", "hits_total", "hits_ambiguous")], list(meth_data_chr_rDNA$read_name), mean) # aggregates the data based on read name
  names(site_stats_rDNA)[1] <- "read_name"
  
  site_stats_rDNA$hits_called <- site_stats_rDNA$hits_total - site_stats_rDNA$hits_ambiguous
  site_stats_rDNA$hit_fraction_methylated <- site_stats_rDNA$hits_methylated / site_stats_rDNA$hits_total
  site_stats_rDNA$hit_call_fraction_methylated <- site_stats_rDNA$hits_methylated / site_stats_rDNA$hits_called
  site_stats_rDNA$hit_fraction_ambiguous <- site_stats_rDNA$hits_ambiguous / site_stats_rDNA$hits_total
  
  length(unique(site_stats_rDNA$read_name))
  
  # write to workspace/global env
  site_stats_rDNA <<- site_stats_rDNA
  
  print(
    paste0("Calculation window from - to: ", rDNA_start, " - ", rDNA_end)
  )
  print(
    paste0("Percentage MORE than ", (100*meth_percentage), "% methylated: ", (100*round(nrow(subset(site_stats_rDNA, hit_call_fraction_methylated >= meth_percentage))/nrow(site_stats_rDNA),4)),"%")
  )
  print(
    paste0("Percentage LESS than ", (100*meth_percentage), "% methylated: ", (100*round(nrow(subset(site_stats_rDNA, hit_call_fraction_methylated <= meth_percentage))/nrow(site_stats_rDNA),4)),"%")
  )
}


calc_txn_overlaps <- function(factorname, barcode) {
  
  txnfactor_t2_all <- data.frame()
  
  nasextension <- "//nas.ads.mwn.de/lmme/bmb/MPG/dverhagen/"
    load(paste0(nasextension, "R_related_stuff/ref_seqs.RData"))
    ref_seqs <- ref_seqs[c(1:11, 13:16)]
    print(ref_seqs)
  
    abf1_coords <- read.csv(paste0(nasextension, "R_related_stuff/external_data/abf1_coords.csv"))
    rap1_coords <- read.csv(paste0(nasextension, "R_related_stuff/external_data/rap1_coords.csv"))
    reb1_coords <- read.csv(paste0(nasextension, "R_related_stuff/external_data/reb1_coords.csv"))
    
    txn_coords <<- paste0(factorname, "_coords")
  
  # print(get0(txn_coords)[get0(txn_coords)$seqnames_arabic == "chr01", ])
  # print(subset(get0(txn_coords), seqnames_arabic == "chr01"))
  
  # meth_data_chr <- load_meth_data(barcode_subset = c("barcode08"),
  #                                 ref_seq_subset = ref_seqs)
  # seqs = "chr01"
  # barcode = "barcode08"
  for (seqs in ref_seqs) {
    print(paste0("currently calculating: ", seqs, " from ", barcode))
    
    meth_data_chr <- load_meth_data(barcode_subset = barcode,
                                    ref_seq_subset = seqs) # slow
    
    # meth_data_chr_temp <- subset(meth_data_chr, ref_seq == seqs)
    
    txn_coords_temp <- subset(get0(txn_coords), 
                              seqnames_arabic == seqs)
    
    gr_meth_data <-  makeGRangesFromDataFrame(data.frame(chr = meth_data_chr$ref_seq,
                                                         start = meth_data_chr$start,
                                                         end = meth_data_chr$end,
                                                         strand = "*"))
    #add meta data
    gr_meth_data$sample <- meth_data_chr$sample
    gr_meth_data$log_lik_ratio <- meth_data_chr$log_lik_ratio
    gr_meth_data$hit_state <- meth_data_chr$hit_state
    gr_meth_data$read_name <- meth_data_chr$read_name
    gr_meth_data$position <- meth_data_chr$position
    gr_meth_data$group_size <- meth_data_chr$group_size
    
    ###
    
    # Make Granges object for GRFs
    
    gr_txn_coords <- makeGRangesFromDataFrame(data.frame(chr = txn_coords_temp$seqnames_arabic,
                                                         start = txn_coords_temp$start - 35,
                                                         end = txn_coords_temp$end + 35,
                                                         strand = "*"))
    #add meta data
    gr_txn_coords$score <- txn_coords_temp$score
    gr_txn_coords$name <- txn_coords_temp$name
    
    # find overlaps between meth_data and gene_infos to find sites that fall within the N1 of a gene
    grf_overlaps <- findOverlapPairs(gr_meth_data, gr_txn_coords)
    grf_overlaps_df <- cbind(as.data.frame(first(grf_overlaps)), as.data.frame(second(grf_overlaps)))
    grf_overlaps_df <- subset(grf_overlaps_df, hit_state != 0) # filter out ambiguous reads
    print(paste0("Rows in dataframe for ", seqs, ": ", nrow(grf_overlaps_df)))
    
    #
    
    ###
    
    # make granges from grf_overlaps to query meth_data against it
    # this identifies a region within a read_name that overlaps with a GRF that is found to overlap between meth_data and grf_coords
    # reads that share a read_name from meth_data and grf_overlaps are the unique reads that i want to have
    gr_grf_overlaps <- makeGRangesFromDataFrame(data.frame(chr = grf_overlaps_df$seqnames,
                                                           start = (grf_overlaps_df$start.1 - 1000),
                                                           end = (grf_overlaps_df$end.1 + 1000),
                                                           strand = grf_overlaps_df$strand.1))
    gr_grf_overlaps$read_name_overlaps <- grf_overlaps_df$read_name
    gr_grf_overlaps$name_overlaps <- grf_overlaps_df$name
    gr_grf_overlaps$hit_state_overlap <- grf_overlaps_df$hit_state
    gr_grf_overlaps$position <- grf_overlaps_df$position
    gr_grf_overlaps$sample <- grf_overlaps_df$sample
    
    ### This step was originally only possible on the cluster, when running this script as a whole you can run in it locally
    grf_t1 <- findOverlapPairs(gr_meth_data, gr_grf_overlaps)
    grf_t1 <- cbind(as.data.frame(first(grf_t1)), as.data.frame(second(grf_t1)))
    print(paste0("Rows in ", seqs, " t1 dataframe: ",nrow(grf_t1)))
    
    txn_t2 <- subset(grf_t1, grf_t1$read_name == grf_t1$read_name_overlaps)
    print(paste0("Rows in ", seqs, " t2 dataframe: ",nrow(txn_t2)))
    
    # append the data in a large dataframe
    txnfactor_t2_all <- rbind(txnfactor_t2_all, txn_t2)
    print(paste0("Rows in appended dataframe after calculating ", seqs, ": ", nrow(txnfactor_t2_all)))
  }
  
  # subset genomic features that are methylated (1) or unmethylated (-1)
  txnfactor_t2_all
  txnfactor_t2_all_neg <- subset(txnfactor_t2_all, hit_state_overlap == -1)
  txnfactor_t2_all_pos <- subset(txnfactor_t2_all, hit_state_overlap == 1)
  
  # chose appropriate columns and columnmanes (required for alignment)
  meth_data_txn_all <- txnfactor_t2_all[,c("seqnames", "start", "end", "sample", "read_name", "position", "group_size", "log_lik_ratio", "hit_state")]
  colnames(meth_data_txn_all) <- c("ref_seq", "start", "end", "sample", "read_name", "position", "group_size", "log_lik_ratio", "hit_state")
  
  meth_data_txn_neg <- txnfactor_t2_all_neg[,c("seqnames", "start", "end", "sample", "read_name", "position", "group_size", "log_lik_ratio", "hit_state")]
  colnames(meth_data_txn_neg) <- c("ref_seq", "start", "end", "sample", "read_name", "position", "group_size", "log_lik_ratio", "hit_state")
  
  meth_data_txn_pos <- txnfactor_t2_all_pos[,c("seqnames", "start", "end", "sample", "read_name", "position", "group_size", "log_lik_ratio", "hit_state")]
  colnames(meth_data_txn_pos) <- c("ref_seq", "start", "end", "sample", "read_name", "position", "group_size", "log_lik_ratio", "hit_state")
  
  # write to global environment (might take out later)
  meth_data_txn_all <<- meth_data_txn_all
  meth_data_txn_neg <<- meth_data_txn_neg
  meth_data_txn_pos <<- meth_data_txn_pos
  
}

#all credits go to https://github.com/stas-g/findPeaks/ for this script
find_peaks <- function (x, m = 3){
  shape <- diff(sign(diff(x, na.pad = FALSE)))
  pks <- sapply(which(shape < 0), FUN = function(i){
    z <- i - m - 1 # only adjustment from i-m+1 to i-m-1
    z <- ifelse(z > 0, z, 1)
    w <- i + m + 1
    w <- ifelse(w < length(x), w, length(x))
    if(all(x[c(z : i, (i + 2) : w)] <= x[i + 1])) return(i + 1) else return(numeric(0))
  })
  pks <- unlist(pks)
  pks
}

