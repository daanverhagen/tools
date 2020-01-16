# it's better to just use the normal get_read_infos_from_bam.py script also for unmapped reads

import pysam
import argparse
import numpy as np

parser = argparse.ArgumentParser(description='Returns read infos from input bam file')
parser.add_argument('-i', '--input', type=str, required=True)
args = parser.parse_args()

if args.input:
  samfile = pysam.AlignmentFile(args.input, "rb")
  print 'read_name' +'\t'+ 'quality_average' +'\t'+ 'quality_std' +'\t'+ 'quality_min' +'\t'+ 'quality_perc_25' +'\t'+ 'quality_median' +'\t'+ 'quality_perc_75' +'\t'+ 'quality_max'
  
  for read in samfile.fetch():
    s1 = read.query_name +'\t'
    alignment_qualities = np.array(read.query_alignment_qualities)  # contains basecalling qualities for the aligned sequence
    
    if alignment_qualities.size > 0:
      try:
        s2 = str(alignment_qualities.mean()) +'\t'+ str(alignment_qualities.std()) +'\t'+ str(alignment_qualities.min()) +'\t'+ str(np.percentile(alignment_qualities,25)) +'\t'+ str(np.median(alignment_qualities)) +'\t'+ str(np.percentile(alignment_qualities,75)) +'\t'+ str(alignment_qualities.max())
      except TypeError:
        s2 = 'NaN' +'\t' + 'NaN' +'\t' + 'NaN' +'\t' + 'NaN' +'\t' + 'NaN' +'\t' + 'NaN' +'\t' + 'NaN'
    else:
      s2 = 'NaN' +'\t' + 'NaN' +'\t' + 'NaN' +'\t' + 'NaN' +'\t' + 'NaN' +'\t' + 'NaN' +'\t' + 'NaN'
      
    print s1 +'\t'+ s2
    
  samfile.close()

