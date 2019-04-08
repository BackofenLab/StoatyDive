
import argparse
import logging
import subprocess as sb
import matplotlib.pyplot as plt
import os
import numpy
import math

#########################
##   NECESSARY TOOLS   ##
#########################

# bedtools

####################
##   ARGS INPUT   ##
####################

tool_description = """
The tool can evalute the profile of peaks.
"""

# parse command line arguments
parser = argparse.ArgumentParser(description=tool_description, usage='%(prog)s [options]',
                                 formatter_class=argparse.RawDescriptionHelpFormatter)

# version
parser.add_argument(
    "-v", "--version", action="version", version="%(prog)s 0.0")

# positional arguments
parser.add_argument(
    "-a", "--input_bed",
    metavar='*.bed',
    required=True,
    help="Paths to the peak file in bed6 format.")
parser.add_argument(
    "-b", "--input_bam",
    metavar='*.bai',
    required=True,
    nargs="+",
    help="List of paths to the read bam files used for the peakcalling.")

# optional arguments
parser.add_argument(
    "-o", "--output_folder",
    metavar='path/',
    default=os.getcwd(),
    help="Write results to this path.")

######################
##   CHECKS INPUT   ##
######################

parser.add_argument(
    "-d", "--debug",
    help="Print lots of debugging information",
    action="store_true")

args = parser.parse_args()
if args.debug:
    logging.basicConfig(level=logging.DEBUG, format="%(asctime)s - %(filename)s - %(levelname)s - %(message)s")
else:
    logging.basicConfig(format="%(filename)s - %(levelname)s - %(message)s")
logging.info("Parsed arguments:")
if args.output_folder:
    logging.info("outfile: enabled writing to file")
    logging.info("outfile: '{}'".format(args.output_folder))
logging.info("")

##############
##   FUNC   ##
##############

def get_line_count(file):
    count = 0
    for line in file:
        count += 1
    return count

##############
##   MAIN   ##
##############

# Generate Coverage files
coverage_files = ["{}/rep{}_coverage.tsv".format(args.output_folder, i+1) for i in range(0, len(args.input_bam))]
for i in range(0, len(args.input_bam)):
    sb.Popen("bedtools coverage -a {} -b {} -d -s > {}".format(args.input_bed, args.input_bam[i], coverage_files[i]),
             shell=True).wait()

peaks_file = open(args.input_bed, "r")
num_peaks = get_line_count(peaks_file)
peaks_file.close()

print("[NOTE] {} peaks will be evaluated.".format(num_peaks))

mean_coverage_peaks = [-1] * num_peaks
variance_coverage_peaks = [-1] * num_peaks
num_bp_peaks = [-1] * num_peaks

# Evaluate peak profiles for each replicate
for file in coverage_files:
    coverage_file = open(file, "r")

    peak_counter = -1
    peak_cov_list = []

    # Calcualte mean and variance of peak coverage profiles
    for line in coverage_file:
        data = line.strip("\n").split("\t")
        bp = int(data[7])
        cov = int(data[8])

        if(bp == 1):
            peak_counter += 1
            if( peak_counter != 0 ):
                mean_coverage_peaks[peak_counter] = numpy.mean(peak_cov_list)
                variance_coverage_peaks[peak_counter] = numpy.std(peak_cov_list)
            peak_cov_list = []
            peak_cov_list.append(cov)
        else:
            peak_cov_list.append(cov)
            num_bp_peaks[peak_counter] = bp

    coverage_file.close()

    # Filter our peaks that are completly uncovered
    filtered_mean_coverage_peaks = []
    filtered_variance_coverage_peaks = []
    filtered_num_bp_peaks = []
    filtered_num_peaks = 0

    for i in range(0, len(mean_coverage_peaks)):
        if ( mean_coverage_peaks[i] != 0 ):
            filtered_mean_coverage_peaks.append(mean_coverage_peaks[i])
            filtered_variance_coverage_peaks.append(variance_coverage_peaks[i])
            filtered_num_bp_peaks.append(num_bp_peaks[i])
            filtered_num_peaks += 1

    varcoeff_coverage_peaks = [-1] * filtered_num_peaks

    print("[NOTE] {} peaks are covered.".format(filtered_num_peaks))

    # Calcualte Variantioncoefficient of peak coverage profile
    for i in range(0, filtered_num_peaks):
        varcoef = filtered_variance_coverage_peaks[i]/filtered_mean_coverage_peaks[i]

        if ( math.isnan(varcoef) ):
            print(varcoef)

        norm_varvoef = varcoef / math.sqrt(num_bp_peaks[i]-1) # Taking the estimation of the standard deviation into account
        varcoeff_coverage_peaks[i] = norm_varvoef

# Make vase plot of variationkoefficients
f = plt.figure()
plt.violinplot(varcoeff_coverage_peaks)
plt.ylabel('Variationcoefficient of the Read Coverage')
plt.show()
f.savefig(args.output_folder + "/Variationcoefficient_Distribution.pdf", bbox_inches='tight')




