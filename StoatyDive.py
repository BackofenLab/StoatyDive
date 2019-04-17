
import argparse
import subprocess as sb
import matplotlib.pyplot as plt
import os
import numpy
import math
import fit_nbinom as fnb
import sys

plt.switch_backend('agg')

#########################
##   NECESSARY TOOLS   ##
#########################

# bedtools

####################
##   ARGS INPUT   ##
####################

tool_description = """
The tool can evalute the profile of peaks. Provide as an input the peaks you want to evalutate
in bed6 format and the reads you used for the peak detection in bed or bam format. The user obtains
a distributions of the variationcoefficient (VC) which can be used to evaluate the profile landscape. 
In addition, the tool generates a list ranking the peaks based on the VC. The VC range from 0.0 for a broad 
to 1.0 for a sharp peak. The ranked list has the following columns: chr, start, end, peakid, VC, strand,
peak length, mean read coverage. 
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
    help="Path to the peak file in bed6 format.")
parser.add_argument(
    "-b", "--input_bam",
    metavar='*.bai',
    required=True,
    help="Path to the read bam file used for the peakcalling in bed or bam format.")

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

# Check if peak file is in bed6 format
bedfile = open(args.input_bed, "r")
firstline = bedfile.readline()
if ( len(firstline.strip("\n").split("\t")) <  6 ):
    sys.exit("[ERROR] Peakfile has to be in bed6 format!")
bedfile.close()

##############
##   FUNC   ##
##############

# Function to obtain the number of lines in a file.
def get_line_count(file):
    count = 0
    for line in file:
        count += 1
    return count

##############
##   MAIN   ##
##############

# Get the outfile name from the input read file.
outfilename = args.input_bam.split("/")
outfilename = outfilename[len(outfilename)-1]
outfilename = outfilename.strip(".bam")
outfilename = outfilename.strip(".bed")

# Generate Coverage file with bedtools
coverage_file_name = "{}/{}_coverage.tsv".format(args.output_folder, outfilename)
sb.Popen("bedtools coverage -a {} -b {} -d -s > {}".format(args.input_bed, args.input_bam, coverage_file_name), shell=True).wait()

# Get the number of peaks from the peak file.
peaks_file = open(args.input_bed, "r")
num_peaks = get_line_count(peaks_file)
peaks_file.close()

print("[NOTE] {} peaks will be evaluated.".format(num_peaks))

# Dictionaries for the algorithm.
mean_coverage_peaks_dict = dict()
prob_success_peaks_dict = dict()
variance_coverage_peaks_dict = dict()
num_bp_peaks_dict = dict()
coordinates_dict = dict()

# Get the number of lines of the coverage file of bedtools.
coverage_file = open(coverage_file_name, "r")
num_coverage_lines = get_line_count(coverage_file)
coverage_file.close()

# Calcualte mean and variance of peak coverage profiles
peak_cov_list = []

coverage_file = open(coverage_file_name, "r")

peak_counter = -1
line_count = 0

# Go over each line of the bedtools coverage file.
for line in coverage_file:
    line_count += 1
    data = line.strip("\n").split("\t")
    bp = int(data[len(data)-2])     # bp of the peak
    cov = int(data[len(data)-1])    # Coverage at that bp

    # If the bp == 1 do the negative binomial estimation an start a new peak entry.
    if(bp == 1):
        if( peak_counter != -1 ):
            # The fit is the alternative version of the NB. But I get the expected number of successes and the
            # probability of success. At least one value needs to be greater than zero, else the estimation
            # makes no sense.
            if ( not all(v == 0 for v in peak_cov_list) ):
                nb_fit = fnb.fit_nbinom(numpy.array(peak_cov_list))
                mean_coverage_peaks_dict[peak_counter] = nb_fit["size"]
                prob_success_peaks_dict[peak_counter] = nb_fit["prob"]
                variance_coverage_peaks_dict[peak_counter] = (nb_fit["size"] * (1-nb_fit["prob"])) / (nb_fit["prob"] * nb_fit["prob"])
            else:
                mean_coverage_peaks_dict[peak_counter] = 0.0
                prob_success_peaks_dict[peak_counter] = 0.0
                variance_coverage_peaks_dict[peak_counter] = 0.0
        peak_cov_list = []
        peak_cov_list.append(cov)
        peak_counter += 1
        coordinates_dict[peak_counter] = [data[0], data[1], data[2]]
    else:
        peak_cov_list.append(cov)
        num_bp_peaks_dict[peak_counter] = bp

        # This condition takes the last line of the coverage file into account. Else I will miss the last entry.
        if ( line_count == num_coverage_lines ):
            if (  not all(v == 0 for v in peak_cov_list) ):
                nb_fit = fnb.fit_nbinom(numpy.array(peak_cov_list))
                mean_coverage_peaks_dict[peak_counter] = nb_fit["size"]
                prob_success_peaks_dict[peak_counter] = nb_fit["prob"]
                variance_coverage_peaks_dict[peak_counter] = (nb_fit["size"] * (1 - nb_fit["prob"])) / (nb_fit["prob"] * nb_fit["prob"])
            else:
                mean_coverage_peaks_dict[peak_counter] = 0.0
                prob_success_peaks_dict[peak_counter] = 0.0
                variance_coverage_peaks_dict[peak_counter] = 0.0

coverage_file.close()

filtered_num_peaks = 0
varcoeff_coverage_peaks_dict = dict()

# Calcualte Variantioncoefficient of peak coverage profile.
for i in range(0, num_peaks):

    # The mean coverage has to be greater than zero or else the VC is not defined.
    if (mean_coverage_peaks_dict[i] > 0):
        varcoef = 1 / math.sqrt(mean_coverage_peaks_dict[i] * (1 - prob_success_peaks_dict[i]))

        # Just a safety condition.
        if ( math.isnan(varcoef) ):
            print(varcoef)

        # Because the standard deviation (variance) was estimated, I have to correct the VC based on the changed number
        # of the degrees of freedom.
        varcoeff_coverage_peaks_dict[i] = varcoef / math.sqrt(num_bp_peaks_dict[i]-1)

        filtered_num_peaks += 1
    else:
        varcoeff_coverage_peaks_dict[i] = -0.01

print("[NOTE] {} peaks are covered.".format(filtered_num_peaks))

# Filter our peaks that are completly uncovered.
filtered_varcoeff_coverage_peaks = []
for i in varcoeff_coverage_peaks_dict.values():
    if (i >= 0):
        filtered_varcoeff_coverage_peaks.append(i)


# Normalize all VC so that the scale goes from 0 to 1 which makes a comparison between
# different profile evaluations easier.
one = numpy.max(filtered_varcoeff_coverage_peaks)

for i in range(0, len(filtered_varcoeff_coverage_peaks)):
    filtered_varcoeff_coverage_peaks[i] = filtered_varcoeff_coverage_peaks[i]/one

print("[NOTE] Generate Plot")

# Make vase plot of variationkoefficients.
f = plt.figure()
plt.violinplot(filtered_varcoeff_coverage_peaks)
plt.ylabel('Normalized Variationcoefficient of the Peak Profiles')
f.savefig(args.output_folder + "/VC_Distribution_{}.pdf".format(outfilename), bbox_inches='tight')

# Generate the output tabular file.
print("[NOTE] Generate Output Tabular")
out_tab_file_name = args.output_folder + "/VC_tab_{}_tmp.bed".format(outfilename)
out_tab_file = open(out_tab_file_name, "w")
for i in range(0, num_peaks):
    coords = coordinates_dict[i]
    out_tab_file.write("{}\t{}\t{}\tpeak_{}\t{}\t.\t{}\t{}\n".format(coords[0], coords[1], coords[2],
                                                             i, varcoeff_coverage_peaks_dict[i],
                                                             num_bp_peaks_dict[i], mean_coverage_peaks_dict[i]))
out_tab_file.close()

# Sort the tabular file.
sb.Popen("sort -r -k 5 {} > {}".format(out_tab_file_name, "{}/VC_tab_{}.bed".format(args.output_folder, outfilename)), shell=True).wait()
sb.Popen("rm {}".format(out_tab_file_name), shell=True).wait()