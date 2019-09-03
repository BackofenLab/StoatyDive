#!/usr/bin/env python

import argparse
import subprocess as sb
import matplotlib.pyplot as plt
import os
import numpy
import math
import sys

from lib import fit_nbinom as fnb
from scipy.stats import norm

#plt.switch_backend('agg')

#########################
##   NECESSARY TOOLS   ##
#########################

# bedtools
# R + umap, data.table

##############
##   FUNC   ##
##############

# Function to obtain the number of lines in a file.
def get_line_count(file):
    count = 0
    for line in file:
        count += 1
    return count

def create_coverage_file_for_R(cov_file_name, out_matrix_file_name, output_folder, num_peaks, num_coverage_lines):
    cov_matrix_for_R = [""] * num_peaks

    coverage_file = open(cov_file_name, "r")

    peak_counter = -1
    line_count = 0

    peak_cov_list = []

    # Go over each line of the bedtools coverage file.
    for line in coverage_file:
        line_count += 1
        data = line.strip("\n").split("\t")
        bp = int(data[len(data) - 2])  # bp of the peak
        cov = str(data[len(data) - 1])  # Coverage at that bp

        # If the bp == 1 do the negative binomial estimation an start a new peak entry.
        if (bp == 1):
            if (peak_counter != -1):
                cov_matrix_for_R[peak_counter] = ",".join(peak_cov_list)
            peak_cov_list = []
            peak_cov_list.append(cov)
            peak_counter += 1
        else:
            peak_cov_list.append(cov)
            # This condition takes the last line of the coverage file into account. Else I will miss the last entry.
            if (line_count == num_coverage_lines):
                cov_matrix_for_R[peak_counter] = ",".join(peak_cov_list)

    ### write data for clustering
    outfile = open(output_folder + "/data_classification_{}.csv".format(out_matrix_file_name), "w")
    for ele in cov_matrix_for_R:
        outfile.write(ele + "\n")
    outfile.close()

def main():

    ####################
    ##   ARGS INPUT   ##
    ####################

    tool_description = """
    The tool can evalute the profile of peaks. Provide the peaks you want to evalutate in bed6 format and the reads
    you used for the peak detection in bed or bam format. The user obtains a distributions of the coefficient of variation (CV)
    which can be used to evaluate the profile landscape. In addition, the tool generates ranked list for the peaks based
    on the CV. The table hast the following columns: Chr Start End ID VC Strand bp r p Max_Norm_VC
    Left_Border_Center_Difference Right_Border_Center_Difference. See StoatyDive's development page for a detailed description.
    """

    # parse command line arguments
    parser = argparse.ArgumentParser(description=tool_description,
                                     usage='%(prog)s [-h] [options] -a *.bed -b *.bam/*bed -c *.txt',
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
        metavar='*.bam/*.bed',
        required=True,
        help="Path to the read file used for the peak calling in bed or bam format.")
    parser.add_argument(
        "-c", "--chr_file",
        metavar='*.txt',
        required=True,
        help="Path to the chromosome length file.")

    # optional arguments
    parser.add_argument(
        "-o", "--output_folder",
        metavar='path/',
        default=os.getcwd(),
        help="Write results to this path. [Default: Operating Path]")
    parser.add_argument(
        "-t", "--thresh",
        metavar='float',
        default=1.0,
        help="Set a CV threshold to divide the peak profiles into specific and unspecific. [Default: 1.0]")
    parser.add_argument(
        "--peak_correction",
        action='store_true',
        help="Activate peak correction. The peaks are recentered (shifted) for the correct sumit.")
    parser.add_argument(
        "--max_translocate",
        action='store_true',
        help="Set this flag if you want to shift the peak profiles based on the maximum value inside the profile "
             "instead of a Gaussian blur translocation.")
    parser.add_argument(
        "--peak_length",
        metavar='int',
        help="Set maximum peak length for the constant peak length.")
    parser.add_argument(
        "--max_norm_value",
        metavar='float',
        help="Provide a maximum value for CV to make the normalized CV plot more comparable.")
    parser.add_argument(
        "--border_penalty",
        action='store_true',
        help="Adds a penalty for non-centered peaks.")
    parser.add_argument(
        "--scale_max",
        metavar='float',
        help="Provide a maximum value for the CV plot.")
    parser.add_argument(
        "--maxcl",
        metavar='int',
        default=15,
        help="Maximal number of clusters of the kmeans clustering of the peak profiles. "
             "The algorithm will be optimized, i.e., the parameter is just a constraint and not absolute. [Default: 15]")
    parser.add_argument(
        "--sm",
        action='store_true',
        help="Turn on the peak profile smoothing for the peak profile classification. It is recommended to turn it on.")
    parser.add_argument(
        "--lam",
        metavar='float',
        default=0.3,
        help="Parameter for the peak profile classification. Set lambda for the smoothing of the peak profiles. "
             "A higher value (> default) will underfit. A lower value (< default) will overfit. [Default: 0.3]")
    parser.add_argument(
        "--turn_off_classification",
        action='store_true',
        help="Turn off the peak profile classification.")
    parser.add_argument(
        "--seed",
        metavar='int',
        help="Set seed for the optimization scheme.")

    ######################
    ##   CHECKS INPUT   ##
    ######################

    print("[START]")

    args = parser.parse_args()

    # Check if peak file is in bed6 format
    bedfile = open(args.input_bed, "r")
    firstline = bedfile.readline()
    if ( len(firstline.strip("\n").split("\t")) <  6 ):
        sys.exit("[ERROR] Peakfile has to be in bed6 format!")
    bedfile.close()

    ##############
    ##   MAIN   ##
    ##############

    # Get the outfile name from the input read file.
    outfilename = args.input_bam.split("/")
    outfilename = outfilename[len(outfilename)-1]
    outfilename = outfilename.replace(".bam", "").replace(".bed", "")

    # Find maximal peak length and get the number of peaks from the peak file.
    peak_file_name = args.input_bed
    peaks_file = open(peak_file_name, "r")
    max_peak_len = 0
    num_peaks = 0
    for line in peaks_file:
        num_peaks += 1
        data = line.strip("\n").split("\t")
        start = data[1]
        end = data[2]
        length = int(end) - int(start)

        if (length > max_peak_len):
            max_peak_len = length
    peaks_file.close()

    if ( args.peak_length ):
        max_peak_len = int(args.peak_length)

    if ( max_peak_len <  10 ):
        sys.exit("[ERROR] Maximal Peak Length has to be at least 10 bases.")

    # Read in chromosome sizes
    chr_sizes_dict = dict()
    chr_sizes_file = open(args.chr_file, "r")
    for line in chr_sizes_file:
        data = line.strip("\n").split("\t")
        if (data[0] not in chr_sizes_dict):
            chr_sizes_dict[data[0]] = int(data[1])
    chr_sizes_file.close()

    ##########################
    ##   Calculate Coverage ##
    ##########################

    print("[NOTE] Caluclate Coverage")

    # Generate Coverage file for central peaks with bedtools
    peak_coverage_file_name = "{}/{}_coverage_central_peaks.tsv".format(args.output_folder, outfilename)
    sb.Popen("bedtools coverage -a {} -b {} -d -s > {}".format(args.input_bed, args.input_bam, peak_coverage_file_name), shell=True).wait()

    # Get the number of lines of the coverage file of bedtools.
    coverage_file = open(peak_coverage_file_name, "r")
    num_coverage_lines = get_line_count(coverage_file)
    coverage_file.close()

    ########################
    ##   Peak Correction  ##
    ########################

    coverage_file = open(peak_coverage_file_name, "r")

    peak_counter = -1
    line_count = 0

    peak_cov_list = []
    start_list = [-100] * num_peaks
    end_list = [-100] * num_peaks

    # Go over each line of the bedtools coverage file.
    for line in coverage_file:
        line_count += 1
        data = line.strip("\n").split("\t")
        bp = int(data[len(data) - 2])  # bp of the peak
        cov = str(data[len(data) - 1])  # Coverage at that bp

        # location data
        chr = data[0]
        start = int(data[1])
        end = int(data[2])

        # If the bp == 1 do the negative binomial estimation an start a new peak entry.
        if (bp == 1 or line_count == num_coverage_lines):

            if ( line_count == num_coverage_lines ):
                peak_cov_list.append(float(cov))

            if (peak_counter != -1):
                y = numpy.array(peak_cov_list)
                x = numpy.arange(-4, 4 , 8/len(y))

                shift = 0

                if ( not args.max_translocate ):

                    dist = norm.pdf(x)
                    dist = dist/max(dist)

                    d = len(dist) - len(y)

                    # Gaussian is a bigger than peak
                    if (d > 0):
                        dist = dist[d:len(dist)]

                    # Gaussian is a bit smaller than peak
                    if (d < 0):
                        zero_vector = numpy.zeros(len(y))
                        zero_vector[1:len(dist)] = dist

                        l = numpy.ceil(abs(d)/2)
                        r = numpy.ceil(abs(d)/2)
                        if ( (l+r+len(dist)) != len(y) ):
                            r += len(y) - (l+r+len(dist))
                        dist = numpy.pad(dist, (l, r), 'constant', constant_values=(0,0))

                    # Check if Gaussian and peak have the same length
                    if ( len(dist) != len(y) ):
                        print("[ERROR] Gaussian has the wrong dimension.")
                        print("Peak " + str(peak_counter))
                        print("Length of Gaussian " + str(len(dist)))
                        print("Length of Peak " + str(len(y)))

                    shift = numpy.argmax(numpy.convolve(dist, y)) - len(y)

                else:

                    shift = numpy.argmax(y) - numpy.floor(len(y)/2)

                new_start = start + shift
                new_end = end + shift

                # check if I dont violate any borders
                if ( new_start < 0 ):
                    new_start = 0
                if ( new_end > chr_sizes_dict[chr] ):
                    new_end = chr_sizes_dict[chr]

                start_list[peak_counter] = new_start
                end_list[peak_counter] = new_end

                if ( peak_counter == num_peaks ):
                    print(y)

            peak_cov_list = []
            peak_cov_list.append(float(cov))
            peak_counter += 1
        else:
            peak_cov_list.append(float(cov))

    ########################
    ##   Peak Extension   ##
    ########################
    #
    # extended_peak_file_name = "{}/peaks_extended.bed".format(args.output_folder)
    #
    # # Extend the peaks to the maximal length if the parameter is set to true.
    # print("[NOTE] Maximal peak length {}.".format(max_peak_len))
    #
    # # Define new coorindate for peaks. Extend to maximal length.
    # peaks_file = open(args.input_bed, "r")
    # extended_peak_file = open(extended_peak_file_name, "w")
    #
    # peak_counter = 0
    # for line in peaks_file:
    #     data = line.strip("\n").split("\t")
    #     start = int(data[1])
    #     end = int(data[2])
    #
    #     # Apply peak shift if the correction is activated
    #     if ( args.peak_correction ):
    #         start = start_list[peak_counter]
    #         end = end_list[peak_counter]
    #
    #     peak_length = end - start
    #     extention_left = numpy.round((max_peak_len - peak_length) / 2)
    #     extentions_right = numpy.round((max_peak_len - peak_length) / 2)
    #
    #     # Check if extention left and right make up the max_peak_length, if not,
    #     # then add or substract randomly either to left or right some extra bases. This happends
    #     # because of the rounding.
    #     current_peak_length = extention_left + extentions_right + peak_length
    #     if (current_peak_length < max_peak_len):
    #         # Set seed if seed is provided.
    #         if (args.seed):
    #             numpy.random.seed(int(args.seed))
    #
    #         if (numpy.random.randint(low=2, size=1) == 0):
    #             extention_left += max_peak_len - current_peak_length
    #         else:
    #             extentions_right += max_peak_len - current_peak_length
    #
    #     if (current_peak_length > max_peak_len):
    #         # Set seed if seed is provided.
    #         if (args.seed):
    #             numpy.random.seed(int(args.seed))
    #
    #         if (numpy.random.randint(low=2, size=1) == 0):
    #             extention_left -= current_peak_length - max_peak_len
    #         else:
    #             extentions_right -= current_peak_length - max_peak_len
    #
    #     # Check if extension goes beyond the borders of the chromosome.
    #     beyond_left = "false"
    #     if ((start - extention_left) < 0):
    #         beyond_left = "true"
    #     beyond_right = "false"
    #     if ((end + extentions_right) > chr_sizes_dict[data[0]]):
    #         beyond_right = "true"
    #
    #     if (beyond_left == "true" and beyond_right == "false"):
    #         extentions_right += extention_left - start
    #         extention_left = start
    #
    #     if (beyond_left == "false" and beyond_right == "true"):
    #         extention_left += (end + extentions_right) - chr_sizes_dict[data[0]]
    #         extentions_right = chr_sizes_dict[data[0]] - end
    #
    #     if (beyond_left == "true" and beyond_right == "true"):
    #         extention_left = start
    #         extentions_right = chr_sizes_dict[data[0]] - end
    #
    #     start = start - extention_left
    #     end = end + extentions_right
    #
    #     # A last checkup if peak length is maximum length.
    #     if ((end - start) != max_peak_len and not (beyond_left == "true" and beyond_left == "true")):
    #         print("[ERROR] Max length of peaks not reached.")
    #         print(data)
    #         print(start)
    #         print(end)
    #         print(end - start)
    #         print(max_peak_len)
    #
    #     # Write extended peak to file.
    #     extended_peak_file.write("{}\t{}\t{}\t{}\n".format(data[0], int(start), int(end), "\t".join(data[3:])))
    #
    # peaks_file.close()
    # extended_peak_file.close()
    # peak_file_name = extended_peak_file_name

    ###################################
    ##   Generating Cov matrix for R ##
    ###################################

    # Central Peaks
    #create_coverage_file_for_R(peak_coverage_file_name, "peaks_" + outfilename, args.output_folder, num_peaks, num_coverage_lines)

    # Windows
    # if (args.peak_correction):
    #     # Get the number of lines of the coverage file of bedtools.
    #     coverage_file_windows = open(window_coverage_file_name, "r")
    #     num_coverage_lines_windows = get_line_count(coverage_file_windows)
    #     coverage_file.close()
    #
    #     create_coverage_file_for_R(window_coverage_file_name, "windows_" + outfilename, args.output_folder, num_peaks, num_coverage_lines_windows)


    ################################
    ##   Peak Correction and uMAP ##
    ################################

    # Add flag in uMAP.R for skipping classifcation


    ####################
    ##   Calculate CV ##
    ####################

    # print("[NOTE] {} peaks will be evaluated.".format(num_peaks))
    #
    # # Dictionaries for the algorithm.
    # size_r_peaks_dict = dict()
    # prob_success_peaks_dict = dict()
    # variance_coverage_peaks_dict = dict()
    # num_bp_peaks_dict = dict()
    # coordinates_dict = dict()
    # peakid_dict = dict()
    # strand_dict = dict()
    # center_border_diff_left_dict = dict()
    # center_border_diff_right_dict = dict()
    #
    # # Calculate mean and variance of peak coverage profiles
    # peak_cov_list = []
    #
    # coverage_file = open(coverage_file_name, "r")
    #
    # peak_counter = -1
    # line_count = 0
    #
    # # Go over each line of the bedtools coverage file.
    # for line in coverage_file:
    #     line_count += 1
    #     data = line.strip("\n").split("\t")
    #     bp = int(data[len(data)-2])     # bp of the peak
    #     cov = int(data[len(data)-1])    # Coverage at that bp
    #     peakid = data[3]
    #     strand = data[5]
    #
    #     # If the bp == 1 do the negative binomial estimation an start a new peak entry.
    #     if(bp == 1):
    #         if( peak_counter != -1 ):
    #             peak_center = int(num_bp_peaks_dict[peak_counter] / 2)
    #             peak_center_ext = int(num_bp_peaks_dict[peak_counter] * 0.1)
    #             # The fit is the alternative version of the NB. But I get the number of successes (r) and the
    #             # probability of success (p). At least one value needs to be greater than zero, else the estimation
    #             # makes no sense.
    #             # The second condition tests for a centered peak. This filters out peaks that have spkied read coverages
    #             # at the peak ends.
    #             center_region = peak_cov_list[(peak_center-peak_center_ext-1):(peak_center+peak_center_ext)]
    #             border_left = peak_cov_list[0:peak_center_ext]
    #             border_right = peak_cov_list[-peak_center_ext:]
    #
    #             if ( not all(v == 0 for v in peak_cov_list) and not all(v <= 10 for v in center_region) ):
    #                 # predefined features
    #                 nb_fit = fnb.fit_nbinom(numpy.array(peak_cov_list))
    #                 size_r_peaks_dict[peak_counter] = nb_fit["size"]
    #                 prob_success_peaks_dict[peak_counter] = nb_fit["prob"]
    #                 variance_coverage_peaks_dict[peak_counter] = (nb_fit["size"] * (1-nb_fit["prob"])) / (nb_fit["prob"] * nb_fit["prob"])
    #                 center_border_diff_left_dict[peak_counter] = numpy.max(center_region) - numpy.max(border_left)
    #                 center_border_diff_right_dict[peak_counter] = numpy.max(center_region) - numpy.max(border_right)
    #             else:
    #                 size_r_peaks_dict[peak_counter] = 0.0
    #                 prob_success_peaks_dict[peak_counter] = 0.0
    #                 variance_coverage_peaks_dict[peak_counter] = 0.0
    #                 center_border_diff_left_dict[peak_counter] = 0.0
    #                 center_border_diff_right_dict[peak_counter] = 0.0
    #
    #         peak_cov_list = []
    #         peak_cov_list.append(cov)
    #         peak_counter += 1
    #         coordinates_dict[peak_counter] = [data[0], data[1], data[2]]
    #     else:
    #         peak_cov_list.append(cov)
    #         num_bp_peaks_dict[peak_counter] = bp
    #         strand_dict[peak_counter] = strand
    #         peakid_dict[peak_counter] = peakid
    #
    #         # This condition takes the last line of the coverage file into account. Else I will miss the last entry.
    #         if ( line_count == num_coverage_lines ):
    #             peak_center = int(num_bp_peaks_dict[peak_counter] / 2)
    #             peak_center_ext = int(num_bp_peaks_dict[peak_counter] * 0.1)
    #             center_region = peak_cov_list[(peak_center - peak_center_ext - 1):(peak_center + peak_center_ext)]
    #             border_left = peak_cov_list[0:peak_center_ext]
    #             border_right = peak_cov_list[-peak_center_ext:]
    #
    #             if (not all(v == 0 for v in peak_cov_list) and not all(v <= 10 for v in center_region)):
    #                 # predefined features
    #                 nb_fit = fnb.fit_nbinom(numpy.array(peak_cov_list))
    #                 size_r_peaks_dict[peak_counter] = nb_fit["size"]
    #                 prob_success_peaks_dict[peak_counter] = nb_fit["prob"]
    #                 variance_coverage_peaks_dict[peak_counter] = (nb_fit["size"] * (1 - nb_fit["prob"])) / (nb_fit["prob"] * nb_fit["prob"])
    #                 center_border_diff_left_dict[peak_counter] = numpy.max(center_region) - numpy.max(border_left)
    #                 center_border_diff_right_dict[peak_counter] = numpy.max(center_region) - numpy.max(border_right)
    #             else:
    #                 size_r_peaks_dict[peak_counter] = 0.0
    #                 prob_success_peaks_dict[peak_counter] = 0.0
    #                 variance_coverage_peaks_dict[peak_counter] = 0.0
    #                 center_border_diff_left_dict[peak_counter] = 0.0
    #                 center_border_diff_right_dict[peak_counter] = 0.0
    #
    # coverage_file.close()
    #
    # filtered_num_peaks = 0
    # varcoeff_coverage_peaks_dict = dict()
    #
    # if (args.border_penalty):
    #     print("[NOTE] Activate border penalty.")
    #
    # # Calcualte Variantioncoefficient of peak coverage profile.
    # for i in range(0, num_peaks):
    #
    #     # The mean coverage has to be greater than zero or else the VC is not defined.
    #     if (size_r_peaks_dict[i] > 0):
    #         # varcoeff_coverage_peaks_dict[i] = 1 / math.sqrt(size_r_peaks_dict[i] * (1 - prob_success_peaks_dict[i]))
    #         varcoeff_coverage_peaks_dict[i] = math.sqrt((1-prob_success_peaks_dict[i])/size_r_peaks_dict[i])
    #
    #         if ( args.border_penalty ):
    #             w1 = 1
    #             w2 = 1
    #             if ( center_border_diff_left_dict[i] < 0 ):
    #                 w1 = 1/abs(center_border_diff_left_dict[i])
    #             if ( center_border_diff_right_dict[i] < 0 ):
    #                 w2 = 1/abs(center_border_diff_right_dict[i])
    #             varcoeff_coverage_peaks_dict[i] = varcoeff_coverage_peaks_dict[i] * w1 * w2
    #
    #         # Just a safety condition.
    #         if ( math.isnan(varcoeff_coverage_peaks_dict[i]) ):
    #             print(varcoeff_coverage_peaks_dict[i])
    #
    #         filtered_num_peaks += 1
    #     else:
    #         varcoeff_coverage_peaks_dict[i] = -0.01
    #
    # #######################
    # ##   Generate Output ##
    # #######################
    #
    # print("[NOTE] {} peaks are covered.".format(filtered_num_peaks))
    #
    # # Filter our peaks that are completly uncovered.
    # filtered_varcoeff_coverage_peaks = []
    # for i in varcoeff_coverage_peaks_dict.values():
    #     if (i >= 0):
    #         filtered_varcoeff_coverage_peaks.append(i)
    #
    # print("[NOTE] Generate CV plot.")
    #
    # scale_max = numpy.max(filtered_varcoeff_coverage_peaks)
    # if ( args.scale_max ):
    #     if ( float(args.scale_max) < 0.0 ):
    #         sys.exit("[ERROR] Wrong value for scale_max!")
    #     scale_max = float(args.scale_max)
    #
    # # Make vase plot of variationkoefficients.
    # plt.violinplot(filtered_varcoeff_coverage_peaks)
    # plt.ylim(0.0, scale_max)
    # plt.xticks([], 'none')
    # plt.ylabel('Coefficient of Variation of the Peak Profiles')
    # plt.savefig(args.output_folder + "/CV_Distribution_{}.pdf".format(outfilename), bbox_inches='tight')
    #
    # # Normalize all VC so that the scale goes from 0 to 1 which makes a comparison between
    # # different profile evaluations easier. Unity-based normalization.
    # one = numpy.max(filtered_varcoeff_coverage_peaks)
    # zero = numpy.min(filtered_varcoeff_coverage_peaks)
    # if ( args.max_norm_value ):
    #     one = float(args.max_norm_value)
    #
    # for i in range(0, len(filtered_varcoeff_coverage_peaks)):
    #     filtered_varcoeff_coverage_peaks[i] = (filtered_varcoeff_coverage_peaks[i]-zero)/(one-zero)
    #
    # print("[NOTE] Generate normalized CV plot.")
    #
    # # Make vase plot of variationkoefficients.
    # f2 = plt.figure()
    # plt.violinplot(filtered_varcoeff_coverage_peaks)
    # plt.ylim(0.0, 1.0)
    # plt.xticks([], 'none')
    # plt.ylabel('Normalized Coefficient of Variation of the Peak Profiles')
    # f2.savefig(args.output_folder + "/Norm_CV_Distribution_{}.pdf".format(outfilename), bbox_inches='tight')
    #
    # # Generate the output tabular file.
    # print("[NOTE] Generate output tabular.")
    # out_tab_file = open(args.output_folder + "/CV_tab_{}.bed".format(outfilename), "w")
    #
    # index_sort = numpy.argsort(list(varcoeff_coverage_peaks_dict.values()))[::-1]
    # keys_list = list(varcoeff_coverage_peaks_dict.keys())
    #
    # for i in index_sort:
    #     k = keys_list[i]
    #     coords = coordinates_dict[k]
    #
    #     max_norma_cv = (varcoeff_coverage_peaks_dict[k]-zero)/(one-zero)
    #
    #     type = 0
    #     if ( varcoeff_coverage_peaks_dict[k] >= float(args.thresh)):
    #         type = 0
    #     else:
    #         type = 1
    #
    #     # Chr  Start   End
    #     # ID_from_File    CV  Strand
    #     # bp    r   p
    #     # Max_Norm_CV
    #     # Left_Border_Center_Difference
    #     # Right_Border_Center_Difference Inter_Index
    #     # Internal_Peak_ID  Specific/Unspecific"
    #     out_tab_file.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\t{12}\t{13}\n".format(coords[0], coords[1], coords[2],
    #                                                              peakid_dict[k], varcoeff_coverage_peaks_dict[k], strand_dict[k],
    #                                                              num_bp_peaks_dict[k], size_r_peaks_dict[k], prob_success_peaks_dict[k],
    #                                                              max_norma_cv,
    #                                                              center_border_diff_left_dict[k],
    #                                                              center_border_diff_right_dict[k],
    #                                                              k+1, type))
    # out_tab_file.close()
    #
    #
    # ### write data for clustering
    # if ( bool_length_norm == 0 ):
    #     numpy.savetxt(args.output_folder + "/data_classification_{}.tsv".format(outfilename), tsne_matrix, delimiter="\t", newline="\n")
    # else:
    #     # Pad profiles with zero to make them the same length if length normalization was tuned of
    #     for v in tsne_matrix:
    #         l = numpy.floor( (max_peak_len - len(v))/2 )
    #         r = numpy.floor( (max_peak_len - len(v))/2 )
    #
    #         if ( (l+r+len(v)) != max_peak_len ):
    #             r += max_peak_len - (l+r+len(v))
    #
    #         v = numpy.pad(v, (l, r), 'constant', constant_values=(0,0))
    #     numpy.savetxt(args.output_folder + "/data_classification_{}.tsv".format(outfilename), tsne_matrix, delimiter="\t", newline="\n")
    #
    # script_path = os.path.join(os.path.dirname(os.path.dirname(fnb.__file__)), "lib/uMAP.R")
    #
    # ### Create R call
    # Rcall = "Rscript {0} {1} {2} {3} {4}".format(script_path, args.output_folder, outfilename, args.lam, args.maxcl)
    #
    # if ( args.sm ):
    #     print("[NOTE] Smoothing turned on")
    #     Rcall += " TRUE"
    # else:
    #     Rcall += " FALSE"
    #
    # if ( args.max_translocate ):
    #     print("[NOTE] Max translocation turned on")
    #     Rcall += " TRUE"
    # else:
    #     Rcall += " FALSE"
    #
    # ### Run Rscript for classification
    # if ( args.turn_off_classification ):
    #     print("[NOTE] Skip Classification")
    # else:
    #     print("[NOTE] Run Classification")
    #     sb.Popen(Rcall, shell=True).wait()

    print("[FINISH]")

if __name__ == '__main__':
    main()
