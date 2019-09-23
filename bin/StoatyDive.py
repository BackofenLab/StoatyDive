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
# R + umap, data.table, zoo

##############
##   FUNC   ##
##############

# Function to obtain the number of lines in a file.
def get_line_count(file):
    count = 0
    for line in file:
        count += 1
    return count

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
        default=0.5,
        help="Set a normalized CV threshold to divide the peak profiles into more specific (0) and more "
             "unspecific (1). [Default: 1.0]")
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
        "-k", "--numcl",
        metavar='int',
        default=-1,
        help="You can forcefully set the number of cluster of peak profiles.")
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

    ######################
    ##   CHECKS INPUT   ##
    ######################

    print("[START]")

    seed = 123

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

    ########################
    ##   Peak Correction  ##
    ########################

    if ( args.peak_correction ):

        print("[NOTE] Activate Peak Correction")

        if (args.max_translocate):
            print("[NOTE] Activate Maximal Value Shifting")

        # Generate Coverage file for central peaks with bedtools
        peak_coverage_file_name = "{}/{}_coverage_central_peaks.tsv".format(args.output_folder, outfilename)
        sb.Popen("bedtools coverage -a {} -b {} -d -s > {}".format(args.input_bed, args.input_bam, peak_coverage_file_name), shell=True).wait()

        # Get the number of lines of the coverage file of bedtools.
        central_peak_coverage_file = open(peak_coverage_file_name, "r")
        central_peak_num_coverage_lines = get_line_count(central_peak_coverage_file)
        central_peak_coverage_file.close()

        central_peak_coverage_file = open(peak_coverage_file_name, "r")

        peak_counter = -1
        line_count = 0

        peak_cov_list = []
        start_list = [-100] * num_peaks
        end_list = [-100] * num_peaks

        chr = -1
        start = -1
        end = -1

        # Go over each line of the bedtools coverage file.
        for line in central_peak_coverage_file:
            line_count += 1
            data = line.strip("\n").split("\t")
            bp = int(data[len(data) - 2])  # bp of the peak
            cov = str(data[len(data) - 1])  # Coverage at that bp

            # If the bp == 1 do the negative binomial estimation an start a new peak entry.
            if ( bp == 1 or line_count == central_peak_num_coverage_lines ):

                if ( line_count == central_peak_num_coverage_lines ):
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

                    # For Debugging the Peak correction
                    # if ( peak_counter <= 20 or peak_counter == (num_peaks-1) ):
                    #     fig = plt.figure()
                    #     plt.plot(peak_cov_list)
                    #     plt.title(str("Shift:{}, Start:{}, NewS:{}, End:{}, NewE:{}".format(str(shift), str(start), str(new_start), str(end), str(new_end))))
                    #     plt.axvline(x=int(len(y)/2), c="b")
                    #     plt.axvline(x=int(len(y) / 2) + shift, c="r")
                    #     fig.savefig(args.output_folder + "/{}.pdf".format(peak_counter+1), bbox_inches='tight')

                peak_cov_list = []
                peak_cov_list.append(float(cov))
                peak_counter += 1

                # location data
                chr = data[0]
                start = int(data[1])
                end = int(data[2])
            else:
                peak_cov_list.append(float(cov))

    ########################
    ##   Peak Extension   ##
    ########################

    extended_peak_file_name = "{}/peaks_extended.bed".format(args.output_folder)

    # Extend the peaks to the maximal length if the parameter is set to true.
    print("[NOTE] Maximal peak length {}.".format(max_peak_len))

    # Define new coorindate for peaks. Extend to maximal length.
    peaks_file = open(args.input_bed, "r")
    extended_peak_file = open(extended_peak_file_name, "w")

    peak_counter = 0
    for line in peaks_file:
        data = line.strip("\n").split("\t")
        start = int(data[1])
        end = int(data[2])

        # Apply peak shift if the correction is activated
        if ( args.peak_correction ):
            start = start_list[peak_counter]
            end = end_list[peak_counter]
            peak_counter += 1

        peak_length = end - start
        extention_left = numpy.round((max_peak_len - peak_length) / 2)
        extentions_right = numpy.round((max_peak_len - peak_length) / 2)

        # Check if extention left and right make up the max_peak_length, if not,
        # then add or substract randomly either to left or right some extra bases. This happends
        # because of the rounding.
        current_peak_length = extention_left + extentions_right + peak_length
        if (current_peak_length < max_peak_len):
            numpy.random.seed(int(seed))

            if (numpy.random.randint(low=2, size=1) == 0):
                extention_left += max_peak_len - current_peak_length
            else:
                extentions_right += max_peak_len - current_peak_length

        if (current_peak_length > max_peak_len):
            numpy.random.seed(int(seed))

            if (numpy.random.randint(low=2, size=1) == 0):
                extention_left -= current_peak_length - max_peak_len
            else:
                extentions_right -= current_peak_length - max_peak_len

        # Check if extension goes beyond the borders of the chromosome.
        beyond_left = "false"
        if ((start - extention_left) < 0):
            beyond_left = "true"
        beyond_right = "false"
        if ((end + extentions_right) > chr_sizes_dict[data[0]]):
            beyond_right = "true"

        if (beyond_left == "true" and beyond_right == "false"):
            extentions_right += extention_left - start
            extention_left = start

        if (beyond_left == "false" and beyond_right == "true"):
            extention_left += (end + extentions_right) - chr_sizes_dict[data[0]]
            extentions_right = chr_sizes_dict[data[0]] - end

        if (beyond_left == "true" and beyond_right == "true"):
            extention_left = start
            extentions_right = chr_sizes_dict[data[0]] - end

        start = start - extention_left
        end = end + extentions_right

        # A last checkup if peak length is maximum length.
        if ((end - start) != max_peak_len and not (beyond_left == "true" and beyond_left == "true")):
            print("[ERROR] Max length of peaks not reached.")
            print(data)
            print(start)
            print(end)
            print(end - start)
            print(max_peak_len)

        # Write extended peak to file.
        extended_peak_file.write("{}\t{}\t{}\t{}\n".format(data[0], int(start), int(end), "\t".join(data[3:])))

    peaks_file.close()
    extended_peak_file.close()

    ##########################
    ##   Calculate Coverage ##
    ##########################

    print("[NOTE] Calculate Coverage")

    # Generate Coverage file for central peaks with bedtools
    coverage_file_name = "{}/{}_coverage.tsv".format(args.output_folder, outfilename)
    sb.Popen("bedtools coverage -a {} -b {} -d -s > {}".format(extended_peak_file_name, args.input_bam, coverage_file_name), shell=True).wait()

    # Get the number of lines of the coverage file of bedtools.
    coverage_file = open(coverage_file_name, "r")
    num_coverage_lines = get_line_count(coverage_file)
    coverage_file.close()

    ####################
    ##   Calculate CV ##
    ####################

    print("[NOTE] {} peaks will be evaluated.".format(num_peaks))

    # Dictionaries for the algorithm.
    size_r_peaks_dict = dict()
    prob_success_peaks_dict = dict()
    variance_coverage_peaks_dict = dict()
    num_bp_peaks_dict = dict()
    coordinates_dict = dict()
    peakid_dict = dict()
    strand_dict = dict()
    center_border_diff_left_dict = dict()
    center_border_diff_right_dict = dict()

    # Calculate mean and variance of peak coverage profiles
    peak_cov_list = []
    cov_matrix_for_R = numpy.empty([num_peaks, max_peak_len])

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
        if( bp == 1 or line_count == num_coverage_lines ):

            if (line_count == num_coverage_lines):
                peak_cov_list.append(float(cov))

            if( peak_counter != -1 ):

                cov_matrix_for_R[peak_counter] = numpy.array(peak_cov_list)

                peak_center = int(num_bp_peaks_dict[peak_counter] / 2)
                peak_center_ext = int(num_bp_peaks_dict[peak_counter] * 0.1)
                # The fit is the alternative version of the NB. But I get the number of successes (r) and the
                # probability of success (p). At least one value needs to be greater than zero, else the estimation
                # makes no sense.
                # The second condition tests for a centered peak. This filters out peaks that have spkied read coverages
                # at the peak ends.
                center_region = peak_cov_list[(peak_center-peak_center_ext-1):(peak_center+peak_center_ext)]
                border_left = peak_cov_list[0:peak_center_ext]
                border_right = peak_cov_list[-peak_center_ext:]

                if ( not all(v == 0 for v in peak_cov_list) and not all(v <= 10 for v in center_region) ):
                    # predefined features
                    nb_fit = fnb.fit_nbinom(numpy.array(peak_cov_list))
                    size_r_peaks_dict[peak_counter] = nb_fit["size"]
                    prob_success_peaks_dict[peak_counter] = nb_fit["prob"]
                    variance_coverage_peaks_dict[peak_counter] = (nb_fit["size"] * (1-nb_fit["prob"])) / (nb_fit["prob"] * nb_fit["prob"])
                    center_border_diff_left_dict[peak_counter] = numpy.max(center_region) - numpy.max(border_left)
                    center_border_diff_right_dict[peak_counter] = numpy.max(center_region) - numpy.max(border_right)
                else:
                    size_r_peaks_dict[peak_counter] = 0.0
                    prob_success_peaks_dict[peak_counter] = 0.0
                    variance_coverage_peaks_dict[peak_counter] = 0.0
                    center_border_diff_left_dict[peak_counter] = 0.0
                    center_border_diff_right_dict[peak_counter] = 0.0

            peak_cov_list = []
            peak_cov_list.append(float(cov))
            peak_counter += 1

            # Get coordination data
            coordinates_dict[peak_counter] = [data[0], data[1], data[2]]
            strand_dict[peak_counter] = data[5]
            peakid_dict[peak_counter] = data[3]
        else:
            peak_cov_list.append(float(cov))
            num_bp_peaks_dict[peak_counter] = bp

    coverage_file.close()

    filtered_num_peaks = 0
    varcoeff_coverage_peaks_dict = dict()

    if (args.border_penalty):
        print("[NOTE] Activate border penalty.")

    # Calcualte Variantioncoefficient of peak coverage profile.
    for i in range(0, num_peaks):

        # The mean coverage has to be greater than zero or else the VC is not defined.
        if (size_r_peaks_dict[i] > 0):
            # varcoeff_coverage_peaks_dict[i] = 1 / math.sqrt(size_r_peaks_dict[i] * (1 - prob_success_peaks_dict[i]))
            varcoeff_coverage_peaks_dict[i] = math.sqrt((1-prob_success_peaks_dict[i])/size_r_peaks_dict[i])

            if ( args.border_penalty ):
                w1 = 1
                w2 = 1
                if ( center_border_diff_left_dict[i] < 0 ):
                    w1 = 1/abs(center_border_diff_left_dict[i])
                if ( center_border_diff_right_dict[i] < 0 ):
                    w2 = 1/abs(center_border_diff_right_dict[i])
                varcoeff_coverage_peaks_dict[i] = varcoeff_coverage_peaks_dict[i] * w1 * w2

            # Just a safety condition.
            if ( math.isnan(varcoeff_coverage_peaks_dict[i]) ):
                print(varcoeff_coverage_peaks_dict[i])

            filtered_num_peaks += 1
        else:
            varcoeff_coverage_peaks_dict[i] = -0.01

    #######################
    ##   Generate Output ##
    #######################

    print("[NOTE] {} peaks are covered.".format(filtered_num_peaks))

    # Filter our peaks that are completly uncovered.
    filtered_varcoeff_coverage_peaks = []
    for i in varcoeff_coverage_peaks_dict.values():
        if (i >= 0):
            filtered_varcoeff_coverage_peaks.append(i)

    print("[NOTE] Generate CV plot.")

    scale_max = numpy.max(filtered_varcoeff_coverage_peaks)
    if ( args.scale_max ):
        if ( float(args.scale_max) < 0.0 ):
            sys.exit("[ERROR] Wrong value for scale_max!")
        scale_max = float(args.scale_max)

    # Make vase plot of variationkoefficients.
    f1 = plt.figure(figsize=(3, 5), dpi=350)
    plt.violinplot(filtered_varcoeff_coverage_peaks)
    plt.ylim(0.0, scale_max)
    plt.xticks([], 'none')
    plt.ylabel('Coefficient of Variation of the Peak Profiles')
    f1.savefig(args.output_folder + "/CV_Distribution_{}.pdf".format(outfilename), bbox_inches='tight')

    # Normalize all VC so that the scale goes from 0 to 1 which makes a comparison between
    # different profile evaluations easier. Unity-based normalization.
    one = numpy.max(filtered_varcoeff_coverage_peaks)
    zero = numpy.min(filtered_varcoeff_coverage_peaks)
    if ( args.max_norm_value ):
        one = float(args.max_norm_value)

    for i in range(0, len(filtered_varcoeff_coverage_peaks)):
        filtered_varcoeff_coverage_peaks[i] = (filtered_varcoeff_coverage_peaks[i]-zero)/(one-zero)

    print("[NOTE] Generate normalized CV plot.")

    # Make vase plot of variationkoefficients.
    f2 = plt.figure(figsize=(3, 5), dpi=350)
    plt.violinplot(filtered_varcoeff_coverage_peaks)
    plt.ylim(0.0, 1.0)
    plt.xticks([], 'none')
    plt.ylabel('Normalized Coefficient of Variation of the Peak Profiles')
    f2.savefig(args.output_folder + "/Norm_CV_Distribution_{}.pdf".format(outfilename), bbox_inches='tight')

    # Generate the output tabular file.
    print("[NOTE] Generate output tabular.")
    out_tab_file = open(args.output_folder + "/CV_tab_{}.bed".format(outfilename), "w")

    index_sort = numpy.argsort(list(varcoeff_coverage_peaks_dict.values()))[::-1]
    keys_list = list(varcoeff_coverage_peaks_dict.keys())

    for i in index_sort:
        k = keys_list[i]
        coords = coordinates_dict[k]

        max_norma_cv = (varcoeff_coverage_peaks_dict[k]-zero)/(one-zero)

        type = 0
        if ( max_norma_cv >= float(args.thresh)):
            type = 0
        else:
            type = 1

        # Chr  Start   End
        # ID_from_File    CV  Strand
        # bp    r   p
        # Max_Norm_CV
        # Left_Border_Center_Difference
        # Right_Border_Center_Difference Inter_Index
        # Internal_Peak_ID  Specific/Unspecific"
        out_tab_file.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\t{12}\t{13}\n".format(coords[0], coords[1], coords[2],
                                                                 peakid_dict[k], varcoeff_coverage_peaks_dict[k], strand_dict[k],
                                                                 num_bp_peaks_dict[k], size_r_peaks_dict[k], prob_success_peaks_dict[k],
                                                                 max_norma_cv,
                                                                 center_border_diff_left_dict[k],
                                                                 center_border_diff_right_dict[k],
                                                                 k+1, type))
    out_tab_file.close()

    ### write data for clustering
    numpy.savetxt(args.output_folder + "/data_classification_{}.tsv".format(outfilename), cov_matrix_for_R, delimiter="\t", newline="\n")

    script_path = os.path.join(os.path.dirname(os.path.dirname(fnb.__file__)), "lib/uMAP.R")

    ### Create R call
    Rcall = "Rscript {0} {1} {2} {3} {4} {5}".format(script_path, args.output_folder, outfilename, args.lam, args.maxcl, args.numcl)

    if ( args.sm ):
        print("[NOTE] Smoothing turned on")
        Rcall += " TRUE"
    else:
        Rcall += " FALSE"

    ### Run Rscript for classification
    if ( args.turn_off_classification ):
        print("[NOTE] Skip Classification")
    else:
        print("[NOTE] Run Classification")
        sb.Popen(Rcall, shell=True).wait()

    print("[FINISH]")

if __name__ == '__main__':
    main()
