# StoatyDive

StoatyDive is a tool to evaluate predicted peak profiles to assess the binding specificity of a protein to its targets. It can be used for sequencing data such as CLIP-seq or ChIP-Seq.

## Installation

You can install StoatyDive via conda or with a direct download
of this repository.

### Conda

```
conda install stoatydive
```

### Repository download

```
git clone https://github.com/heylf/StoatyDive.git
cd StoatyDive
python setup.py install
python ./bin/StoatyDive.py -h
```

or

```
https://github.com/heylf/StoatyDive/archive/v1.0.1.tar.gz
```

Requirements: python >= 3.6, bedtools >= 2.27.1, numpy>=1.13.3, matplotlib>=2.1, scipy >= 0.19.1

## Usage

```
StoatyDive.py [-h] [options] -a *.bed -b *.bam/*bed -c *.txt
```

```
optional arguments:
  -h, --help            show this help message and exit
  -v, --version         show program's version number and exit
  -a *.bed, --input_bed *.bed
                        Path to the peak file in bed6 format.
  -b *.bam/*.bed, --input_bam *.bam/*.bed
                        Path to the read file used for the peak calling in bed
                        or bam format.
  -c *.txt, --chr_file *.txt
                        Path to the chromosome length file.
  -o path/, --output_folder path/
                        Write results to this path.
  --length_norm         Set length normalization. StoatyDive will expand every
                        peak to the maximal length.
  --length_norm_value int
                        Set length normalization value (maximum peak length).
  --max_norm_value float
                        Provide a maximum value for CV to make the normalized
                        CV plot more comparable.
  --border_penalty      Adds a penalty for non-centered peaks.
  --scale_max float     Provide a maximum value for the CV plot.
  --seed int            Set seed for the optimization scheme.
```

### Recommendations

#### Border Penalty and Length Normalization
It is recommended to use StoatyDive with `--border_penalty` and `--length_norm`.
Adding the border penalty takes care of peaks that are not correctly centered and might
just overlap with a short appendage of a read stack. The length normalization takes
care of different sized peaks. All peak are extended to a certain length with `--length_norm`. The user can either provide a peak length with `length_norm_value` or
StoatyDive just takes the maximal peak length of the given peak set.

#### Seed
Please set a seed if you use StoatyDive for any publications to warrant reproducible
results.

#### Other options
- You can set a maximal value for the normalized CV distribution plot with `max_norm_value`. This option helps, if you want to compare several  normalized
CV distribution plots from different experiments. Take the highest CV from all experiments as a maximal value.
- You can set a maximal value for the  CV distribution plot with `scale_max`. This option helps, if you want to compare several normalized CV distribution plots from different experiments. Take the highest CV from all experiments as a maximal value.

## Quick Example

Example 1:
`StoatyDive.py -a /test/broad_peaks/peaks.bed -b /test/broad_peaks/reads.bed -c /test/UCSC_edu_goldenPath_dm6_bigZips_dm6.chrom.sizes.txt --length_norm --border_penalty -o /test/broad_peaks/`

Example 2:
`StoatyDive.py -a /test/sharp_peaks/peaks.bed -b /test/sharp_peaks/reads.bed -c /test/UCSC_edu_goldenPath_dm6_bigZips_dm6.chrom.sizes.txt --length_norm --border_penalty -o /test/sharp_peaks/`

## Output

### CV distribution plot

| A        | B           |
| :-------------: |:-------------:|
| <img src="/test/broad_peaks/CV_Distribution_reads.png" width="400"> | <img src="/test/sharp_peaks/CV_Distribution_reads.png" width="400"> |

The diagram will give you a first impression of the binding specificity of your
protein of interest. The diagram also tells you about the performance/quality
of your experiment. An experiment with lots of unspecific binding sites will have
a CV distribution close to zero, as in our example A. An experiment with lots
of specific binding sites will have a CV distribution with a high expected CV, as
in our example B.

### Normalized CV distribution plot

| A        | B           |
| :-------------: |:-------------:|
| <img src="/test/broad_peaks/Norm_CV_Distribution_reads.png" width="400"> | <img src="/test/sharp_peaks/Norm_CV_Distribution_reads.png" width="400"> |

The normalized CV distribution helps to identify specific and unspecific sites within
an experiment. The normalized CV is in a range [0,1]. A specific site will have a value of 1. An unspecific site will have a value of 0.

## CV_tab file
The CV tabular file is a ranked, tab separated list of your predicted binding sites:
  1. Chromosome
  2. Start of Peak
  3. End of peak
  4. Peak ID/Name
  5. CV
  6. Strand
  7. Peak length
  8. r (hyperparameter of negative binomial)
  9. p (hyperparameter of negative binomial)
  10. Normalized CV
  11. Difference between the maximal value of the left border and the maximal
  value of the center region of the peak. (Penalty used for the border penalty.)
  12. Difference between the maximal value of the right border and the maximal
  value of the center region of the peak. (Penalty used for the border penalty.)
