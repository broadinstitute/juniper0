# JUNIPER

Specht et al. 2025, https://doi.org/10.1101/2025.03.02.25323192

JUNIPER: Large-Scale Outbreak Reconstruction

## Installation

### Installation of R package

Installing `juniper0` currently requires installation through GitHub. To
do this, first you will need to install (if you have not already) and
load the package `devtools` in R:

```         
install.packages("devtools")
library(devtools)
```

After loading `devtools`, the `juniper0` package can be installed and
attached to the current R script by running:

```         
install_github("broadinstitute/juniper0")
library(juniper0)
```
The typical installation time for the `juniper0` is under one minute. This does not include the installation time of dependencies (listed below). It has been tested on Mac OS and Linux.

### Dependencies

`juniper0` requires `igraph` (version >= 2.1.4), `ape` (version >= 5.8.1), `ggplot2` (version >= 3.5.1), `ggraph` (version >= 2.2.1), `extraDistr` (version >= 1.10.0), and `TransPhylo` (version >= 1.4.5) as dependencies.

## Assembling and formatting data

After you have installed `juniper0`, ensure that the data are in the
correct format and stored in the file structure that `juniper0`
requires. The only *required* inputs are a metadata table and aligned sequences (either as one combined file or with each genome in a separate file). It is
*strongly recommended* that iSNV data in the form of Variant Call Format
(VCF) files also be included when available.

### Metadata table

The metadata table must contain a column named `sample` containing the 
sequence names to be reconstructed and a column named `date` containing the 
date of sample collection, in `YYYY-MM-DD` format. The names provided in the 
`sample` column must be substrings the sequence names if using a single aligned FASTA, 
or substrings of the file names if providing one FASTA per sequence. If VCF files are provided,
the sequence names in the metadata table must be substrings of the corresponding VCF file names.

### Aligned FASTA

For consensus-level genomes, `juniper0` accepts either a single FASTA file which must be named `aligned.fasta`, or 
one FASTA file per genome. As mentioned above, each sequence name in the metadata table must appear as either a 
sequence name in `aligned.fasta` or as a file name for one of the per-individual FASTAs. It is critical that each of the sequences in this
`aligned.fasta`, or each individual FASTA file, has been aligned to the same reference sequence.

By default, all files should be stored in a user-created directory
named `input_data` which is inside of your current working directory in
R. For instance the following R code would set the working directory and
create a sub-folder with the required name for a Mac user with profile
named `myname`:

```         
setwd("/Users/myname/Desktop/projects/")
mkdir("/Users/myname/Desktop/projects/input_data")
```

### VCF files

Note that `juniper0` does not require iSNV data and additionally does
not require iSNV data from every sequence. It is generally recommended,
however, to include this data whenever available. To include iSNV data
in the reconstruction process, you must include a Variant Call Format
(VCF) file per sequence.

To run `juniper0` with iSNV data, place all available VCF files in the same directory as the metadata file and FASTA file(s).
The name of each VCF file **must contain the name of the consensus sequence
to which is corresponds**. As long as the sequence name (ex: `SEQ001a`)
is fully present in the VCF name (even as a sub-string, ex:
`cleaned_final_SEQ0001a_lofreq.vcf`) then `juniper0` will be able to
match the VCF to the consensus sequence.

## Running `juniper0`

### Running with default parameters

Once your data is in the format shown above and you have created an R
script in the directory that also contains your `input_data` folder,
`juniper0` can be run using default settings as follows:

```         
init <- initialize()
results <- run_mcmc(init)
```

### Running with custom parameters

If you prefer to give your input directory a non-standard name or file
path (ex: `my_path/my_dir`) then use the following command (**note
that** **the sub-directory `VCF` must still be present in this folder**
if you intend to use iSNV data in the reconstruction):

```         
init <- initialize(indir="my_path/my_dir")
```

While `juniper0` provides joint estimation of many epidemiological
parameters with limited prior information, results will generally
improve when pathogen-specific priors are provided. Estimates of the
generation interval and mutations rate (usually readily available in
epidemiological literature) can be provided. For instance, when
reconstructing Hepatitis A Virus, one might set:

```         
init <- initialize(a_g = 26.9,
                   init_mu = 1.99e-4)
```

Run `?initialize` to see more details about optional parameters.

### Viewing the results

Once you have run the above with an appropriate number of iterations and
filters, you can view the summary of results by running the command:

```         
res_summary <- summarize(results)
res_summary
```

## Examples

For detailed instructions on reproducing the JUNIPER runs in the manuscript, see the repository https://github.com/broadinstitute/juniper-analyses.

## Troubleshooting

Please report any bugs to the lead and corresponding author, Ivan Specht, at `ispecht@broadinstitute.org`.

