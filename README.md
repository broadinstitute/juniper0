---
editor_options: 
  markdown: 
    wrap: 72
---

# JUNIPER

Specht et al. 2025

Preliminary Version of JUNIPER: Large-Scale Outbreak Reconstruction

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

## Assembling and formatting data

After you have installed `juniper0`, ensure that the data are in the
correct format and stored in the file structure that `juniper0`
requires. The only *required* input is an aligned FASTA file which
contains the names and dates of sequences to be reconstructed. It is
*strongly recommended* that iSNV data in the form of Variant Call Format
(VCF) files also be included when available. There is also optional
functionality to include contact tracing data when available.

### Aligned FASTA

The primary input that `juniper0` requires is a FASTA alignment of all
the sequences you intend to reconstruct, **which must be named
`aligned.fasta`**. It is critical that each of the sequences in this
FASTA has been aligned to the same reference sequence and has been named
as `>sequence_name|YYYY-MM-DD.` This date represents the date that will
be used to reconstruct transmission between sequences; typically either
the specimen collection date or the date of symptom onset, depending on
the pathogen and availability of data.

By default, this FASTA should be stored in a user-created directory
named `input_data` which is inside of your current working directory in
R. For instance the following R code would set the working directory and
create a sub-folder with the required name for a Mac user with profile
named `myname`:

```         
setwd("/Users/myname/Desktop/projects/")
mkdir("/Users/myname/Desktop/projects/input_data")
```

See [Running with custom parameters] for more details on file paths. If
you do not have an aligned FASTA in this format, see [Other
considerations] for help creating this file.

### VCF files

Note that `juniper0` does not require iSNV data and additionally does
not require iSNV data from every sequence. It is generally recommended,
however, to include this data whenever available. To include iSNV data
in the reconstruction process, you must include a Variant Call Format
(VCF) file per sequence.

To run `juniper0` with iSNV data, create a sub-folder called `vcf`
within the `input_data` folder. The `vcf` sub-folder should contain VCF
files for any sequence included in `aligned.fasta` that you wish to
reconstruct using iSNV data.

The name of each file **must contain the name of the consensus sequence
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

See the documentation of `initialize()` for more details optional
parameters.

### Viewing the results

Once you have run the above with an appropriate number of iterations and
filters, you can view the summary of results by running the command:

```         
res_summary <- summarize(results)
res_summary
```

