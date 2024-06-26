# juniper0
Preliminary Version of JUNIPER: Large-Scale Outbreak Reconstruction

## Installation
### Installation of R package
Installing `juniper0` currently requires installation through GitHub directly. To do this, first you will need to install (if you have not already) and load the package `devtools` in R:
```
install.packages("devtools")
library(devtools)
```

After loading `devtools`, the `juniper0` package can be installed by running: 
```
install_github("broadinstitute/juniper0")
```

## Assembling and formatting data
After you have installed `juniper0`, you will need to assemble the data that you would like to run and ensure that all of the data are in the correct format. There is also a particular file structure that `juniper0` currently requires. It is important that your data is formatted exactly as shown in this guide and placed in folders exactly as shown. First, within the directory that you will be doing analysis, create a folder called `input_data`. This folder needs to include: 
 - `ref.fasta`
    - This should be the reference sequence that was used to align all other genomic samples that will be included. The actual name of the sequence in the FASTA file may be anything, but the file itself must be named exactly as above.
 - `aligned.fasta`
    - This should be an aligned FASTA file that contains all of the genomic samples you would like to reconstruct. 
 - `date.csv`
    - This should be a csv file where the first column contains the exact name of all genomic sequences in the `aligned.fasta` file and the second column contains the number of days between the collection date of each sample and the collection date of `ref.fasta`. 
    - For example if `aligned.fasta` contained three genomes named `genome_1`, `genome_2`, `genome_3` which were tested 1, 3, and 5 days after the collection date of the reference genome, respectively, then `date.csv` should look like:
      ```
      date = data.frame(col1=c("genome_1","genome_2","genome_3"), col2=c(1,3,5))
      print(date)
      ```
 - A folder called `vcf`
    - This folder should contain Variant Call Format (vcf) files for each of the genomic sequences included in `aligned.fasta`. The name of each file should match the exact name of the sequence. Each `.vcf` file should look like:
      ```
      file = read.table('demo.vcf')
      colnames(file) <- c('#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO')
      print(head(file))
      ```

Finally, in the directory that you created `input_data` create a new R script, which is where you will run `juniper0`.


## Running `juniper0`
Once your data is in the format shown above and you have created an R script in the directory that also contains your `input_data` folder, `juniper0` can be run using default settings as follows: 
```
init <- initialize()
results <- run_mcmc(init)
```

To run `juniper0` with custom settings, see the optional arguments available to 'initialize()' in the documentation.


Once you have run the above with an appropriate number of iterations and filters, you can view the results by running
```
summarize(results)
```
which returns matrices of direct and indirect transmissions and their associated posterior probabilities, and posterior samples of the global mutation rate (substitutitons/site/day), the within-host mutation rate (mutations/site/cycle), and the probability of an incomplete bottleneck.