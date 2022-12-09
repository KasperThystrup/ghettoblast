# ghettoblast

## Dependencies
This workflow have been developped with the following versions, previous versions have not been tested
* Snakemake 6.13.1
* Python 3.6.5
* conda 4.6.14

## Setup
In order to make the pipeline work, fill following text field, copy it into a plain text file and save it as `config/config.yaml`
```
query_file:
        /path/to/target_seqs.fasta
sample_path:
        /path/to/sample_dir
outdir:
        /path/to/output_dir
exclude_seqs:
        FALSE
top_only:
        TRUE

```
When exclude_seqs is TRUE, query and subject sequences will be extracted from the final blast matrix, and the subject sequences saved as individual fasta files.
When top_only is TRUE, the final blast matrix only contains tophits for each subject sequences against each query sequences.

The workflow automatically scans the `sample_dir` for `.fasta` files and will attempt to include these in the workflow, using their file name as **sample names** (excluding the  `.fasta` extension. 

Example:

```
> ls /path/to/sample_dir/
S_aureus1.fasta  S_aureus2.fasta  C_freundii1.fasta
```

Their **sample name** will be `S_aureus1` `S_aureus2`, and `C_freundii1`

### Sample exclusion
If you wish to exclude specific samples within the `sample_dir/`, provide sample names into the `config/blacklist.tsv` file. Make sure to maintain `sample` as the header, and write one **sample name**  per line (note don't use their full file names!).
Example: to exclude `C_freundii1` the `config/blacklist.tsv` file should look like this

```
sample
C_freundii1
```

# Execution
To execute the workflow first activate the snakemake environment, then specify the amounts of cores to use (for parallization) with `--cores X` and add the `--use-conda` so that the correct environments can be generated and activated. 
```
> cd ghettoblast
> conda activate snakemake
> snakemake --cores 3 --use-conda
```
