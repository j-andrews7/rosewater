# rosewater

**rosewater** assigns genes to (super-) enhancer output from [ROSE](https://bitbucket.org/young_computation/rose/src/master/) in an expression-aware manner. It allows users to set a TPM threshold to filter genes that are not expressed on a sample-by-sample basis.

## Installation

`rosewater` can be installed via pip. For use, it **requires `bedtools` be available on your PATH**.

```pip install rosewater```

## Usage

`rosewater` is fairly simple to use. It requires an annotation GTF file, a TSV file of TPMs with the gene name column named 'gene' (these should match the 'gene_name' attributes in the annotation GTF), the name of the sample column in the TPM file, and an output file from [ROSE](https://bitbucket.org/young_computation/rose/src/master/). Optionally, users can set a TPM threshold (set to 5 by default) for filtering out lowly/non-expressed genes prior to assignment.

```
Usage: rosewater [OPTIONS]

  rosewater assigns genes to ROSE output in an expression-aware manner.

Options:
  -a, --annotation PATH       Path to annotation GTF file.  [required]
  -s, --sample TEXT           Sample name that should match a column in the
                              TPM file.  [required]

  -e, --enh_file PATH         Output from ROSE ending with
                              'ENHANCER_TO_GENE.txt'.  [required]

  -t, --tpm_file PATH         A file containing a matrix of TPMs with genes as
                              rows and samples as columns. The gene label
                              column should be named 'gene'.  [required]

  -th, --tpm_threshold FLOAT  The minimum TPM to retain genes for assignment.
                              [default: 5]

  -o, --output_dir PATH       The output directory.  [default:
                              ./EnhancerGeneAssignments]

  --version                   Show the version and exit.
  -h, --help                  Show this message and exit.

```

## Output

Two output files will be generated, named after the ROSE enhancer input file appended with either `.rosewater.GeneAssignment.log` or `.rosewater.GeneAssignment.bed`. The **log file** will contain useful stats such as how many TSSes are filtered by the TPM threshold, how many TSSes overlap each enhancer, the average enhancer size, and how many assignments change from the original ROSE assignments. 

The **BED-like file** will contain the assignments for each enhancer. Two assignments are made for each enhancer - one utilizing all TSSes in the annotation file that meet the TPM threshold and another utilizing only the protein-coding TSSes. These assignments are the last 4 columns of the file. The additional columns are fairly self-explanatory. In short, they contain the overlapping TSSes, the closest TSS using all transcripts that meet the TPM threshold, the closest TSS using only protein-coding transcripts that meet the TPM threshold, and the TPMs for each of those.

## Assignment Logic

The original ROSE gene mapper just assigns the TSS that is closest to the center of the enhancer. `rosewater` takes a more sophisticated (and therefore complicated approach):

- If the enhancer overlaps no TSSes for a gene that meets the TPM threshold:
	- The "final_allgene_assignment" will be set to the gene that meets the TPM threshold for the closest TSS while "final_proteincoding_assignment" will be set to the gene that meets the TPM threshold for the closest 'protein_coding' TSS.
- If the enhancer overlaps a single TSS for a gene that meets the TPM threshold:
	- If the 'gene_type' of the gene is `protein_coding`, the "final_allgene_assignment" and "final_proteincoding_assignment" will both be set to that gene.
	- If the 'gene_type' of the gene is **not** `protein_coding`, the "final_allgene_assignment" will be set to that gene while the "final_proteincoding_assignment" will be set to the gene for the closest 'protein_coding' TSS.
- If the enhancer overlaps two or more TSS for a gene that meets the TPM threshold:
	- If the 'gene_type' of the most highly-expressed gene is `protein_coding`, the "final_allgene_assignment" and "final_proteincoding_assignment" will both be set to that gene.
	- If the 'gene_type' of the most highly-expressed gene is **not** `protein_coding`, the "final_allgene_assignment" will be set to that gene while the "final_proteincoding_assignment" will be set to the most highly-expressed overlapping 'protein_coding' gene. If there are no overlapping TSSes for 'protein_coding' genes, the "final_proteincoding_assignment" will be set to the gene for the closest 'protein_coding' TSS.

Users are free to use whichever assignment they feel is most appropriate for their use case.

## Known Issues

Users may get a warning like `RuntimeWarning: line buffering (buffering=1) isn't supported in binary mode` depending on their version of python. This can be safely ignored. It stems from `pybedtools` and should be fixed in their next release.

## Contributing

Feel free to submit a [pull request](https://github.com/j-andrews7/rosewater/pulls) or open an [issue](https://github.com/j-andrews7/rosewater/issues) if you have an idea to enhance this tool.

## License

`rosewater` is available under the [GNU-GPLv3 license](https://github.com/j-andrews7/rosewater/blob/master/LICENSE). It is provided as-is, with no warranty or guarantees. 