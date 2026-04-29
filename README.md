## Benchmarking variant calling and phasing methods using GIAB datasets.

### Setup 
This workflow was implemented using **Snakemake v8.25.1**, **Python v3.12.10**, **conda v25.5.1**, and **R v4.5.1**.  
To run the workflow, modify the R version and data directory in the `config.yaml` file. You also need to specify the 
directory for reference genome, transcriptome annotation files in `config.yaml`. 


You can then execute the Snakemake workflow with:

```bash
snakemake --cores 8 --use-conda
```

All GIAB raw sequencing files (BAM or FASTQ) are publicly available. Please refer the paper for the links. 
In addition, all figures can be generated from the Snakemake.
The output VCFs can be found on [Zenodo](https://zenodo.org/records/19857089).

### Contact
If you have any question regarding the code, please contact jiayi.wang2@uzh.ch.


