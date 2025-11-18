## Hybrid untargeted and targeted RNA sequencing facilitates genotype-phenotype associations at  single-cell resolution

### Setup 
This workflow was implemented using **Snakemake v8.25.1**, **Python v3.12.10**, **conda v25.5.1**, and **R v4.5.1**.  
To run the workflow, modify the R version and data directory in the `config.yaml` file.
You can then execute the Snakemake workflow with:

```bash
snakemake --cores 8 --use-conda
```

### Figures
All figures in the manuscript were generated using files located in the `docs` directory.

### Data Availability
Expression count table is available at https://zenodo.org/records/17494517. The raw sequencing data is under the submission process under EGA. 

### Citation
```bash
@ARTICLE{Wang2025-mv,
  title    = "Hybrid untargeted and targeted {RNA} sequencing facilitates
              genotype-phenotype associations at single-cell resolution",
  author   = "Wang, Jiayi and Maldifassi, Maria Constanza and
              Bratus-Neuenschwander, Anna and Zhang, Qin and Beuschlein, Felix
              and Penton, David and Robinson, Mark D",
  journal  = "bioRxiv",
  month    =  nov,
  year     =  2025
}

```

### Contact
If you have any question regarding the code, please contact jiayi.wang2@uzh.ch
