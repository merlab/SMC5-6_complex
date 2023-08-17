# SMC5/6 complex

Code associated with the bioinformatics analysis presented in the manuscript:

**Large-scale phenogenomic analysis of human cancers uncovers frequent alterations affecting SMC5/6 complex components in breast cancer**

Shamayita Roy<sup>1</sup>, Arvin Zaker<sup>2</sup>, Arvind Mer<sup>2,#</sup> & Damien D’Amours<sup>1,#</sup>

1 Ottawa Institute of Systems Biology, Department of Cellular and Molecular Medicine, University of Ottawa, Canada

2 Department of Biochemistry, Microbiology & Immunology, University of Ottawa, Canada

\#Corresponding author

**ABSTRACT**

> *Cancer cells often experience large-scale alterations in genome architecture because of DNA damage and replication stress. Whether mutations in core regulators of chromosome structure can also lead to cancer-promoting loss in genome stability is not fully understood. To address this question, we conducted a systematic analysis of mutations affecting a global regulator of chromosome biology –the SMC5/6 complex– in cancer genomics cohorts. Analysis of 64,959 cancer samples spanning 144 tissue types and 199 different cancer genome studies revealed that the SMC5/6 complex is frequently altered in breast cancer patients. Patient-derived mutations targeting this complex associate with strong phenotypic outcomes such as loss of ploidy control and reduced overall survival. Remarkably, the phenotypic impact of several patient mutations can be observed in a heterozygous context, hence providing an explanation for a prominent role of SMC5/6 mutations in breast cancer pathogenesis. Overall, our findings suggest that genes encoding global effectors of chromosome architecture can act as key contributors to cancer development in humans.*


## Steps to reproduce the results:

1. Download the manuscript associated data files from Zenodo at https://doi.org/10.5281/zenodo.8256823
2. Place the files in the data folder in the data folder of the repository.
3. Runs the scripts in the `./R/` folder from `s01_...` to `s12_...` in order. 
4. Run the script pretaining to the figure you would like.

## Example

To reproduce figure 2, run the `fig2.R` script. Ensure that you are executing the script.
From the root directory of the repository (from the `./SMC5-6_complex` folder).

```
Rscript ./R/fig2.R
```

Results:

![](readme-img.png)
