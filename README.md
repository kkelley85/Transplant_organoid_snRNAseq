# Maturation and circuit integration of transplanted human cortical organoids
Repository for code used to analyze single-nuclei RNA sequencing data from transplanted human brain cortical organoids for accompanying paper:

<p align="center">
<img src="/Fig1_snRNAseq.png" width="600"/>
</p>

<p>Revah R<sup>*</sup>, Gore F<sup>*</sup>, Kelley KW<sup>*</sup>, Andersen J, Sakai N, Chen X, Li MY, Birey F, Yang X, Saw NL, Baker SW, Amin ND, Kulkarni S, Mudipalli R, Cui B, Nishino S, Grant GA, Knowles JK, Shamloo M, Huguenard JR, Deisseroth K, Pașca SP. Maturation and circuit-level integration of transplanted human cortical organoids. Nature (2022)<p>
  
**_Nature:_** https://www.nature.com/
## Data links
Raw and processed data available through GEO: [GSE190815](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE190815)
| GEO_ID        | Condition     | Diff. days   | Line |
| ------------- | ------------- | ------------ | ---- |
| [GSM5732392](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM5732392)    | t-hCS | 276 | 8119-1 |
| [GSM5732393](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM5732393)    | t-hCS | 224 | 2242-1 |
| [GSM6225773](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM6225773)    | t-hCS | 227 | 2242-1 |
| [GSM5732394](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM5732394)    |  hCS  | 243 | 1208-2 |
| [GSM6225774](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM6225774)    |  hCS  | 257 | 8119-1 |
| [GSM6225775](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM6225775)    |  hCS  | 224 | 2242-1 |

  
Processed and integrated Seurat binarized R objects in **.rds** format: 
* [transplanted hCS .rds file 2GB](https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE190815&format=file&file=GSE190815%5Ft%2DhCS%5Fprocessed%5FSeuratObject%2Erds%2Egz)
* [hCS .rds file 2GB](https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE190815&format=file&file=GSE190815%5FhCS%5Fprocessed%5FSeuratObject%2Erds%2Egz)
## Preprocessing
Combined human and rat genome references were created usng 10x Cellranger software [suite](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/what-is-cell-ranger). Reference notes can be found here: [Combined_HumanRat_reference.sh](code/Combined_HumanRat_reference.sh)
  
Count matrices were created using Cellranger count v6.1.2 and processed on the 10x cloud analysis [server](https://www.10xgenomics.com/products/cloud-analysis). Following paramters were used:`include_introns: true`
  
To identify putative high-quality human nuclei from transplanted cortical organoids the following processing code was used: [Preprocess.R](code/Preprocess.R)
 
## Analysis scrips
Code to generate QC plots from Extended Data Figure 4: [QC_plots.R](code/QC_plots.R)
  
Code to generate t-hCS plots in Figure 1 and Extended Data Figure 4: [t-hCS_plots.R](code/t-hCS_plots.R)
  
Code to generate hCS plots in Extended Data Figure 4: [hCS_plots.R](code/hCS_plots.R)

Code to perform pseudobulk differential expression analysis using [Libra package](https://github.com/neurorestore/Libra) for Supplementary Table 5: [GluN_DE_analysis.R](code/GluN_DE_analysis.R)

Code to perform gene set significance testing and heatmap plots from activity-dependent gene sets displayed in Figure 2: [ActivityGenePlots.R](code/ActivityGenePlots.R)

Code to perform cluster overlap mapping to reference datasets as shown in Extended Data Figure 5: [Dataset_clusterMapping.R](code/Dataset_clusterMapping.R)

Code to perfrom glutamatergic neuron transfer labeling from adult Allen brain references as shown in Extended Data Figure 5: [GluN_Allen_mapping.R](code/GluN_Allen_mapping.R)

Example outputs from above code: [output](output)
