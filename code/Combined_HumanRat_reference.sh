#Example code and paramaters used to make combined human rat genome reference using 10x CellRanger software

#Human reference assembled following 10x instructions: https://support.10xgenomics.com/single-cell-gene-expression/software/release-notes/build#GRCh38_2020A
#source files:
#fasta_url="http://ftp.ensembl.org/pub/release-98/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
#gtf_url="http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/gencode.v32.primary_assembly.annotation.gtf.gz"

#Rat reference 
#source files: ftp://ftp.ensembl.org/pub/release-100/fasta/rattus_norvegicus/dna/Rattus_norvegicus.Rnor_6.0.dna.toplevel.fa.gz
#software version: cellranger-3.1.0

cellranger mkgtf Rattus_norvegicus.Rnor_6.0.100.chr.gtf Rattus_norvegicus.Rnor_6.0.100.chr.filtered.v2.gtf \
                 --attribute=gene_biotype:protein_coding \
                 --attribute=gene_biotype:lincRNA \
                 --attribute=gene_biotype:antisense \
                 --attribute=gene_biotype:IG_LV_gene \
                 --attribute=gene_biotype:IG_V_gene \
                 --attribute=gene_biotype:IG_V_pseudogene \
                 --attribute=gene_biotype:IG_D_gene \
                 --attribute=gene_biotype:IG_J_gene \
                 --attribute=gene_biotype:IG_J_pseudogene \
                 --attribute=gene_biotype:IG_C_gene \
                 --attribute=gene_biotype:IG_C_pseudogene \
                 --attribute=gene_biotype:TR_V_gene \
                 --attribute=gene_biotype:TR_V_pseudogene \
                 --attribute=gene_biotype:TR_D_gene \
                 --attribute=gene_biotype:TR_J_gene \
                 --attribute=gene_biotype:TR_J_pseudogene \
                 --attribute=gene_biotype:TR_C_gene

#######################################
#Combined human rat reference
#software version: cellranger-5.0.0

cellranger mkref --ref-version=v5 \
    --genome=GRCh38_2020A_v5 --fasta=Homo_sapiens.GRCh38.dna.primary_assembly.fa.modified --genes=gencode.v32.primary_assembly.annotation.gtf.filtered \
    --genome=Rnor_6.0 --fasta=Rattus_norvegicus.Rnor_6.0.dna.toplevel.fa --genes=Rattus_norvegicus.Rnor_6.0.100.chr.filtered.v2.gtf \
    --nthreads=30 \
    --memgb=190
