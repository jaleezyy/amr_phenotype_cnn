# Download reference `Neisseria gonorrhoeae` reference using NCBI Datasets

```
datasets download genome accession GCA_013030075.1 --include gff3,rna,cds,protein,genome,seq-report
```

# Get AMRGen (Euro-GASP) using R script `scripts/amrgen-download-pheno.R`

```
Rscript --vanilla amrgen-download-pheno.R
```

This saves a file called `NeisseriaGono-CIP.tsv`, which looks like so:

```
head ../data/NeisseriaGono-CIP.tsv 
        id      drug_agent      pheno_eucast
1       ECDC_BE18_497   CIP     S
2       ERR11171225     CIP     R
3       ERR11171451     CIP     R
4       ERR11171618     CIP     R
5       ERR11171778     CIP     S
6       ERR11171779     CIP     R
7       ERR11171780     CIP     S
8       ERR11171781     CIP     R
9       ERR11171782     CIP     S
```