# Download reference `Neisseria gonorrhoeae` reference using NCBI Datasets

```
datasets download genome accession GCA_013030075.1 --include gff3,rna,cds,protein,genome,seq-report
```

# Download NCBI SRA data using `download_sra.sh` script

```
./download_sra.sh -s "ERR11171225" > "ERR11171225".log 2>&1 &
```

# Get AMRGen (Euro-GASP) using R script `amrgen-download-pheno.R`

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

# Align short-reads to NG reference and get VCF file 

1. Index the reference

```
kma index \
-i Neisseria-gonorrhoeae-TUM19854.fasta \
-o index_kma_db > index.log 2>&1 &
```

2. Create jobs

```
while read accesion; do
    if [ -e "$file" ]; then
        R1=${accesion}_R1.fastq.gz
        R2=${accesion}_R2.fastq.gz
        echo "kma -ipe ${R1} ${R2} -o ${accesion}_output -t_db index_kma_db -t 4 -nc -a -vcf > $accesion.kma.log 2>&1 &"
    fi
done < accessions.txt
```

3. Align reads to index database


```
kma -ipe ERR11174636_R1.fastq.gz ERR11174636_R2.fastq.gz -o ERR11174636_output -t_db index_kma_db -t 4 -nc -a -vcf > ERR11174636.kma.log 2>&1 &
```

