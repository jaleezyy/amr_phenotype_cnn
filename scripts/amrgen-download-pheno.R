# Load libraries

library(dplyr)
library(tidyr)
library(ggplot2)

# ack - Kara
# set path
setwd("/Users/amos/Documents/repos/amr_phenotype_cnn/scripts")

# install.packages("AMRgen")
library(AMRgen)
if (!requireNamespace("AMRgen", quietly = TRUE)) {
  devtools::load_all()
}

# this works, for installing AMRgen
#install.packages("remotes") # if you haven't already
#remotes::install_github("AMRverse/AMRgen")


# reference
# https://amrgen.org/articles/NeisseriaGonoExamples.html#importing-genotype-and-phenotype-data-into-amrgen

eurogasp_pheno <- eurogasp_pheno_raw %>%
  mutate(across(c(Azithromycin, Ciprofloxacin, Cefixime, Ceftriaxone), as.character)) %>%
  pivot_longer(
    cols = c(Azithromycin, Ciprofloxacin, Cefixime, Ceftriaxone),
    names_to = "antibiotic",
    values_to = "mic"
  )

eurogasp_ast <- format_ast(
  input = eurogasp_pheno,
  sample_col = "id",
  species = "Neisseria gonorrhoeae",
  ab_col = "antibiotic",
  mic_col = "mic",
  interpret_eucast = TRUE,
  interpret_ecoff = TRUE
)

res <- eurogasp_ast %>% filter(!is.na(pheno_eucast)) %>% filter(pheno_eucast != 'I') %>%
  filter(drug_agent == 'CIP') %>%
  select(id,drug_agent,pheno_eucast)

res <- res[order(res$id),]

write.table(res,file = "/Users/amos/Documents/repos/amr_phenotype_cnn/data/NeisseriaGono-CIP.tsv", sep="\t")
