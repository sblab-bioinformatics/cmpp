Authors: Kamal Kishore (CRUK-CI) and Sergio Martinez Cuesta (CRUK-CI and Deparment of Chemistry, Unversity of Cambridge)

## Filtering and imputation

```r
library(qPLEXanalyzer)
library(gridExtra)
library(pander)
library(readxl)
library(data.table)

data(human_anno)

metadata <- read_excel("../data/labelfree.xlsx", sheet = "metadata")
intensities <- read_excel("../data/labelfree.xlsx", sheet = "peptide_intensities")

MSnset_data <- convertToMSnset(intensities, metadata=metadata, indExpData=c(7:30), Sequences=2, Accessions=6, rmMissing = FALSE)
```

Continue incorporating `PR1239.rmd` and `20200526_proteomics.md`
