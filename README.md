# LuscombeU_annotationpolish

## Description
R package for annotation polishing. This package can perform the following on base gff3 files:
1. Add 5' UTR based on CAGE bed files (`add_five_prime_utr()`)
2. Add 3' UTR based on canonical polyA site using information from BSgenome (`add_three_prime_utr_wrapper()`)
3. Clean manual gff3 annotation written by hand (`clean_manual_anno()`)
4. Add unique identifier column `ID` (`add_unique_id()`)
5. Create a list of operons from gene lists and CAGE non-SL peaks (`make_operons()`)
6. Fix liftoff gff files previously fixed wrongly by AGAT (`fix_liftoff_oikobase()`)
7. Prioritize one transcript model over another (`pick_transcripts()`)
8. Split chimeric transcripts based on CAGE and liftoff (`split_transcripts()`)


## Installation
```r
devtools::install_github("oist/LuscombeU_annotationpolish")
```



