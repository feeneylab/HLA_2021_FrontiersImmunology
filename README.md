# HLA_2021_FrontiersImmunology
HLA alleles B*53:01 and C*06:02 are associated with higher risk of P. falciparum parasitemia in a cohort in Uganda

Stata do files:
Malaria_HLA_1_Variables_2020-10-31_git.do
  Creates all variables necessary for analysis. Creates both .dta files.
Malaria_HLA_2_Analysis_2020-10-31_git.do
  Primary analysis file.
Malaria_HLA_3_Analysis_FDR_2020-12-21_git.do
  Merges primary analysis and batch permutation results. Calculates false discovery rate q-values and creates tables.
BatchPermute_allmodels_2020-07-07_git.do
  Permutation analysis file.

Stata dta files:
Malaria_HLA_Variables_2021-02-26_git.dta
  Primary .dta file.
Malaria_HLA_Variables_PSI_2021-02-26_git.dta
  Above file subsetted for probability of symptoms if infected analyses.
