This is the code for my BIOL 499 Project. It includes scripts, raw data and extra csvs. 
---
* fs = field survey
* gc = growth chamber
---
**I only notated when I thought it was neccessary**
---
Here is the table of contents: 

`1_scripts/`
- [biol499script1_clean_descript.qmd](./1_scripts/biol499script1_clean_descript.qmd) Cleaning germination data
- [biol499script2_cleanseedtraitforuse.R](./1_scripts/biol499script2_cleanseedtraitforuse.R) Cleaning seed trait data
- [biol499script3_traits.qmd](./1_scripts/biol499script3_traits.qmd) Statistical analysis of plant traits
- [biol499script4_univariate8.qmd](./1_scripts/biol499script4_univariate8.qmd) GLMM models for diversity indices - Abundance, Richness, Evenness
- *[biol499script6_multivariate.qmd](./1_scripts/biol499script6_multivariate.qmd) Multivariate analysis (NMDS/PERMANOVA)*
- [biol499script5_micrositeabovebelow.qmd](./1_scripts/biol499script5_micrositeabovebelow.qmd) Above-ground vs. below-ground comparisondata
- [spec_accum_curve.R](./1_scripts/spec_accum_curve.R) Species accumulation curve generation
  
`2_data_raw/`
- *[gc_seedlingdata.csv](./2_data_raw/gc_seedlingdata)*
- *[fs_summerveg.csv](./2_data_raw/fs_summerveg.csv)*

`4_outputfigures`

`5_outputtables`: output dataframes from above scripts
- *[gc_wcategories.csv](./5_outputtables/gc_wcategories.csv)* OUTPUT from biol499script1_gccleandata.R; used for scripts 3 and 4

`6_manualediting` : not raw data; manually edited excel
- [gc_seedtraitsfinal.csv](./6_manualediting/gc_seedtraitsfinal.csv) Generation of trait-related plots and figures
- [gcspeciestraits_edited.csv](./6_manualediting/gcspeciestraits_edited.csv) Generation of trait-related plots and figures
---
*Last updated: March 23, 2026*
