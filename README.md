This is the code for my BIOL 499 Project. It includes scripts, raw data and extra csvs. 
---
* fs = field survey
* gc = growth chamber
---
**I only notated when I thought it was neccessary**
---
Here is the table of contents: 

`1_scripts/`
- [biol499script1_gccleandata.R](./1_scripts/biol499script1_gccleandata.R) Cleaning germination data
- [biol499script2_cleanseedtraitforuse.R](./1_scripts/biol499script2_cleanseedtraitforuse.R) Cleaning seed trait data
- [biol499script3_v2_traitanalysis.R](./1_scripts/biol499script3_v2_traitanalysis.R) Statistical analysis of plant traits
- [biol499script3_v2_traitvisualization.R](./1_scripts/biol499script3_v2_traitvisualization.R) Generation of trait-related plots and figures
- [biol499script4a_v2_glmm_indices.R](./1_scripts/biol499script4a_v2_glmm_indices.R) GLMM models for diversity indices - Abundance, Richness, Evenness
- [biol499script4a_visualizations.R](./1_scripts/biol499script4a_visualizations.R) Visualizing GLMM models 
- *[biol499script4b_v2_multivariate.R](./1_scripts/biol499script4b_v2_multivariate.R) Multivariate analysis (NMDS/PERMANOVA)*
- [biol499script4c_abovebelowground.R](./1_scripts/biol499script4c_abovebelowground.R) Above-ground vs. below-ground comparison
- [biol499script4c_visualization.R](./1_scripts/biol499script4c_visualization.R) Visualization for above/below-ground data
- [biol499scripte1_TRYdataextraction.R](./1_scripts/biol499scripte1_TRYdataextraction.R) Extraction of trait data from the TRY database
- [biol499scripte2_mergingseed.R](./1_scripts/biol499scripte2_mergingseed.R) Merging seed datasets for final analysis
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
*Last updated: March 10, 2026*
