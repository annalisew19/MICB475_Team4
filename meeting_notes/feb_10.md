# February 10, 2025

## Agenda
- Discuss R code and the age binning we decided on
- Start planning team proposal and dividing work

### What we have done
R Portion:
- Loaded IVF dataset into R and made a histogram to show the number of samples for each age
<img src="../R_project/age_samplesizes/age_samplesizes.png" height="600" width="650">
- Decided to bin ages as follows: 25-30, 35-40, 40-45, 45-50. Will filter out the age 20-25 groups.
<img src="../R_project/age_samplesizes/agegroup_samplesizes.png" height="600" width="650">

Qiime2 Processing:
- Imported dataset and demultiplexed sequences. Sample size before demultiplexing = 824 samples
- Denoised sequences using DADA2. Set truncation length at 303 base pairs based on the interactive quality plot from demultiplexing sequences below:
<img src="../qiime2_files/qiime2view_screenshots/ivf_demux_qiime2_interactive_quality_plot.png" height="600" width="650">

### Questions to ask


## Meeting Notes


## Next week

