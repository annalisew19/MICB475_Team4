# March 17 2025

### What we have done

##### Aim 1 (Annalise):
Ran Beta diversity metrics
- Weighted Unifrac:
  
  <img src="../R_project/aim1/weighted_unifrac_pcoa_try2.png" height="300" width="400">
  
  Statistical Test:

  <img src="../R_project/aim1/wunifrac_stat.png" height="100" width="300"> 

- Bray-Curtis:

  <img src="../R_project/aim1/bray_curtis_pcoa_try1.png" height="300" width="400">

  Statistcial Test:

  <img src="../R_project/aim1/braycurt_stat.png" height="100" width="300">

Taxonomic Composition

<img src="../R_project/aim1/tax_composition.png" height="600" width="350">

Alpha diveristy metrics 
- Shannon's Diversity:

<img src ="../R_project/aim1/shannon_diversity.png" height="300" width="400">

- Faith's PD:

<img src="../R_project/aim1/faithpd_boxplot.png" height="300" width="400">

Still have to do:
- Linear regression model
- Taxonomic composition statistical analysis

##### Aim 2 (Michelle):
Core microbiome based on outcome (detecion= 0.001, prevalence = 0.1): only 1 ASV that is unique to unsuccessful core microbiome, and 6 ASVs unique to successful core microbiome. 
- <img src="../images/core_microbiome_outcome.png" height="300" width="400">
- <img src="../images/core_mic_successful.png" height="130" width="800">
- <img src="../images/core_mic_unsuccessful.png" height="50" width="800">
Core microbiome based on outcome (detecion= 0.001, prevalence = 0.2)
- <img src="../images/core_microbiome_outcome2.png" height="300" width="400">
ISA analysis based on age group and outcome (p= 0.05)
- <img src="../images/ISA_age_outcome_p0.05.png" height="70" width="1200">
ISA analysis based on age group and outcome (p= 0.1)
- <img src="../images/ISA_age_outcome_p0.1.png" height="200" width="1200">

Still have to do:
- Statistical analysis

### Questions to ask/Issues
- CONTRADICTS PAPER FINDING: core microbiome analysis showed that Lactobacillus was the only unique ASV in the unsuccessful group, and streptomyces (along with 5 unassigned) is among the 6 unique ASVs in the successful group.
- 
### Meeting Notes
- linear regression plots can replace shannon diversity graphs
- for unassigned species, find the corresponding file (reps seq file) for blast
- complete result analysis by next monday, so will have time to complete additional analysis
- can start working on slides this week
- slides need to be submitted on March 30th (with a background and why we did this study), and have one day to prepare.
- have 10 min to ask the other team on the day of presentation, 10 min of presentation time, no speaker notes
- presentation details will be released soon

- taxonomic composition bar plots:
    - go down one more level in order to tell difference
    - don't need statistical test
- alpha diversity:
    - use linear regression code, one for shannon diversity, one for faith's pd
    - also add shannon to meta data like what was done with pd
- for core microbiome: 
    - can look into the literature 
    - can try to see if there are differences based on different outcomes (not just successful and unsuccessful)
      

### Next Week
- evelyn joins meeting next week
- need to present result next week during the meeting, have all three aims ready

