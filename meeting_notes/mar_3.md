
# March 3, 2025

## Agenda

### What we have done
R Portion:
- Aim 1: (still have to do statistical tests (linear regression model))
  - Shannon Diversity boxplots
    - no clear trend of changes in Shannon diversity with age
    - and perhaps no major differences between successful and unsuccessful groups
  <img src="../R_project/aim1/shannon_diversity.png" height="300" width="400">

  - Faith's PD boxplots
 
  <img src="../R_project/aim1/faithpd_boxplot.png" height="300" width="400">

- Weighted Unifrac PCoA plot
  
- Aim 2:ISA based on reproductive outcome and age-group
  
    - ISA based on outcome: 
      - expected lactobacillus dominance in successful outcome, but found streptococcus (from the order lactobacillales)
      - maybe other taxa within the lactobacillales order may also play role in reproductive success?
      - ISA at p=0.05, 3 ASVs are associated with successful outcomes.
        <img src="../images/ISA_outcome_p0.05.png" height="80" width="1100">
      - ISA at p=0.1, 7 ASVs are associated with successful outcomes.
        <img src="../images/ISA_outcome_p0.1.png" height="110" width="1100">

    - ISA based on age group: 
      - indicated taxa are not the expected non-LD microbes
      - lactobacillus was not detected as an indicator taxon in the younger age groups
      - ISA at p=0.05: 6 ASVS are assocaited with 46-50 age group.
        <img src="../images/ISA_age_group_p0.05.png" height="110" width="1100">
      - ISA at p=0.1: 7 ASVS are assocaited with 46-50 age group, 1 with 26-30, and 1 with 31-35 and 46-50
        <img src="../images/ISA_age_group_p0.1.png" height="140" width="1100">

### Questions to ask/Issues
- When creating phyloseq object, lots of taxonomy is "Unassigned" so when we put it in proper format to be a phyloseq object lots of samples come up as "Unassigned" or "NA":
  
    <img src="../images/phyloseq_tax_q.png" height="300" width="400">
- When creating phyloseq object, is filtering it and rarefaction needed? Because losing lots of samples:
<img src="../images/filter-rarefy_Q.png" height="300" width="400"> 
- Phyloseq object rarefaction parameter (thinking of choosing 2500 as rarefaction parameter)
  <img src="../images/rare_curve.png" height="400" width="600">


## Meeting Notes

Manuscript feedback:
- While most of the feedback can be addressed and integrated, Ritu confirmed we are not expected to address every single suggestion onto our revision
  - Be more explicit in addressing certain kinds of feedback (the lactobaccilus part: "Yet, studies have not yet fully explored whether other microbial species, beyond Lactobacillus, may also contribute to reproductive success")
- Ritu would like us to explain in-depth how we are investigating the nuances in age group differences (e.g. how this is novel and more profound compared to established/published data)
  - Describe how we established our age groups (she encouraged us to attach the HISTOGRAM!)
- Ritu brought up the loss in sample size (e.g. after rarefaction) -> please adjust to minimize loss(?)

Additional Topics Discussed:
- Ritu recommended getting rid of the microbial diversity part(?) -> [if anyone knows exactly what she meant, please address this?]
- Re: regarding references from the matrix, please refer back to the provided reference page.
- Discussed TA Hans' feedback (about how infertility impacts the endometrial microbiome)
- Statistical testing: none technically completed -> asked about sample size loss

Clarified ISA
- Michelle discussed analysis results and rationale with Ritu -> Ritu suggested following through, and if Lactobaccilus doesn't show up, that's still a result
  - potential to analyze with DESeq to yield counts?
  - could expand on this during next meeting

## Week tasks

TA (pending - may take 2-3 days from today):
- Ritu will email us the statistical/correlation test she referred to during last meeting
- Ritu will also look into methods to minimize sample size loss

Re: assigned tasks
- Annalise: Aim 1 and double check Aim 3
- Michelle: Aim 2 and double check Aim 1
- Carleton: Aim 2
- Wren: Aim 3

## Future Reference
- Dr. Evelyn Sun will be around for the last 2 meetings to help clarify/troubleshoot
  - Be prepared for these last 2 meetings 
