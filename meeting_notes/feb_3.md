# February 3, 2025

## Agenda
- Chosen topic: IVF and the endometrial microbiome
  
[Paper we got dataset from](https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-021-01184-w)

Metadata variables within dataset: age, biomaterial provider, biosample model, disease, isolate, library name, organism, sex, tissue

Screenshot of dataset:

<img src="../images/dataset.png" height="650" width="650">

### Potential Topic Questions: 
Age related: data set includes age 21-49
- How does the distribution of endometrial microbiota sample types (EF vs. EB) vary across different age groups?
- Does age correlate with reproductive success (live birth vs. no pregnancy)?

Biomaterial provider: metadata includes 19 different biomaterial providers (1-19)
- Do microbiota profiles differ between samples from different biomaterial providers?
- Could sample processing techniques at different centers influence microbiota composition and affect the study’s conclusions?
- paper mentioned 13 different centres on three continents. howver did not specify which centres corresponds to provider on metadata.

Other potential questions: the paper mentioned the ethnic distribution: caucasian 57.3%, east asian 14.0%, hispanic 11.4%, others 17.3%. If possible to obtain data for metadata category, could also look at potential differences in microbiota composition in different ethnicities.

## Meeting Notes

- Evelyn putting paper and dataset onto server and test it for us, will send this back to us

Chosen Research Question:
- How does distribution of endometrial microbiota sample types (EF vs EB) vary across different age groups, and does age correlate with reproductive success?
  
To Do:
- Use given metadata and load into R to find sample sizes of age
  - Make histogram to see how many samples we have of each age
  - If sample sizes are small, can bin ages together (ex: 20-25, 26-30, etc) 
- Within each age, measure the frequency of each pregnancy outcome
- Need a back-up plan in case the metadata doesn't work within the server (could maybe use the paper 9 data and just work with a smaller sample size)
- Process metadata using Qiime before proposal (once on the server) using the team server
- Download project 1 qza files onto GitHub as reference and easy access
- Keep copy of metadata file, and other files generated as backups

Next week
- Discuss writing an outline for the proposal

Proposal assignment
- Posted on CANVAS
- Expected: data processing already done before proposal is due


