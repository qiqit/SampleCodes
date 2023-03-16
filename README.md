# Sample Codes

This repo contains a few sample codes for secondary metabolism analysis on plant-asssociated bacteria.  

**BGC**: biosynthetic gene clusters (predicted by antiSMASH software)

### antiSMASH_to_gff.py
Converts gbk format output of antiSMASH analysis to gff format. Used for scoary/roary analysis. Only CDS genes are considered.  
```
    *** convert antiSMASH generated gbk file to gff format ***
    -F folder that contains the antiSMASH result from a bacetria genome
    -O output foldername
```
### heatmap_data.py
Conver Biom file to three dataframes containing whole BGC coverage score, core BGC coverage score and a normalized RPKM value for heatmap plotting. Normalization is based on the **metagenomeSeq R package**. 
```
    *** process the result.ALL file from metagenome.mapping.py ***
    *** filter BGCs with low coverage score and re-arrange the data ***
    -B1 biom file for whole gene cluster mapping (required)
    -O output folder name (required): output folder name
    -B2 biom file for core gene cluster mapping (optional): include only when you want to include a filtration based on core coverage
    -c1 cov_cutoff (optional): default= 0.5
    -c2 corecov_cutoff (optional): default= 0.5
```
### hypergeometric.R
- Perform hypergeometric test on two bacterial groups (plant-associated and non-plant-associated in the current script) based on the abundance each BGC family.  
- Need output file "BGCF_matrix_0.2_phylum.csv" from ` Generate_BGCF_matrix.R `. 
- The null hypothesis is that there is no significant difference in the abundance of BGCs between the plant-associated and non-plant-associated groups, while the alternative hypothesis is that there is a significant difference.  
- Two versions of tests are performed: (1) raw counts of BGCs belongs to each BGC type (copy_number) (2) presence or absence of each BGC type in a bacteria genome (bin).

### mapping_heatmap.R
- Draw a complex heatmap of mapping results.  
- Need output file "BGC_mapping_result_all.csv" from ` heatmap_data.py `. 
- The heatmap is based on abundance of every BGC in every metagenome samples.  
- Whole coverage, core coverage, marks to known BGC are included in the annotation. 
