# Details for the QC of WGS data after joint calling

This pipeline will performed some QC steps on a complete dataset considering the following parameters:

Sample related:

* PCA
* Singleton distribution
* Missingness
* Coverage
* Heterozygosity rate

Variant related:

* HWE
* Heterozygosity rate
* Call rate

Samples and variants comparison with SNP array data (available only for samples of the FVG cohort):

* MAF mismatch with SNP array data
* NRDR

