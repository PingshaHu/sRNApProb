#################################
### posterior prob estimation ###	
### random_sampling_v1.02.r   ###
################################
random_sampling_v1.02.r estimate the posterior probability for smRNA loci based on its coverage and length distribution.  
Usage example:
R --vanilla --slave --args smRNA_mergecounts.txt BW0,BW0,BW0,BWM13,BWM13,BWM13,BWM21,BWM21,BWM21,BWM4,BWM4,BWM4,BWM9,BWM9,BWM9,BWZt13,BWZt13,BWZt13,BWZt21,BWZt21,BWZt21,BWZt4,BWZt4,BWZt4,BWZt9,BWZt9,BWZt9,CDB,CDB,CDB < random_sampling_v1.02.r

BW0,BWM13 etc are different treatments, which is corresponding to the smRNA_mergecounts.txt file columns.