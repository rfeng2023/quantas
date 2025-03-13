##################################
#For each sam file with name '$f'#
##################################
#Annotation files
intron_file=hg19.intron.hmr.inclusive.plus.Lister.bed #hg38.intron.hmr.inclusive.plus.Lister.bed
gene_file=refGene_knownGene_hg19.bed #refGene_knownGene_hg38.bed

#1.Sam2bed
perl sam2bed.pl -v -u $f.sam $f.bed

#2.Count
perl summarize_splice_site_usage.pl -big --anchor 5 -v -gene $gene_file $intron_file $f.bed $f.ss.count.txt


#3.Splicing site usage
#file information in conf file
#An example of input.onegroup.conf file: for one sam/bam file as input, only the first line is needed:
# $f.ss.count.txt	cortex 
# hs_fc_53yr_run2.ss.count.txt	cortex
# hs_fc_53yr_run3.ss.count.txt	cortex
# hs_fc_53yr_run4.ss.count.txt	cortex
# hs_mfg_25yr_run1.ss.count.txt	cortex
# hs_mfg_25yr_run2.ss.count.txt	cortex
perl gen_splicing_matrix.pl -type ss --min-cov 20 --max-std 0.1 --na-string NaN --print-info -v input.onegroup.conf $f.ss.ssu.txt