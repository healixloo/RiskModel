cd /Users/jlu/Desktop/Pro_Mine/Paper_TCGA_2/script
Rscript step3_moduleGene.Diff.R

cd /Users/jlu/Desktop/Pro_Mine/Paper_TCGA_2/data_in
perl step4_moduleGene_mergeExpTime.pl

cd /Users/jlu/Desktop/Pro_Mine/Paper_TCGA_2/script
Rscript step4pl_moduleGene.uniCox.R 
Rscript step5_moduleGene.multiCox.R 
Rscript step6_moduleGene.survial.R

