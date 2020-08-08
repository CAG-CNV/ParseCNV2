# ParseCNV2: CNV GWAS Tool
Parse Copy Number Variation from Array and Sequencing  
#Make sure compiles and Print Usage Message  
perl ParseCNV2.pl  
#Specify a sequencing CNV call VCF, affected samples, and genome build  
perl ParseCNV2.pl -i v4.3.vcf -c Cases_vcf.list -b hg19 -m 1  
#Specify an array CNV call text file, affected samples, and genome build  
perl ParseCNV2.pl -i Cases.rawcnv -c Cases.list -b hg19  
#Specify an array CNV call text file, quantitative trait for samples, genome build, and batch to track  
perl ParseCNV2.pl -i Cases.rawcnv -q Cases.qt -m 1 -d 2 -p 1 -b hg19 -batch ID_Category.txt  
#Specify an array CNV call text file, perform quality control, plink format SNP genotype bed file, PennCNV detect_cnv log, affected samples, and genome build  
perl ParseCNV2.pl -i Cases.rawcnv -qc -bfile plinkBed -log penncnv.log -c Cases.list2 -b hg19  

![SNPArrayGenotypingData](/images/SNPArrayGenotypingData.png)  
  
![SequencingData](/images/SequencingData.png)  

perl ParseCNV2.pl -i CNV_Calls[.vcf/.rawcnv] [-c Cases.txt/-q Samples_QuantitativeTrait.txt] [arguments]  

Required arguments:  

-i      input[vcf/rawcnv/txt]  
-c      cases  
-b      build (if vcf can be inferred)  

Optional arguments:  

-o      output  
-p      merge p variation  
-d      merge distance  
-t      transmission disequilibrium test  
-m      max p inclusion  
-bfile  plink bed file of SNP genotypes  
-qc     quality control upfront  
-log    penncnv log  
-q      quantitative trait  
-b      batch  
-h      help  
