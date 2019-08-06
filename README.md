# ParseCNV2: CNV GWAS Tool
Parse Copy Number Variation from Array and Sequencing  
perl ParseCNV2.pl  
perl ParseCNV2.pl -i v4.3.vcf -c Cases_vcf.list -b hg19 -m 1  
perl ParseCNV2.pl -i Cases.rawcnv -c Cases.list -b hg19  
perl ParseCNV2.pl -i Cases.rawcnv -q Cases.qt -m 1 -d 2 -p 1 -b hg19 -batch ID_Category.txt  
perl ParseCNV2.pl -i Cases.rawcnv -qc -bfile plinkBed -log penncnv.log -c Cases.list2 -b hg19  
