# ParseCNV2: CNV GWAS Tool
Parse Copy Number Variation from Array and Sequencing  

Cite: Glessner JT, Li J, Liu Y, Khan M, Chang X, Sleiman PMA, Hakonarson H. ParseCNV2: efficient sequencing tool for copy number variation genome-wide association studies. Eur J Hum Genet. 2022 Nov 1. doi: 10.1038/s41431-022-01222-7. Epub ahead of print. PMID: 36316489.  
Original: Glessner JT, Li J, Hakonarson H. ParseCNV integrative copy number variation association software with quality tracking. Nucleic Acids Res. 2013 Mar 1;41(5):e64. doi: 10.1093/nar/gks1346. Epub 2013 Jan 4. PMID: 23293001; PMCID: PMC3597648.  

#Test Commands Updated to include 2 options needed for test data and latest options -stat fisher -no_freq  
#Make sure compiles and Print Usage Message  
perl ParseCNV2.pl  
#Specify a sequencing CNV call VCF, affected samples, and genome build  
perl ParseCNV2.pl -i v4.3.vcf -c Cases_vcf.list -b hg19 -m 1 -stat fisher -no_freq  
#Specify an array CNV call text file, affected samples, and genome build  
perl ParseCNV2.pl -i Cases.rawcnv -c Cases.list -b hg19 -stat fisher -no_freq  
#Specify an array CNV call text file, quantitative trait for samples, genome build, and batch to track  
perl ParseCNV2.pl -i Cases.rawcnv -q Cases.qt -m 1 -d 2 -p 1 -b hg19 -batch ID_Category.txt -stat fisher -no_freq  
#Specify an array CNV call text file, perform quality control, plink format SNP genotype bed file, PennCNV detect_cnv log, affected samples, and genome build  
perl ParseCNV2.pl -i Cases.rawcnv -qc -bfile plinkBed -log penncnv.log -c Cases.list2 -b hg19 -stat fisher -no_freq  
#Specify an array CNV call text file, affected samples, and genome build, covariates file, and statistic as logistic  
perl ParseCNV2.pl -i Cases.rawcnv -c Cases.list -b hg19 -covar SamplesAgeSexRace.txt -stat logistic -no_freq  

![SNPArrayGenotypingData](/images/SNPArrayGenotypingData.png)  
  
![SequencingData](/images/SequencingData.png)  

Update: ParseCNV2_6-28-22.pl is the last version of the original model of convert all inputs to rawcnv2 format  
ParseCNV2.pl now has several upgrades including the new model of convert all inputs to vcf format and integrating Plink2 and RvTests  

perl ParseCNV2.pl -i CNV_Calls[.vcf/.rawcnv] [-c Cases.txt/-q Samples_QuantitativeTrait.txt] [arguments]  

Required arguments:  

-i      input[vcf/rawcnv/txt]  
-c      cases  
-b      build (if vcf can be inferred)  

Optional arguments:  

-o      output  
-covar  covariates file  
-stat   statistic (fisher, logistic, or linear) or RvTests:  
        Single variant: score, wald, exact, dominantExact, famLRT, famScore, famGrammarGamma, firth  
        Burden: cmc, zeggini, mb, fp, exactCMC, cmcWald, rarecover, cmat, famcmc, famzeggini  
        Variable threshold: price, analytic, famAnalytic  
        Kernel: skat, skato, kbac, famSkat  
        Meta-Analysis: score, dominant, recessive, cov, bolt, boltCov  
-no_freq        No CNV frequency filtering, default is rare (minor allele frequency <0.01) recurrent (minor allele count > 2)  
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
