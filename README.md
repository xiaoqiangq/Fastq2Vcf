# Fastq2Vcf
variantcalling.pl -fastq2gvcf -qplot -combinegvcf -genotypeGVCFs -vqsr -annovar -input bwa.txt -output_dir ./ -where local -sequencing_stype exome  -nthreads 14


*************************

	run -fastq2gvcf in -local:
		Example 1.1: variantcalling.pl -fastq2gvcf -input bwa.txt -output_dir ab -where local -nthreads 8
		Example 1.1.1: variantcalling.pl -fastq2gvcf -input bwa.txt -output_dir ab -where local -nthreads 8 -nthreads_perjob 4 
		Example 1.2: variantcalling.pl -fastq2gvcf -input bwa.txt -get_config calling.config -output_dir ab -where local -nthreads 8 -jobname_pre bwa
		Example 1.3: variantcalling.pl -fastq2gvcf -input bwa.txt -gvcf_intervals 22 -output_dir ab -where local -nthreads 8 
		Example 1.4: variantcalling.pl -fastq2gvcf -input bwa.txt -steps 1,2,4 -output_dir ab -where local -nthreads 8 
		Example 1.6: variantcalling.pl -fastq2gvcf -input bwa.txt -output_dir ab -where local -nthreads 8 -email *********@********* #will send email if there is an error
		Example 1.8: variantcalling.pl -fastq2gvcf -input bwa2.txt -steps 2,3,4 -output_dir ab -where local -nthreads 8 
		Example 1.9: variantcalling.pl -fastq2gvcf -input bwa2.txt -steps 4 -output_dir ab -where local -nthreads 8 
		Example 1.10: variantcalling.pl -fastq2gvcf -input bwa2.txt -steps 2,4 -output_dir ab -where local -nthreads 8 
		Example 1.11: variantcalling.pl -fastq2gvcf -qplot -combinegvcf -genotypeGVCFs -vqsr -annovar -input bwa.txt -output_dir ab -where local -nthreads 8 -sequencing_stype exome -jobname_pre 100  -steptime '13,2,11,15' -gvcf_intervals 22 -sleep 3600
		Example 1.12: variantcalling.pl -fastq2gvcf -combinegvcf -genotypeGVCFs -vqsr -input bwa.txt -output_dir ab -where local -nthreads 8
		Example 1.13: variantcalling.pl -fastq2gvcf -qplot -input bwa.txt -output_dir ab -where local -nthreads 8 -sequencing_stype exome 
		Example 1.14: variantcalling.pl -fastq2gvcf -combinegvcf -genotypeGVCFs -annovar -input bwa.txt -output_dir ab -where local -nthreads 8

	run -fastq2gvcf in -qsub or -sbatch:
		Example 2.1: variantcalling.pl -fastq2gvcf -input bwa.txt -output_dir ab -where sbatch -email *********@********* 
		Example 2.2: variantcalling.pl -fastq2gvcf -input bwa.txt -get_config calling.config -output_dir ab -where sbatch -email *********@*********  -jobname_pre bwa
		Example 2.3: variantcalling.pl -fastq2gvcf -input bwa.txt -gvcf_intervals 22 -output_dir ab -where sbatch -email *********@*********  
		Example 2.4: variantcalling.pl -fastq2gvcf -input bwa.txt -steps 1,2,4 -output_dir ab -where sbatch -email *********@*********  
		Example 2.6: variantcalling.pl -fastq2gvcf -input bwa.txt -output_dir ab -where sbatch -email *********@*********  #will send email if there is an error
		Example 2.8: variantcalling.pl -fastq2gvcf -input bwa2.txt -steps 2,3,4 -output_dir ab -where sbatch -email *********@*********  
		Example 2.9:  variantcalling.pl -fastq2gvcf -input bwa2.txt -steps 4 -output_dir ab -where sbatch -email *********@*********  
		Example 2.10:  variantcalling.pl -fastq2gvcf -input bwa2.txt -steps 2,4 -output_dir ab -where sbatch -email *********@*********  
		Example 2.11:  variantcalling.pl -fastq2gvcf -qplot -combinegvcf -genotypeGVCFs -vqsr -annovar -input bwa.txt -output_dir ab -where sbatch -email *********@*********  -sequencing_stype exome -jobname_pre 100 -steptime '13,2,11,15' -sleep 30
		Example 2.12:  variantcalling.pl -fastq2gvcf  -combinegvcf -genotypeGVCFs -vqsr -input bwa.txt -output_dir ab -where sbatch -email *********@********* 
		Example 2.13:  variantcalling.pl -fastq2gvcf -qplot -input bwa.txt -output_dir ab -where sbatch -email *********@*********  -sequencing_stype exome 
		Example 2.14:  variantcalling.pl -fastq2gvcf  -combinegvcf -genotypeGVCFs -annovar -input bwa.txt -output_dir ab -where sbatch -email *********@********* 

	run -genotypeGVCFs in -local:
		Example 3.1:  variantcalling.pl -genotypeGVCFs -input combine.g.vcf.gz -output_dir ab -where local -nthreads 20 
		Example 3.2:  variantcalling.pl -genotypeGVCFs -input combine.g.vcf.gz -output_dir ab -gvcf_intervals 22 -where local -nthreads 20 
		Example 3.3:  variantcalling.pl -genotypeGVCFs -vqsr -input combine.g.vcf.gz -output_dir ab -where local -nthreads 20 
		Example 3.4:  variantcalling.pl -genotypeGVCFs -vqsr -annovar -input combine.g.vcf.gz -output_dir ab -where local -nthreads 20  

	run -vqsr in -local:
		Example 4.1:  variantcalling.pl -vqsr -input **.vcf.gz -output_dir ab -where local -nthreads 5 
		Example 4.2:  variantcalling.pl -vqsr -annovar -input **.vcf.gz -output_dir ab -where local -nthreads 5  

	run -qplot in -local:
		Example 5.1:  variantcalling.pl -qplot -sequencing_stype exome -input bwa2.txt -output_dir ab -where local -nthreads 5 
		Example 5.2:  variantcalling.pl -qplot -sequencing_stype whole -input bwa2.txt -output_dir ab -where local -nthreads 5 
		Example 5.2:  variantcalling.pl -qplot -sequencing_stype cpg -input bwa2.txt -output_dir ab -where local -nthreads 5 
		Example 5.2:  variantcalling.pl -qplot -sequencing_stype *.bed -input bwa2.txt -output_dir ab -where local -nthreads 5 ## -dup_keep '--dup_keep' 

	run -annovar in -local:
		Example 6.1:  variantcalling.pl -annovar -input **.vcf.gz -output_dir ab -where local -nthreads 5 

	run -beagle4 in -local:

	run -combinegvcf in -local:
		Example 8.1:  variantcalling.pl -combinegvcf -input bwa3.txt -output_dir ab -where local -nthreads 1  
		Example 8.2:  variantcalling.pl -combinegvcf -input bwa3.txt -output_dir ab -where local -nthreads 1 -gvcf_intervals 1:100-200 
		X-man: variantcalling.pl -fastq2gvcf -combinegvcf -genotypeGVCFs -annovar -input bwa3.txt -output_dir ./ -where local -sequencing_stype exome -nthreads 20 -steps 4 -gvcf_intervals X -ploidy 1 

	the whole process:
		local:  variantcalling.pl -fastq2gvcf -qplot -combinegvcf -genotypeGVCFs -vqsr -annovar -input bwa.txt -output_dir ./ -where local -sequencing_stype exome  -nthreads 14 
		local: variantcalling.pl -fastq2gvcf -qplot -combinegvcf -genotypeGVCFs -vqsr -annovar -input bwa.txt -output_dir ./ -where local -sequencing_stype exome  -nthreads 14 -steps 1,2,4 
		local: variantcalling.pl -fastq2gvcf -qplot -combinegvcf -genotypeGVCFs -vqsr -annovar -input bwa.txt -output_dir ./ -where local -sequencing_stype exome  -nthreads 14 -buildver hg19 # default -buildver hg38
		cluster: variantcalling.pl -fastq2gvcf -qplot  -combinegvcf -genotypeGVCFs -vqsr -annovar -input bwa.txt -output_dir ./ -where sbatch -sequencing_stype exome -email *********@********* ## -nthreads_perjob 2 

*************************
*******************************************************
