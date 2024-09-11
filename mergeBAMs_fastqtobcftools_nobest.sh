#you should start in your fastq file and have only the fastq you want to analyse
#make the directories that will sort through as the fastq becomes a summary file

#!/bin/bash
mkdir -p {process_fastq,process_bam,summary_files,comparison_files}

#processing the fastq into a sorted.bam and moving to the other file where they are stored

#minimap2 step
for i in `ls *.fastq.gz | sed 's/.fastq.gz//g'`; 
    do
        echo -e "Processing sample \e[1;31m${i}\e[0m now";
        minimap2 -R "@RG\tID:${i}\tSM:${i}" -ax map-ont /mnt/storageG1/dan/refseq/PlasmoDB-49_Pfalciparum3D7_Genome.fasta ${i}.fastq.gz | samtools view -S -b | samtools sort - -o process_fastq/${i}.sorted.bam;
    done

cd process_fastq/

echo -e "Fastqs have been processed. Moving onto sorted.bam with files:\n$(ls *)\n\e[1;32m ################################################## \e[0m"

#mkdir of each sample and mv relevant files into the new directory
#for i in `ls *.sorted.bam | sed 's/_.*//g'`;
    #do
       #mkdir -p "$i";
       #mv "$i"* "$i"/;
    #done

#jumping into directory of each sample and merging the bams
#for dir in *;
    #do
        #if [ -d "$dir" ]; then 
            #cd $dir || { echo "Error: Unable to change to directory $dir" >&2; exit 1; }
for i in `ls *.sorted.bam | sed 's/.sorted.bam//g'`; do
            echo -e "Working in: \e[1;31m ${i}\e[0m now";
            ls *.sorted.bam > bam_list.txt;
            samtools merge -f -b bam_list.txt ../process_bam/"${i}_merged.bam";
        #fi
        #cd ../;
    done

cd ../process_bam/

echo -e "In the process_bam file directory now; creating summary files now. \n\e[1;32m ################################################## \e[0m"

#processing the sorted.bam into vcf files and using bcftools to merge them (they will need to be bgzip first to allow for indexing); to facilitate process, made each a new folder and shift them into the folder and cd into them to be processed (INSTALL NEEDED: freebayes)

for i in `ls *_merged.bam | sed 's/_merged.bam//g'`;
    do
        echo -e "Working with: \e[1;31m${i}_merged.bam\e[0m";
        freebayes -f /mnt/storage8/mtan/ref_genomes/PlasmoDB-49_Pfalciparum3D7_Genome.fasta --genotype-qualities -b ${i}_merged.bam -C 10 -q 10 --min-coverage 10 -E -1 > ${i}.vcf; 
        echo -e "\e[1;31m${i}.vcf\e[0m produced";
        bgzip ${i}.vcf;
        bcftools index -f ${i}.vcf.gz;
    done

echo -e "BAMs have been processed. Producing summary files with:\n\e[1;31m$(ls *)\e[0m\n\e[1;32m #################################################### \e[0m"

#processing vcf into tsv files

for i in `ls *.vcf.gz | sed 's/.vcf.gz//g'`;
    do
        echo -e "Working with: \e[1;31m ${i}.vcf.gz\e[0m now";
        bcftools query --print-header -f '%FILTER\t%CHROM\t%POS\t%TYPE\t%REF\t%ALT[\t%GT;\t%AD;\t%DP]\n' ${i}.vcf.gz > ../summary_files/${i}_allmut.tsv;
        echo -e "Created: \e[1;30m ${i}_removedmut.tsv\e[0m - This is the master SNP list, check this if you intend on getting other snps";
        bcftools view -V "indels,mnps" ${i}.vcf.gz > ${i}_removedmut.vcf.gz;
        bcftools query --print-header -f '%FILTER\t%CHROM\t%POS\t%TYPE\t%REF\t%ALT[\t%GT;\t%AD;\t%DP]\n' ${i}_removedmut.vcf.gz > ../summary_files/${i}_removedmut.tsv;
        echo -e "Created: \e[1;31m ${i}_removedmut.tsv\e[0m";
    done

cd ../summary_files/

echo -e "VCFs have been processed. Now in the final steps in summary_files directory. Moving onto tsv with files to produce the comparison files:\n$(ls *)\n\e[1;32m ################################ \e[0m"

#processing the merged vcf files into summary files (tsv format) and then awk out the discrepancies

for i in `ls *.tsv | sed 's/.tsv//g'`;
    do 
        head -n 1 ${i}.tsv > header.tmp;
        awk -F'\t' '$7=="."||$9=="."||$11=="."' ${i}.tsv >> header.tmp;
        mv header.tmp ../comparison_files/${i}_discrep.tsv;
        echo -e "\e[1;31m${i}_discrep.tsv\e[0m has been created";
    done
cd ../
echo -e "All done, remember to always modify this script if you want to change the conditions."

