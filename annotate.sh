#!bin/bash

mkdir -p annotate

conda deactivate
    for i in `ls *_core_removedmut.vcf.gz | sed 's/_core_removedmut.vcf.gz//g'`;
        do mv ${i}_core_removedmut.vcf.gz ${i}_core_removedmut.vcf;
        bcftools annotate --rename-chrs /mnt/storage8/mtan/ref_genomes/pf3d7_resistance/chrs.txt ${i}_core_removedmut.vcf -o ${i}_renamed.vcf
        done
    
conda activate snpeff
    for i in `ls *_renamed.vcf | sed 's/_renamed.vcf//g'`;
        do snpEff ann -v Plasmodium_falciparum ${i}_renamed.vcf > annotate/${i}_ann.vcf;
        echo "${i}_ann.vcf created";
        done
conda deactivate
cd annotate/
    for i in `ls *_ann.vcf | sed 's/_ann.vcf//g'`;
    do bcftools query -f '%CHROM\t%POS\t%ALT\t%ANN\n' ${i}_ann.vcf > ${i}_ann_vars.tab;
    echo "${i}_ann_var.tab created";
    done

cd ../