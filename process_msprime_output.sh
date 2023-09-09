## Processing msprime sim output ##
source /share/hennlab/progs/miniconda3/etc/profile.d/conda.sh
conda activate msprime

# Set up
filename="decliningNe_Ne1-4000_Ne2-450_n120_gen6"
outdir="/share/hennlab/projects/himba/himba_msprime/${filename}/"
recomb_map="/share/hennlab/reference/recombination_maps/genetic_map_AfricanAmerian/AAmap.superChr.map"
scripts_dir="/share/hennlab/projects/himba/himba_msprime/decliningNe_Ne1-4000_Ne2-450_n120_gen60/scripts/"
MAF=0.1
cd $outdir

# Update chromosome numbers in vcf file
python ${scripts_dir}format_msprime_vcf.py -v ${filename}.vcf -g $recomb_map -o ${filename}_updated.vcf
# Update tsk_0 individual
sed -i "0,/tsk_0/s//tsk_O/" ${filename}_updated.vcf
module load vcftools
bgzip ${filename}_updated.vcf
tabix ${filename}_updated.vcf.gz

# Remove duplicate variants
module load bcftools
bcftools norm -d both ${filename}_updated.vcf.gz -o ${filename}_updated_noDups.vcf

# Rename variants in vcf file
grep -n "^#" ${filename}_updated_noDups.vcf > temp
N=`cut -f1 -d":" temp |tail -n1`
cut -f2 -d":" temp > header.txt
sed "1,${N}d" ${filename}_updated_noDups.vcf | awk -F"\t" '{OFS=FS}{$3=$1":"$2; print}' > columns.txt

cat header.txt columns.txt > ${filename}_updated_noDups.vcf
rm -f header.txt columns.txt temp

# Convert to plink and calculate frequency
module load plink2/2.00a
plink2 --vcf ${filename}_updated_noDups.vcf --recode --make-bed --out ${filename}_updated_noDups
plink2 --bfile ${filename}_updated_noDups --freq --out ${filename}_updated_noDups_freq --silent
plink2 --bfile Himba_H3Africa_unrel4th --freq --out Himba_H3Africa_unrel4th_freq --silent
awk -v OFS="\t" '$1=$1' Himba_H3Africa_unrel4th_freq.afreq > test.afreq ; mv test.afreq Himba_H3Africa_unrel4th_freq.afreq
awk -v OFS="\t" '$1=$1' ${filename}_updated_noDups_freq.afreq > test.afreq ; mv test.afreq ${filename}_updated_noDups_freq.afreq

# Thin SNPs to approximately match a template SNP array
python ${scripts_dir}SNP_array_match.py --data ${filename}_updated_noDups --target Himba_H3Africa_unrel4th --outdir $outdir --drop_target_mono True --norm_target True --add_mono_vars False --fold_target True
awk '{gsub(":","\t",$0);  print;}' vars_to_extract.list > temp; mv temp vars_to_extract.list
module load vcftools
vcftools --vcf ${filename}_updated_noDups.vcf --positions vars_to_extract.list --recode --out ${filename}_updated_noDups_SNP-array-thinned

# Add the colon back into the rsid column
#awk 'BEGIN { OFS=FS } {gsub(/ /,":", $3);  print}' ${filename}_updated_noDups_SNP-array-thinned.recode.vcf > temp_recode.vcf ; mv temp_recode.vcf ${filename}_updated_noDups_SNP-array-thinned.recode.vcf

# Calculate MAX_MAF for next step
MAX_MAF=`echo 1 - $MAF|bc`

# Break out into chromosomes
#for chr in `seq 1 22`
#do
 #   vcftools --vcf ${filename}_updated_noDups_SNP-array-thinned.recode.vcf --chr $chr --maf $MAF --max-maf $MAX_MAF --recode --out ${filename}_chr${chr}
    # Convert to haps/sample format
  #  bcftools convert --hapsample --vcf-ids ${filename}_chr${chr}.recode.vcf -o ${filename}.chr${chr}.phased
   # gunzip ${filename}.chr${chr}.phased.hap.gz; mv ${filename}.chr${chr}.phased.hap ${filename}.chr${chr}.phased.haps
    #rm -f ${filename}_chr${chr}.recode.vcf ${filename}_chr${chr}.log
#done

# Make binary format for IBDNe pipeline
#plink2 --vcf ${filename}_updated_noDups_SNP-array-thinned.recode.vcf --recode --make-bed --out ${filename}_final

# Fix rsID values in bim file - give each a unique ID
#awk -F"\t" '{OFS=FS}{$2=$1":"$4; print}' ${filename}_final.bim > temp
#mv temp ${filename}.bim

