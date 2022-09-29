## code for code run to impute samples
## rwdavies

## cd to this directory and get environment variables
cd ~/proj/QUILT-wrap
source activate

## source previously installed conda
. ~/.conda_activate
## otherwise, see README, or try code like the below
## wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
## bash Miniconda3-latest-Linux-x86_64.sh

## install snakemake and QUILT from CONDA
conda create --name quilt-wrap
conda install -n quilt-wrap snakemake
conda install -n quilt-wrap r-quilt
conda activate quilt-wrap

## install missing R packgae
## not needed if using R -e '()' 
## R -e 'install.packages("optparse", repos="http://cran.us.r-project.org")'

## install liftOver
cd ${ANALYSIS_DIR}
wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/liftOver
chmod +x liftOver
wget http://hgdownload.cse.ucsc.edu/goldenpath/hg19/liftOver/hg19ToHg38.over.chain.gz

## Manually set RECOMB_POP="KHV" in main.smk
## TODO - set this more algorithmically, or just choose an African / African American population

## Download ref information
./run.sh download local 4 ## done locally as only head nodes have internet

## Convert recombination rate files and reference files
./run.sh convert cluster 100

## Determine per-chromosome chunks
## this requires the conversion done in convert above
## this could be done using Snakemake programatically, but we do using R for simplicity
R -f determine_chunks.R ## done locally on head node

## Prepare per-region imputation files
./run.sh prep cluster 1000













##
## Azim specific work - though note - recombination map above is specific 
## 


## make bamlist
## manually remove 4 high coverage bams
BAMLIST=${QUILT_WRAP_HOME}bamlist.txt
ls /well/ansari/shared/lcwgs/data/raw_data_1x/*bam | \
    grep -v 'WTCHG_904268_72135285.bam\|WTCHG_904268_72125284.bam\|WTCHG_904268_72115283.bam\|WTCHG_904268_72105282.bam' > ${BAMLIST}

## Impute samples
./run.sh impute cluster 1000






exit




## scratch


## Why am I not just doing 1000 Genomes NYGC VCFs? Then converting?


## Manually download 1000 Genomes data in correct format
cd ${ANALYSIS_DIR}
wget https://mathgen.stats.ox.ac.uk/impute/1000GP_Phase3.tgz
tar -xzvf 1000GP_Phase3.tgz
wget https://mathgen.stats.ox.ac.uk/impute/1000GP_Phase3_chrX.tgz
tar -xzvf 1000GP_Phase3_chrX.tgz
mv *X*.gz 1000GP_Phase3/
mv genetic_map_chrX_nonPAR_combined_b37.txt genetic_map_chrX_NONPAR_combined_b37.txt
mkdir refs
mkdir recomb

for chr in `echo $(seq 1 22) X_PAR1 X_PAR2 X_NONPAR`
do
    echo transform ${chr}
    if [ ! -e refs/oneKG.chr${chr}.hap.filtered.gz ]
    then
        mv 1000GP_Phase3/1000GP_Phase3_chr${chr}.legend.gz refs/oneKG.chr${chr}.legend.gz
        mv 1000GP_Phase3/1000GP_Phase3_chr${chr}.hap.gz refs/oneKG.chr${chr}.hap.gz
        gunzip -c refs/oneKG.chr${chr}.legend.gz | awk '{if($5 == "Biallelic_SNP") {print NR-1}}' > refs/temp.chr${chr}.txt
        awk '{if(NR==FNR) {a[$1]1;next} {if((FNR==1) || ((FNR-1) in a)) {print $0}}}' refs/temp.chr${chr}.txt <(gunzip -c refs/oneKG.chr${chr}.legend.gz) | gzip -1 > refs/oneKG.chr${chr}.legend.filtered.gz
        awk '{if(NR==FNR) {a[$1]1;next} {if((FNR==1) || ((FNR) in a)) {print $0}}}' refs/temp.chr${chr}.txt <(gunzip -c refs/oneKG.chr${chr}.hap.gz) | gzip -1 > refs/oneKG.chr${chr}.hap.filtered.gz
    fi
done


## fix recomb
(cd ~/proj/QUILT/snakemake && ./run.sh download_only_recomb local 4)
for chr in $(seq 1 22)
do
    echo "position COMBINED_rate.cM.Mb. Genetic_Map.cM." > KHV/KHV-${chr}-final.fixed.txt
    gunzip -c KHV/KHV-${chr}-final.txt.gz | awk 'NR>1' | cut -f1-3 | sed -e 's/\s\+/ /g' | sed -e 's/^ *//g'  >> KHV/KHV-${chr}-final.fixed.txt
    gzip -1 -f KHV/KHV-${chr}-final.fixed.txt
done


cd ~/proj/QUILT/snakemake

## Determine per-chromosome chunks
## this requires the conversion done in convert above
## this could be done using Snakemake programatically, but we do using R for simplicity
R -f determine_chunks.R

## Prepare per-region imputation files
./run.sh prep cluster 1000

## Impute samples
./run.sh impute cluster 1000
