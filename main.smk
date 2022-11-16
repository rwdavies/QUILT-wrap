from os.path import exists
import json


CHRLIST=[i for i in range(1,23)]
##CHRLIST=[i for i in range(21,23)]
##CHRLIST.append("X") # no recombination rate
RECOMB_POP="KHV" ## Vietnam
NGEN=100 ## probably about accurate



rule download_recomb:
    output:
        expand(RECOMB_POP + "/" + RECOMB_POP + "-{chr}-final.txt.gz", chr=CHRLIST)
    params:
        N='download_recomb',
        threads=1
    shell:
        """
            wget ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/working/20130507_omni_recombination_rates/{RECOMB_POP}_omni_recombination_20130507.tar
            tar -xvf {RECOMB_POP}_omni_recombination_20130507.tar
        """

rule convert_recomb:
    input:
        f"{RECOMB_POP}/{RECOMB_POP}-{{chr}}-final.txt.gz"
    output:
        f"{RECOMB_POP}/{RECOMB_POP}-chr{{chr}}-final.b38.txt.gz"
    params:
        N='convert_recomb',
        threads=1
    wildcard_constraints:
        chr='\d{1,2}'
    shell:
        """
            R -f ${{QUILT_HOME}}scripts/make_b38_recomb_map.R --args "./" {RECOMB_POP} {wildcards.chr}
        """


rule recomb:
    input:
        [
	expand(RECOMB_POP + "/" + RECOMB_POP + "-chr{chr}-final.b38.txt.gz", chr=CHRLIST)
	]





rule download_ref:
    output:
        vcf = f"refs/oneKG.chr{{chr}}.vcf.gz",
        tbi = f"refs/oneKG.chr{{chr}}.vcf.gz.tbi"
    params:
        N='download_ref',
        threads=1
    wildcard_constraints:
        chr='\d{1,2}'
    shell:
        """
            mkdir -p refs
            cd refs
            oneKG_vcf_name=CCDG_14151_B01_GRM_WGS_2020-08-05_chr{wildcards.chr}.filtered.shapeit2-duohmm-phased.vcf.gz
            wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/${{oneKG_vcf_name}}
            wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/${{oneKG_vcf_name}}.tbi
            mv ${{oneKG_vcf_name}} oneKG.chr{wildcards.chr}.vcf.gz
            mv ${{oneKG_vcf_name}}.tbi oneKG.chr{wildcards.chr}.vcf.gz.tbi
        """


rule convert_ref:
    input:
        vcf = f"refs/oneKG.chr{{chr}}.vcf.gz",
        tbi = f"refs/oneKG.chr{{chr}}.vcf.gz.tbi"
    output:
        hap = f"refs/oneKG.chr{{chr}}.hap.gz",
        legend = f"refs/oneKG.chr{{chr}}.legend.gz"
    params:
        N='convert_ref',
        threads=1
    wildcard_constraints:
        chr='\d{1,2}'
    shell:
        """
        bcftools view --output-file {input.vcf}.temp.gz --output-type z --min-alleles 2 --max-alleles 2 --types snps {input.vcf}
        tabix {input.vcf}.temp.gz
        bcftools convert --haplegendsample refs/oneKG.chr{wildcards.chr} {input.vcf}.temp.gz
        rm {input.vcf}.temp.gz
        """


rule download:
    input:
        [
        expand(RECOMB_POP + "/" + RECOMB_POP + "-{chr}-final.txt.gz", chr=CHRLIST),
        expand("refs/oneKG.chr{chr}.vcf.gz", chr = CHRLIST)
	]

rule download_only_recomb:
    input:
        [
        expand(RECOMB_POP + "/" + RECOMB_POP + "-{chr}-final.txt.gz", chr=CHRLIST)
	]


rule convert:
    input:
        [
	expand(RECOMB_POP + "/" + RECOMB_POP + "-chr{chr}-final.b38.txt.gz", chr=CHRLIST),
	expand("refs/oneKG.chr{chr}.hap.gz", chr = CHRLIST)
	]

rule convert_refs:
    input:
        [
	expand("refs/oneKG.chr{chr}.hap.gz", chr = CHRLIST)
	]

        



## ideally this should be checkpoint but I don't have the bandwidth to figure that out right now
## also ideally this would intersect with CHRLIST above but right now it's independent
REGIONS={}
for chr in CHRLIST:
    start=[10000001, 15000001]
    end=[  15000000, 20000000]
    REGIONS[str(chr)]={"start":start, "end":end}


## if file exists, load up and replace
file="regions.json"
if exists(file):
    print("Replacing regions to impute with derived file")
    with open(file) as json_file:
        REGIONS = json.load(json_file) ## python is dumb


## basically, how to determine start and end of each chromosome? want to use that to define, I think?
## do 5 Mbp chunks with 500 Mbp buffers
## cheat - pre-make chunks
## old version	
rule prep_ref:
    input:
        hap = f"refs/oneKG.chr{{chr}}.hap.gz",
        legend = f"refs/oneKG.chr{{chr}}.legend.gz",
	recomb = f"{RECOMB_POP}/{RECOMB_POP}-chr{{chr}}-final.b38.txt.gz"
    output:
        RData = f"refs/RData/ref_package.chr{{chr}}.{{regionStart}}.{{regionEnd}}.RData"
    params:
        N='prep_ref',
        threads=1
    wildcard_constraints:
        chr='\w{1,2}',
        regionStart='\d{1,9}',
        regionEnd='\d{1,9}'
    shell:
        """
            mkdir -p refs/RData/other/
            ## buffer set to 0, here we impute the buffers, get rid of them later
            ## using concat
            R -e 'library("data.table"); library("QUILT"); QUILT_prepare_reference( \
            outputdir="refs/RData/other/", \
            chr="chr{wildcards.chr}", \
            nGen={NGEN}, \
            reference_haplotype_file="{input.hap}" ,\
            reference_legend_file="{input.legend}", \
            genetic_map_file="{input.recomb}", \
            regionStart={wildcards.regionStart}, \
            regionEnd={wildcards.regionEnd}, \
            buffer=0, \
            output_file="{output.RData}")'
        """

regions_to_prep=[]
vcfs_to_impute=[]
for chr in CHRLIST:
    start=REGIONS[str(chr)]["start"]
    end=REGIONS[str(chr)]["end"]
    for i in range(0, start.__len__()):
        regionStart=start[i]
        regionEnd=end[i]
        file="refs/RData/ref_package.chr" + str(chr) + "." + str(regionStart) + "." + str(regionEnd) + ".RData"
        regions_to_prep.append(file)
        file="vcfs/regions/quilt.chr" + str(chr) + "." + str(regionStart) + "." + str(regionEnd) + ".vcf.gz"
        vcfs_to_impute.append(file)
	



rule prep:
    input:
        [
	regions_to_prep
	]

##        expand("refs/oneKG.chr{chr}.hap.gz", chr = CHRLIST),
##         regions_to_prep


rule quilt:
    input:
        bamlist = "bamlist.txt",
        RData = f"refs/RData/ref_package.chr{{chr}}.{{regionStart}}.{{regionEnd}}.RData"
    output:
        vcf = f"vcfs/regions/quilt.chr{{chr}}.{{regionStart}}.{{regionEnd}}.vcf.gz"
    params:
        N='impute',
        threads=1
    wildcard_constraints:
        chr='\w{1,2}',
        regionStart='\d{1,9}',
        regionEnd='\d{1,9}'
    shell:
        """
            ## set a seed here, randomly, so can try to reproduce if it fails
            SEED=`echo $RANDOM`
            mkdir -p vcfs/regions/
            R -e 'library("data.table"); library("QUILT"); QUILT( \
            outputdir="refs/RData/other/", \
            chr="chr{wildcards.chr}", \
            regionStart={wildcards.regionStart}, \
            regionEnd={wildcards.regionEnd}, \
            buffer=0, \
            bamlist="{input.bamlist}", \
            prepared_reference_filename="{input.RData}", \
            output_filename="{output.vcf}", \
            seed='${{SEED}}')'
        """

rule quilt_info:
    input:
        vcf = f"vcfs/regions/quilt.chr{{chr}}.{{regionStart}}.{{regionEnd}}.vcf.gz"
    output:
        vcf = f"vcfs/regions/quilt.chr{{chr}}.{{regionStart}}.{{regionEnd}}.vcf.gz.output.RData"
    params:
        N='info',
        threads=1
    wildcard_constraints:
        chr='\w{1,2}',
        regionStart='\d{1,9}',
        regionEnd='\d{1,9}'
    shell:
        """
            R -f ${{QUILT_WRAP_HOME}}info.R --args {output.vcf}
        """


rule impute_per_chr:
    input:
        [
	vcfs_to_impute
	]


vcfs_to_concat={}
final_vcfs=[]
for chr in CHRLIST:
    start=REGIONS[str(chr)]["start"]
    end=REGIONS[str(chr)]["end"]
    per_chr_vcfs=[]
    for i in range(0, start.__len__()):
        regionStart=start[i]
        regionEnd=end[i]
        file="vcfs/regions/quilt.chr" + str(chr) + "." + str(regionStart) + "." + str(regionEnd) + ".vcf.gz"
        per_chr_vcfs.append(file)
    vcfs_to_concat[str(chr)]=per_chr_vcfs
    final_vcfs.append("vcfs/quilt.chr" + str(chr) + ".vcf.gz")



def get_input_vcfs_as_list(wildcards):
    return(vcfs_to_concat[str(wildcards.chr)])

def get_input_vcfs_as_string(wildcards):
    return(" ".join(map(str, vcfs_to_concat[str(wildcards.chr)])))


rule concat:
    input:
        vcfs = get_input_vcfs_as_list
    output:
        vcf = f"vcfs/quilt.chr{{chr}}.vcf.gz"
    params:
        N='impute',
        threads=1,
        input_string=get_input_vcfs_as_string
    wildcard_constraints:
        chr='\w{1,2}',
        regionStart='\d{1,9}',
        regionEnd='\d{1,9}'
    shell:
        """
            bcftools concat \
            --ligate \
            --output-type z \
            --output {output.vcf}.temp1.vcf.gz \
            {params.input_string}

            ## remove HD entry as it is now unstable after ligations
            gunzip -c {output.vcf}.temp1.vcf.gz | grep '#' > {output.vcf}.temp2.vcf
            bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\t%FILTER\t%INFO\tGT:GP:DS:PS\t[%GT:%GP:%DS:%PS\t]\n' {output.vcf}.temp1.vcf.gz  >> {output.vcf}.temp2.vcf
            bgzip {output.vcf}.temp2.vcf
            tabix {output.vcf}.temp2.vcf.gz

            mv {output.vcf}.temp2.vcf.gz {output.vcf}
            mv {output.vcf}.temp2.vcf.gz.tbi {output.vcf}.tbi
            rm {output.vcf}.temp1.vcf.gz
        """




rule impute:
    input:
        [
        final_vcfs
	]
