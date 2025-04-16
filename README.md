full workflow to sequence RNA from Oxford Nanopore Flowcells using Dorado/Minimap2 and running differential analysis using salmon & DEseq2, or differential modification using m6anet & Xpore
**Currently set for RNA004 Datasets on experimental models**

**https://github.com/nanoporetech/dorado**

**for me, the following modules seemed to be needed to run dorado without errors**
    
    module load GCC
    module load PyTorch
    module load FlexiBLAS/3.3.1-GCC-12.3.0  
    module load FFmpeg/4.4.2-GCCcore-11.3.0 
    module load HTSlib
    module load protobuf
    
run the following command for the basecall step, ensure the correct model and directory that data is pulling from.

    /path/to/updated/bin/dorado \
    basecaller \
    /path/to/desired/version/dorado-0.8.2/dorado-0.8.2-linux-x64/models/rna004_130bps_sup@v5.1.0 \
    -r "/path/to/nanopore/rawdata/pod5/" \
    --min-qscore 10 \
    > /path/to/output/file.bam
    
#add --modified-bases tag for context-specific basecalling.

Alignment can be done with dorado, but minimap2 is preferred for downstream analysis purposes.
MarginAlign has built-in minimap, and provides more statistical data (MarginAlign will need to be run on Python 2.7)
**https://github.com/benedictpaten/marginAlign**

Prep the data to run statistical tests
        
    module load SAMtools/1.3-intel-2016a 

    samtools fastq path/to/basecalled.bam > path/to/basecalled.fastq
    
**added path to my marginalign download & necessary modules to prevent additional downloads**

    module load QUAST/5.0.2-foss-2018b-Python-2.7.15
 
    "/ourdisk/hpc/rnafold/gjandebeur/dont_archive/software/marginAlign/submodules/minimap2/minimap2"
    -ax splice -uf -k14 "/path/to/reference/genome.fa"     "/path/to/basecalled.fastq" \
    > "path/to/minimap/aligned.sam"

**for -ax splice, it uses spliced alignment. Change to -ax map-ont if an unspliced alignment is needed (for m6anet and Xpore)**
  **from this point, there are multiple next steps, not one direct workflow**

 Next steps will sort the bam file, filter out secondary/supplementary reads, as well as run an index and flagstat on the files.

    samtools sort -o /path/to/sorted.bam /path/to/minimap/aligned.sam
    samtools view -h -F 0x800 /path/to/sorted.bam | samtools view -b -o /path/to/sorted/filtered.bam
    samtools index /path/to/sorted/filtered.bam
    samtools flagstat /path/to/sorted/filtered.bam #this will give the basecalled #, mapped #, and mapping %


**this next step looks into the basecalled data for statistics, used in creating the substitution matrix and identity plots** 
To run marginStats you will need Python2.7.15, create a virtual environment/conda with this before running.

         marginStats \
        --localAlignment /path/to/filtered.bam/ \
        /path/to/uniqueheaders.fastq/ \
        /path/to/reference/genome.fa \
        --readIdentity --alignmentIdentity \ 
        --readCoverage \
        --mismatchesPerAlignedBase \
        --deletionsPerReadBase \
        --insertionsPerReadBase \
        --readLength --printValuePerReadAlignment 
        
With the marginstats output, run substitutionplot.R and Identity.R scripts

**to prep for Differential expression analysis**
First use salmon to quantify abundance of the expression in each sample's transcript (TPM = transcript abundance)
**https://github.com/COMBINE-lab/salmon**
**You can install salmon completely into a conda env, needing no real installation**

    salmon index -t "/path/to/reference/transcriptome.fa" -i         
    "/path/to/reference/transcriptome_index" -k 31 

    salmon quant -t \
     "/path/to/reference/transcriptome.fa" -l A -a          
     "/path/to/aligned/sorted/filtered.bam" -p 8 --seqBias --gcBias -o         
     "/path/to/output/directory"

This output can then be used to run DEseq2 in an R script. I'll attach one of mine, but it is a work in progress.

**The next section uses the aligned data to produce a consensus sequence for the entire script, (transitions from .bam to a readable transcript (bam forces ACTG, but should be ACUG)**


     samtools mpileup -uf /path/to/reference/genome.fa \
     -Q 10 -q 10 \
     /path/to/the/aligned.bam | bcftools call -c -o /path/to/variants/file.vcf

zip up the vcf file

        bgzip -c /path/to/variants/file.vcf/ > /path/to/variants.vcf.gz

to create the consensus sequence, 

    bcftools index /path/to/variants.vcf.gz/

    bcftools consensus -f /path/to/reference/genome.fa \
    /path/to/variants.vcf.gz/ \
    > /path/to/output/consensus/sequence.fq

Last step is converting from T to U since this focuses on RNA (bam doesn't write Uracil, only Thymine)

      sed 's/T/U/g' /path/to/output/consensus/sequence.fq > /path/to/output/rna/consensus/sequence.fq

**The following section is the dataprep for running differential modification using m6anet and Xpore**

**https://github.com/hasindu2008/f5c**

**f5c is built to be the updated version of nanopolish that allows RNA004 modification detection**

        source /path/to/conda/env.sh
        conda activate $ENV_PATH
        conda new (set up nanopolish environment for dependencies, should be on bioconda)



 install pod5 to convert from pod5 to fast5 (f5c only reads the old fast5 or new slow5, not dorado's pod5 format)

    pod5 convert to_fast5 "/path/to/the/pod5/files/pod5/" \
    --output "/path/to/new/directory/for/fast5/"
    
With library prep done, perform f5c index on the sample to rewrite the original ionic current from the flowcell for differential analysis

    /path/to/f5c \
    index -d "/path/to/the/file/fast5/flowdata.fast5 \
    /path/to/fastq/basecalled.fastq

**f5c also allows for m5C detection, but false positive rate seems to be higher than preferred (showing high m5C detection on IVT samples)**

    /path/to/f5c \
    call-methylation --min-mapq 20 \
    -r /path/to/basecalled.fastq \
    -b /path/to/filtered/and/sorted.bam \ 
    -g /path/to/reference/genome.fa \
    --pore rna004 > /path/to/output_results.tsv

    /path/to/f5c \
    meth-freq -i /path/to/output_results.tsv > /path/to/output_frequency.tsv

**use f5c to run eventalign on the rna004 files**
#This will generate a VERY large file, always send to /scratch/

    /path/to/f5c eventalign --rna \
    -b /path/to/aligned/and/sorted.bam \
    -r /path/to/basecalled.fastq \
    -g /path/to/reference.fa \
    --kmer-model /path/to/F5C/experimentalrna004/5mer.model \
    --signal-index --scale-events > /path/to/output/f5cevent.txt

**once complete the event file can be run with m6anet**
**https://github.com/GoekeLab/m6anet**
**best bet is to use my m6anet directory, I've had to slightly change code for there brand new model, the glori model isnt built into the software last I checked**
**also use a python3.8 conda environment**

        m6anet dataprep --eventalign /path/to/event.txt \
        --out-dir /path/to/output/directory/ --n_processes 4 

        m6anet inference --input_dir /path/to/input/from/dataprep/ \
        --out-dir /path/to/m6anet/inference/directory/ \
        --pretrained_model HEK293T_RNA004_glori --n_processes 4 --num_iterations 1000

**using the same event.txt file from f5c, run Xpore**
**https://github.com/GoekeLab/xpore**

        xpore dataprep --eventalign /path/to/event.txt \
        --out_dir /path/to/output/xpore/directory/

To run the differential expression you will need a .yml config file, here's an example:
        notes: Pairwise comparison without replicates with default parameter setting.

    data:
        KO:
            rep1: #path/to/the/modification/xporedataprep
        IVT:
            rep1: #path/to/the/IVT/xporedataprep
        
    out: /path/to/output/location

    method:
        prefiltering:
            method: t-test
            threshold: 0.1

and then run Xpore using
        
        xpore diffmod --config /path/to/the/config/file.yml
        
