full workflow to sequence RNA from Oxford Nanopore Flowcells using Dorado/Minimap2 and running differential analysis using nanopolish and nanocompore.
**Currently set for RNA004 Datasets on experimental models**
**recently updated to include m6anet and xpore workflows for RNA004**

**https://github.com/nanoporetech/dorado**

run the following command for the basecall step, ensure the correct model and directory that data is pulling from.
change to the correct working directory (use cd .. to go up a level, cd /path/to/file to go down)

    /path/to/updated/bin/dorado \
    basecaller \
    /path/to/desired/version/dorado-0.8.2/dorado-0.8.2-linux-x64/models/rna004_130bps_sup@v5.1.0 \
    -r "/path/to/nanopore/rawdata/pod5/" \
    --min-qscore 10 #statistical test to rid bad basecalls \
    > /path/to/output/file.bam

alignment can be done with dorado, but minimap2 is preferred for downstream analysis purposes.
marginAlign has built in minimap, and provides more statistical data (marginalign will need to be run on python2.7)
**https://github.com/benedictpaten/marginAlign**

    -ax splice -uf -k14 "/path/to/reference/genome.fa"     "/path/to/basecalled.fastq" \
    > "path/to/minimap/aligned.sam"

**for -ax splice, it uses spliced alignment. Change to -ax map-ont if an unspliced alignment is needed (for m6anet and Xpore)**

prep the data to run statistical tests
       
        samtools view -h \
        /path/to/minimap/aligned.sam \
        | awk '$10 != "*"'  | samtools view -b -o             /path/to/filtered/file.bam 
  
To run marginStats you will need Python2.7.15 
create a virtual environment/conda with this before running.

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
        
With marginstats output, run substitutionplot.R and Identity.R scripts

With Dorado completed, use the following to move from unreadable bam into a readable transcript (AUCGs)

     samtools mpileup -uf /path/to/reference/genome.fa \
     -Q 10 -q 10 \
     /path/to/the/aligned.bam | bcftools call -c -o /path/to/variants/file.vcf

zip up the vcf file

        bgzip -c /path/to/variants/file.vcf/ > /path/to/variants.vcf.gz

to create consensus sequence, 

    bcftools index /path/to/variants.vcf.gz/

    bcftools consensus -f /path/to/reference/genome.fa \
    /path/to/variants.vcf.gz/ \
    > /path/to/output/consensus/sequence.fq

Last step is converting from T to U since this focuses on RNA (bam doesn't write Uracil, only Thymine)

      sed 's/T/U/g' /path/to/output/consensus/sequence.fq > /path/to/output/rna/consensus/sequence.fq

Next steps are to take the Dorado files and run using nanopolish to index and eventalign for modifications.

**https://github.com/hasindu2008/f5c**
**f5c is built to be the updated version of nanopolish that allows RNA004 modification detection**

        source /path/to/conda/env.sh
        conda activate $ENV_PATH
        conda new (set up nanopolish environment for dependencies)

Start the f5c nanopolish with initial prep,

    samtools sort -O bam -T "label" -o /path/to/aligned/sorted.bam \
    /path/to/aligned/file/wanting/to/be/sorted.bam

    samtools index /path/to/aligned/and/sorted.bam
    also make sure to remove the secondary and supplmentary reads ("filter" them)

last, install pod5 to convert from pod5 to fast5 (f5c only reads the old fast5, not pod5)

    pod5 convert to_fast5 "/path/to/the/pod5/files/pod5/" \
    --output "/path/to/new/directory/for/fast5/"
    
With library prep done, perform f5c index on the sample to rewrite the original ionic current from the flowcell for differential analysis

    /path/to/f5c \
    index -d "/path/to/the/file/fast5/flowdata.fast5 \
    /path/to/fastq/basecalled.fastq

**f5c also allows for m5C detection, but false positive rate seems to be higher than preferred**
**this is additional but can be skipped in workflow**
    /path/to/f5c \
    call-methylation --min-mapq 20 \
    -r /path/to/basecalled.fastq \
    -b /path/to/filtered/and/sorted.bam \ 
    -g /path/to/reference/genome.fa \
    --pore rna004 > /path/to/output_results.tsv

    /path/to/f5c \
    meth-freq -i /path/to/output_results.tsv > /path/to/output_frequency.tsv

**use f5c to run eventalign on the rna004 files**

    /path/to/f5c eventalign --rna \
    -b /path/to/aligned/and/sorted.bam \
    -r /path/to/basecalled.fastq \
    -g /path/to/reference.fa \
    --kmer-model /path/to/F5C/experimentalrna004/5mer.model \
    --signal-index --scale-events > /path/to/output/f5cevent.txt

**once complete the event file can be run with m6anet**
**https://github.com/GoekeLab/m6anet**

        m6anet dataprep --eventalign /path/to/event.txt \
        --out-dir /path/to/output/directory/ --n_processes 4 

        m6anet inference --input_dir /path/to/input/from/dataprep/ \
        --out-dir /path/to/m6anet/inference/directory/ \
        --pretrained_model HEK293T_RNA004 --n_processes 4 --num_iterations 1000

**using the same event.txt file from f5c, run Xpore**
**https://github.com/GoekeLab/xpore**

        xpore dataprep --eventalign /path/to/event.txt \
        --out_dir /path/to/output/xpore/directory/

To run the differential expression you will need a .yml config file, heres an example
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
        
