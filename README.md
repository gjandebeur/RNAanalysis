full workflow to sequence RNA from oxford nanopore using dorado and run differential analysis using nanopolish and nanocompore

download the up to date dorado from the dorado github, pull onto OSCER using wget.

**https://github.com/nanoporetech/dorado**

run the following command for the basecall step, ensure correct model and directory that data is pulling from.
change to the correct working directory (use cd .. to go up a level, cd /path/to/file to go down)

    /path/to/updated/bin/dorado \
    basecaller \
    /path/to/desired/version/dorado-0.8.2/dorado-0.8.2-linux-x64/models/rna004_130bps_sup@v5.1.0 \
    -r "/path/to/nanopore/rawdata/pod5/" \
    --min-qscore 10 #statistical test to rid bad basecalls \
    > /path/to/output/file.bam


        #dependencies potentially needed (input into virtual env for ease)
        module load GCC
        module load PyTorch
        module load FlexiBLAS/3.3.1-GCC-12.3.0  
        module load FFmpeg/4.4.2-GCCcore-11.3.0 
        module load HTSlib
        module load protobuf

next demux samples to label unique reads with barcode

    /path/to/dorado/bin/dorado \
    demux \
    --kit-name SQK-RPB004 --output-dir /path/to/output/directory/demux/ \ 
    "/path/to/basecalled/file.bam"
    if [ $? -ne 0 ]; then
    echo "Dorado demux failed."
    exit 1
    fi
    echo "Dorado demux completed!"

next align the basecalled/demux'd data with the reference genome from gencode

    /path/to/dorado/bin/dorado \
    aligner \
    "/path/to/reference/gencome.fa" \
    "/path/to/demuxed/file/directory/pick/largest/one.bam" > /output/path/aligned.bam

next run dorado summary to create summary file for downstream analysis.

        dorado summary \
        /path/to/aligned/output.bam \
        > /path/to/output/summary.txt

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

**https://github.com/jts/nanopolish**

        source /path/to/conda/env.sh
        conda activate $ENV_PATH
        conda new (set up nanopolish environment for dependencies)

Start the nanopolish with initial prep,

    samtools sort -O bam -T "label" -o /path/to/aligned/sorted.bam \
    /path/to/aligned/file/wanting/to/be/sorted.bam

    samtools index /path/to/aligned/and/sorted.bam

    samtools fasta -@ 4 /path/to/original/dorado/basecalled.bam \
    > /path/to/new/basecalled.fa

last, install pod5 to convert from pod5 to fast5 (nanopolish only reads the old fast5 currently)

    pod5 convert to_fast5 "/path/to/the/pod5/files/pod5/" \
    --output "/path/to/new/directory/for/fast5/"
    
With library prep done, perform nanopolish index on the sample to rewrite the original ionic current from flowcell for differential analysis

    nanopolish index -s /path/to/summary/file/summary.txt \
    -d /path/to/directory/with/rawdata/fast5/ \
    /path/to/the/basecalled/original/file/basecalled.fasta

After nanopolish indexes (takes a while) it will generate .fai and .gzi files for the .fasta,

    nanopolish eventalign --reads /path/to/created/basecalled.fasta \
    --bam /path/to/the/aligned/and/sorted.bam \
    --genome /path/to/reference/genome/transcript.fa \
    --print-read-names --scale-events --samples --min-mapping-quality 10 >
    /path/to/new/file/eventalign.txt

Install nanocompore & dependencies into a virtual environment or conda

**https://github.com/tleonardi/nanocompore**

additional assistance and syntax

**https://nanocompore.rna.rocks/data_preparation/**

    nanocompore eventalign_collapse -t 6 -i \
    /path/to/new/file/eventalign.txt \
    -o /path/to/output/collapsed.eventalign.tsv

Once the following is complete for different variables, ideally with replicates, run sampcomp

    nanocompore sampcomp -1 /path/to/first/collapsed.eventalign.tsv \
    -2 /path/to/other/variable/collapsed.eventalign.tsv \
    --fasta /path/to/reference/genome/transcript.fa \
    --label1 "label" \
    --label2 "2ndlabel" \
    --outpath /path/to/directory/for/nanocompore/output/ \
    --overwrite #this forces it to write over the directory! \
    --min_coverage 10 


    
