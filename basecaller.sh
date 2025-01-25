#!/bin/sh
#SBATCH --partition=rnafold
#SBATCH --job-name=hbec24r1a_sequencing
#SBATCH --output=hbec24r1a_sequencing_output.txt
#SBATCH --error=hbec24r1a_sequencing_debug.txt
#SBATCH --gpus-per-node=1
#SBATCH --time=48:00:00           
#SBATCH --mem=32G         # good memory allocation?

/ourdisk/hpc/rnafold/gjandebeur/dont_archive/software/dorado-0.8.2/dorado-0.8.2-linux-x64/bin/dorado \
basecaller \
/ourdisk/hpc/rnafold/gjandebeur/dont_archive/software/dorado-0.8.2/dorado-0.8.2-linux-x64/models/rna004_130bps_sup@v5.1.0 \
-r "/ourdisk/hpc/rnafold/sschroed/dont_archive/HBEC24r1a/20241025_2020_P2S-00590-A_PAY65891_bb185c9d/pod5/" \
--min-qscore 10 \
> /ourdisk/hpc/rnafold/gjandebeur/dont_archive/hbec_24r1a/rawdata/hbec24r1a_basecalled.bam

/ourdisk/hpc/rnafold/gjandebeur/dont_archive/software/dorado-0.8.2/dorado-0.8.2-linux-x64/bin/dorado \
basecaller \
/ourdisk/hpc/rnafold/gjandebeur/dont_archive/software/dorado-0.8.2/dorado-0.8.2-linux-x64/models/rna004_130bps_sup@v5.1.0 \
-r "/ourdisk/hpc/rnafold/sschroed/dont_archive/HBEC24r1a/20241025_2020_P2S-00590-A_PAY65891_bb185c9d/pod5/" \
--modified-bases m6A_DRACH \
--min-qscore 10 \
> /ourdisk/hpc/rnafold/gjandebeur/dont_archive/hbec_24r1a/rawdata/hbec24r1a_basecalled_m6ADRACH.bam

echo "m6A_DRACH sequenced!"

# the modified bases past this point are high false positive, don't blindly trust results

echo "starting m5C sequencing next!"

/ourdisk/hpc/rnafold/gjandebeur/dont_archive/software/dorado-0.8.2/dorado-0.8.2-linux-x64/bin/dorado \
basecaller \
/ourdisk/hpc/rnafold/gjandebeur/dont_archive/software/dorado-0.8.2/dorado-0.8.2-linux-x64/models/rna004_130bps_sup@v5.1.0 \
-r "/ourdisk/hpc/rnafold/sschroed/dont_archive/HBEC24r1a/20241025_2020_P2S-00590-A_PAY65891_bb185c9d/pod5/" \
--modified-bases m5C \
> /ourdisk/hpc/rnafold/gjandebeur/dont_archive/hbec_24r1a/rawdata/hbec24r1a_basecalled_m5C.bam

echo "m5C sequencing finished"


echo "starting pseU basecaller"
/ourdisk/hpc/rnafold/gjandebeur/dont_archive/software/dorado-0.8.2/dorado-0.8.2-linux-x64/bin/dorado \
basecaller \
/ourdisk/hpc/rnafold/gjandebeur/dont_archive/software/dorado-0.8.2/dorado-0.8.2-linux-x64/models/rna004_130bps_sup@v5.1.0 \
-r "/ourdisk/hpc/rnafold/sschroed/dont_archive/HBEC24r1a/20241025_2020_P2S-00590-A_PAY65891_bb185c9d/pod5/" \
--modified-bases pseU \
 > /ourdisk/hpc/rnafold/gjandebeur/dont_archive/hbec_24r1a/rawdata/hbec24r1a_basecalled_pseU.bam


echo "starting inosine modification sequencing"

/ourdisk/hpc/rnafold/gjandebeur/dont_archive/software/dorado-0.8.2/dorado-0.8.2-linux-x64/bin/dorado \
basecaller \
/ourdisk/hpc/rnafold/gjandebeur/dont_archive/software/dorado-0.8.2/dorado-0.8.2-linux-x64/models/rna004_130bps_sup@v5.1.0 \
-r "/ourdisk/hpc/rnafold/sschroed/dont_archive/HBEC24r1a/20241025_2020_P2S-00590-A_PAY65891_bb185c9d/pod5/" \
--modified-bases inosine_m6A \
 > /ourdisk/hpc/rnafold/gjandebeur/dont_archive/hbec_24r1a/rawdata/hbec24r1a_basecalled_inosine.bam

echo "inosine sequencing finished"

echo "program complete!"
