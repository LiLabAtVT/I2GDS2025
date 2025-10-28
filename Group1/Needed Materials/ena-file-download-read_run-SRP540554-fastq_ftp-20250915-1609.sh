#!/bin/bash
#SBATCH -J g1-16s
#SBATCH --account=leaph
#SBATCH --partition=normal_q
#SBATCH --time=0-20:00:00
#SBATCH --mem=20G 
#SBATCH --mail-user=chanifah23@vt.edu
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR311/072/SRR31100672/SRR31100672_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR311/075/SRR31100675/SRR31100675_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR311/079/SRR31100679/SRR31100679_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR311/069/SRR31100669/SRR31100669_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR311/076/SRR31100676/SRR31100676_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR311/078/SRR31100678/SRR31100678_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR311/069/SRR31100669/SRR31100669_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR311/070/SRR31100670/SRR31100670_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR311/071/SRR31100671/SRR31100671_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR311/075/SRR31100675/SRR31100675_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR311/077/SRR31100677/SRR31100677_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR311/077/SRR31100677/SRR31100677_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR311/079/SRR31100679/SRR31100679_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR311/076/SRR31100676/SRR31100676_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR311/073/SRR31100673/SRR31100673_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR311/074/SRR31100674/SRR31100674_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR311/070/SRR31100670/SRR31100670_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR311/071/SRR31100671/SRR31100671_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR311/068/SRR31100668/SRR31100668_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR311/072/SRR31100672/SRR31100672_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR311/073/SRR31100673/SRR31100673_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR311/074/SRR31100674/SRR31100674_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR311/068/SRR31100668/SRR31100668_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR311/078/SRR31100678/SRR31100678_1.fastq.gz
