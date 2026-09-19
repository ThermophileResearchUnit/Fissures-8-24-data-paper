## read trimming (bbduk v39.52)
parallel -j 10 -a samples.list 'bbduk.sh ktrim=r ordered minlen=50 mink=11 tbo rcomp=f k=21 ow=t ftm=5 zl=4 qtrim=rl trimq=20 in1=raw/{}_R1_001.fastq.gz in2=raw/{}_R2_001.fastq.gz ref=~/adapters/adapters.fa out1=trimmed/{}_R1.trimmed.fastq.gz out2=trimmed/{}_R2.trimmed.fastq.gz'

## metaSPAdes (v3.15.4) assembly sbatch script
#!/bin/bash
#SBATCH -o tr.%j.out
#SBATCH -e tr.%j.err
#SBATCH -D .
#SBATCH --mem-per-cpu=8G
#SBATCH --time=1-0:00:00
#SBATCH --cpus-per-task=24

for i in `cat samples.list`;do
    spades.py \
        --meta \
        -t 24 \
        -m 192 \
        --pe1-1 trimmed/${i}_R1.trimmed.fastq.gz \
        --pe1-2 trimmed/${i}_R2.trimmed.fastq.gz \
        -k 21,33,55,77,101 \
        -o assemblies/metaspades_${i}
done

## binning using Metabat2 (v2.15) sbatch script

#!/bin/bash
#SBATCH -o bn.%j.out
#SBATCH -e bn.%j.err
#SBATCH -D .
#SBATCH --mem-per-cpu=2G
#SBATCH --time=1-0:00:00
#SBATCH --cpus-per-task=16
#SBATCH -J binning

for i in `cat samples.list`;do
    cd assemblies/metaspades_${i}
        seqtk comp contigs.fasta | awk '{if($2 >= 1000) print $1}' > contigs_gt1kb.list
	seqtk subseq contigs.fasta contigs_gt1kb.list > contigs_gt1kb.fasta
	mkdir binning
	cd binning
	bbwrap.sh \
	    ref=../contigs_gt1kb.fasta \
	    in1=../../trimmed/${i}_R1.trimmed.fastq.gz \
	    in2=../../trimmed/${i}_R2.trimmed.fastq.gz \
	    out=${i}.bam \
	    t=24 \
	    nodisk
        samtools sort -O bam -@ 16 ${i}.bam > ${i}.sorted.bam
	rm ${i}.bam
	jgi_summarize_bam_contig_depths \
	    --outputDepth depth_min1500.txt \
	    --pairedContigs paired_min1500.txt \
	    --minContigLength 1500 \
	    --minContigDepth 2 *.sorted.bam
        metabat2 -i ../contigs_gt1kb.fasta -a depth_min1500.txt -o bins/bin -t 16
    cd ../../../
done


## CheckM (v1.2.4) sbatch script

#!/bin/bash
#SBATCH -o cm.%j.out
#SBATCH -e cm.%j.err
#SBATCH -D .
#SBATCH -p large-gpu
#SBATCH --cpus-per-task=32
#SBATCH --mem=100G
#SBATCH --time=1-0:00:00
#SBATCH -J checkm


checkm lineage_wf -f CheckM.txt -t 32 -x fa mags_found checkm_out

## GTDB-Tk (v2.1.0) sbatch script

#!/bin/bash
#SBATCH -o gt.%j.out
#SBATCH -e gt.%j.err
#SBATCH -D .
#SBATCH -p viz
#SBATCH --cpus-per-task=32
#SBATCH --mem=500G
#SBATCH --time=2-0:00:00
#SBATCH -J gtdbtk

gtdbtk classify_wf --genome_dir all_bins --out_dir gtdb-tk_out --cpus 30 --pplacer_cpus 1 -x fa

