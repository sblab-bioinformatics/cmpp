

## G4-ChIP-seq analysis
The detection of high confidence G4 structures in K562 has previously been described in [Mao et. al 2018](https://www.nature.com/articles/s41594-018-0131-8).
A more detailed description of our in-house G4 ChIP-seq pipeline can be found [here](https://github.com/sblab-bioinformatics/dna-g4-methylation-dnmt1). 
A bedfile containing BG4 ChIP-seq high confidence sites is available from GEO [GSE107690](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE107690)
Sites that have the biophysical potential to form G4s can also be downloaded from GEO [GSE110582]{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE110582}

## SMARCA4 ChIP-seq analysis
```
mkdir bam
mkdir clean_bam
mkdir macs2_output
mkdir tdf
mkdir bigWig
mkdir logs
mkdir analysis
mkdir analysis/confidence_peaks
mkdir fastq_trimmed
```


### adapter trimming for single-end sequencing

```
cd fastq

for f in *1.fq.gz ; do sbatch --mem 8G --time 120 -o $f.out -e $f.err -J cutadapt --wrap "~/.local/bin/cutadapt -m 10 -q 20 -O 3 -a AGATCGGAAGAGC -g GCTCTTCCGATCT -o ../fastq_trimmed/${f%%.H2FK2BGXG.s_1.r_1.fq.gz}.trimmed.fq.gz $f"; done

```


### alignment single-end to hg19
Align to hg19 reference genome.

Align:
```
cd ../fastq_trimmed

genome='~/reference_data/hg19/Sequence/Bowtie2Index/bowtie2-build/genome.fa'; for f in *.trimmed.fq.gz; do sbatch --mem 16G -o $f.bwa.out -e $f.bwa.err -J bwa --wrap "~/bwa mem -t 8 -M $genome $f | ~/bin/samtools view -Sb -F2308 -q 10 > ../bam/${f%%.trimmed.fq.gz}.hg19.bam"; done
```


### Sort bam files, mark and remove duplicates
```
for f in *.bam; do sbatch --mem 16G -o $f.out -e $f.err -J sort --wrap "~/bin/samtools sort -@ 8 $f > ${f%%.hg19.bam}.sort.bam"; done
```


Mark and remove duplicates
```
for f in *sort.bam; do sbatch --mem 8G -o $f.out -e $f.err -J dups --wrap "java -Xmx7g -jar B/picard/picard-2.14.0/picard.jar MarkDuplicates I=$f O=../clean_bam/${f%%.sort.bam}.nodup.bam M=${f%%.sort.bam}.md.txt AS=true REMOVE_DUPLICATES=true"; done
```

Some basic dup and size stats
```
grep 'Unknown Library' *txt > dups_info.txt

srun --mem 32G --pty /usr/bin/bash
cd ../clean_bam/
for f in *nodup.bam; do echo $f;echo $f >> all_counts && ~/bin/samtools view -c $f >> all_counts; done
exit
```

**RESULTS**
Duplication rate
```
K562_Rep1_Input_SLX-19726.NEBi708_NEBi508.bam.md.txt:Unknown Library	50500555	0	0	0	8068346	0	0	0.159767	
K562_Rep1_SMARCA4_SLX-19726.NEBi709_NEBi502.bam.md.txt:Unknown Library	64296462	0	0	0	13840513	0	0	0.215261	
K562_Rep2_Input_SLX-19727.NEBi701_NEBi501.md.txt:Unknown Library	33064335	0	0	0	1331326	0	0	0.040265	
K562_Rep2_SMARCA4_SLX-19727.NEBi707_NEBi507.md.txt:Unknown Library	53983605	0	0	0	4937001	0	0	0.091454	
K562_Rep3_Input_SLX-19727.NEBi702_NEBi502.md.txt:Unknown Library	84253366	0	0	0	4157543	0	0	0.049346	
K562_Rep3_SMARCA4_SLX-19727.NEBi708_NEBi508.md.txt:Unknown Library	44722139	0	0	0	3028800	0	0	0.067725	

```
Duplication rate a bit high rep1, but overall ok.


Read count per bam file:
```
K562_Rep1_Input_SLX-19726.NEBi708_NEBi508.bam 42432209
K562_Rep1_SMARCA4_SLX-19726.NEBi709_NEBi502.bam 50455949
K562_Rep2_Input_SLX-19727.NEBi701_NEBi501.bam 31733009
K562_Rep2_SMARCA4_SLX-19727.NEBi707_NEBi507.bam 49046604
K562_Rep3_Input_SLX-19727.NEBi702_NEBi502.bam 80095823
K562_Rep3_SMARCA4_SLX-19727.NEBi708_NEBi508.bam 41693339

```
>30 Mio reads in all cases. 


### create index bam files and generate tdf files for visual inspection
Generate bai files using samtools index
```
for f in *nodup.bam; do sbatch -o %j.$f.log --mem 16G -J index.$f --wrap "~/bin/samtools index -b $f"; done

```

Generate tdfs
```
for f in *nodup.bam; do sbatch -o %j.$f.log --mem 8G -J igv.$f --wrap "~/bin/igvtools/igvtools-2.3.91/igvtools count -w 15 $f ../tdf/${f%%.nodup.bam}.tdf ~/bin/igvtools/igvtools-2.3.91/genomes/hg19.chrom.sizes"; done
```

### create bigWig files
```
for FILE in *nodup.bam
do
bname=`basename $FILE .bam`
sbatch -o %j.$f.log -J BCoverage --mem=8G \
--wrap "bamCoverage -b $FILE \
-o ../bigWig/${FILE%%.nodup.bam}.bs50.bl.RPKM.bw \
--binSize 50 \
--blackListFileName ~/bin//reference_data/hg19/hg19.blacklist_merge1000nt.bed \
--numberOfProcessors max \
--normalizeUsing RPKM"
done
```



### call peaks
(using macs2 2.1.1.20160309)
Call peaks vs input file.

```
sbatch --mem 8G -o Rep1.macs2.out -e Rep1.macs2.err -J macs2 --wrap "macs2 callpeak --name ../macs2_output/K562_Rep1_SMARCA4_SLX-19726.NEBi709_NEBi502.nodup.q005.all -t K562_Rep1_SMARCA4_SLX-19726.NEBi709_NEBi502.bam -c K562_Rep1_Input_SLX-19726.NEBi708_NEBi508.bam --format=BAM --gsize 'hs' --bw=300 --qvalue 0.05"; done

sbatch --mem 8G -o Rep2.macs2.out -e Rep2.macs2.err -J macs2 --wrap "macs2 callpeak --name ../macs2_output/K562_Rep2_SMARCA4_SLX-19727.NEBi707_NEBi507.nodup.q005.all -t K562_Rep2_SMARCA4_SLX-19727.NEBi707_NEBi507.bam -c K562_Rep1_Input_SLX-19726.NEBi708_NEBi508.bam --format=BAM --gsize 'hs' --bw=300 --qvalue 0.05"; done

sbatch --mem 8G -o Rep3.macs2.out -e Rep3.macs2.err -J macs2 --wrap "macs2 callpeak --name ../macs2_output/K562_Rep3_SMARCA4_SLX-19727.NEBi708_NEBi508.nodup.q005.all -t K562_Rep3_SMARCA4_SLX-19727.NEBi708_NEBi508.bam -c K562_Rep3_Input_SLX-19727.NEBi702_NEBi502.bam --format=BAM --gsize 'hs' --bw=300 --qvalue 0.05"; done

```


#### Generate high confidence peaks
```
cd ../macs2_output

wc -l *narrowPeak
  21487 K562_Rep1_SMARCA4_SLX-19726.NEBi709_NEBi502.narrowPeak
  27356 K562_Rep2_SMARCA4_SLX-19727.NEBi707_NEBi507.narrowPeak
  40707 K562_Rep3_SMARCA4_SLX-19727.NEBi708_NEBi508.narrowPeak

mkdir confidence_peaks

multiIntersectBed -i *.narrowPeak > confidence_peaks/K562_SMARCA4_rep1-3_multintersect.bed

cd confidence_peaks

for ((X=1; X<=3; X++))
do
FILENAME="K562_SMARCA4_rep1-3_mult.$X""of3.bed"
awk -v VARIABLE="$X" '$4 >= VARIABLE'< K562_SMARCA4_rep1-3_multintersect.bed | sortBed -i - | mergeBed -i - > $FILENAME
echo $FILENAME
done


wc -l *.bed

42367 K562_SMARCA4_rep1-3_mult.1of3.bed
28265 K562_SMARCA4_rep1-3_mult.2of3.bed
18889 K562_SMARCA4_rep1-3_mult.3of3.bed
```

Plot overlap of individual replicates
```
sbatch -J Intervene -o pw.log --mem 2g --wrap " intervene pairwise -i \
K562_Rep1_SMARCA4_SLX-19726.NEBi709_NEBi502.narrowPeak \
K562_Rep2_SMARCA4_SLX-19727.NEBi707_NEBi507.narrowPeak \
K562_Rep3_SMARCA4_SLX-19727.NEBi708_NEBi508.narrowPeak \
--names \
Rep1,Rep2,Rep3 \
--figsize 7 4 \
--compute frac --htype tribar \
-o SMARCA4_ChIP_replciates_PW_triang"

sbatch -J Venn -o pw.log --mem 2g --wrap " intervene venn -i \
K562_Rep1_SMARCA4_SLX-19726.NEBi709_NEBi502.narrowPeak \
K562_Rep2_SMARCA4_SLX-19727.NEBi707_NEBi507.narrowPeak \
K562_Rep3_SMARCA4_SLX-19727.NEBi708_NEBi508.narrowPeak \
--names \
Rep1,Rep2,Rep3 \
-o SAMRCA4_venn"
```
Focus on peaks present in at least 2 of 3 replicates.

#### Overlap with BG4 ChIP

```
#with BG4 signal
bedtools intersect -a K562_SMARCA4_rep1-3_mult.2of3.bed -b /reference_data/Quadruplex/GSE107690_K562_High_confidence_peaks.bed > SMARCA4_with_BG4.bed
wc -l SMARCA4_with_BG4.bed 
7565

#without BG4 signal
bedtools intersect -v -a K562_SMARCA4_rep1-3_mult.2of3.bed -b /scratchb/sblab/spiege01/ENCODE_K562/reference_data/Quadruplex/GSE107690_K562_High_confidence_peaks.bed > SMARCA4_no_BG4.bed
 wc -l SMARCA4_no_BG4.bed
```




### SMARCA4 ChIP Singal in G4 ChIP
```
cd ../../bigWig


sbatch -o ${f%%.bs50.bl.RPKM.bw}.log -J DT_Profile --mem=10G \
--wrap "computeMatrix reference-point \
--referencePoint center \
--missingDataAsZero \
-b 400 -a 400 \
-S \
K562_Rep1_SMARCA4_SLX-19726.NEBi709_NEBi502.bw \
K562_Rep2_SMARCA4_SLX-19727.NEBi707_NEBi507.bw \
K562_Rep3_SMARCA4_SLX-19727.NEBi708_NEBi508.bw \
K562_Rep1_Input_SLX-19726.NEBi708_NEBi508.bw \
K562_Rep2_Input_SLX-19727.NEBi701_NEBi501.bw \
K562_Rep3_Input_SLX-19727.NEBi702_NEBi502.bw \ 
-R \
/reference_data/Quadruplex/GSE107690_K562_High_confidence_peaks.bed \
/reference_data/Quadruplex/GSE110582_G4-Seq_K_PDS_plus_and_minus.bed.gz \
--skipZeros \
-o SMARCA4_in_stranded_BG4_5o8_vs_OQS.mat.gz" 


sbatch -o profile.log -J DT_Profile --mem=5000 \
--wrap "plotProfile -m SMARCA4_in_stranded_BG4_5o8_vs_OQS.mat.gz \
-out SMARCA4_in_stranded_BG4_5o8_vs_OQS.pdf \
--refPointLabel PeakCenter \
--regionsLabel BG4 \
--samplesLabel Rep_1 Rep_2 Rep_3 Input_1 Input_2 Input_3 \
--perGroup \
--dpi 300 \
--plotHeight 10 \
--plotWidth 10 \
--colors black blue red grey lightblue orange" 

```

### PAVIS analysis 
The distribution of SMARCA4 and BG4 ChIP confidence sites across functional genomic regions was explored using the [PAVIS online tool](https://manticore.niehs.nih.gov/pavis2/).


### MEME ChIP analysis
Motif analysis for SMARCA4 confidence binding sites with/without BG4 was performed using the [MEME ChIP online tool](https://meme-suite.org/meme/tools/meme-chip) 

*MEME ChIP parameters*:
- count of motifs: 5
- width 6-40











