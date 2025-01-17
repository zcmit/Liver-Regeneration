### macs2 calls H3K27me3 broad peaks
macs2 callpeak -t /scratch/cz21/Chip_seq_Mouse/data/analysis/Sample_WTB-H3K27me3/samtools/Sample_WTB-H3K27me3_aligned.bam -c /scratch/cz21/Chip_seq_Mouse/data/analysis/Sample_Input/samtools/Sample_Input_aligned.bam -g mm -p 1e-3 --nomodel --shift 37 --extsize 73 -n WT_H3K27me3_broad
### SICER calls H3K27me3 broad peaks
sh /SICER1.1/SICER/SICER.sh . WT_K27_S26_mdup_Addchr.bed Input_S23_mdup_Addchr.bed . mm10 1 1000 75 0.75 4000 .01


### macs2 calls H3K4me3 narrow peaks
#### 1st batch
macs2 callpeak -t /scratch/cz21/Chip_seq_Mouse/data/analysis/Sample_KOB-H3K4me3/samtools/Sample_KOB-H3K4me3_aligned.bam -g mm -n KO_H3K4me3_dup -c /scratch/cz21/Chip_seq_Mouse/data/analysis/Sample_Input/samtools/Sample_Input_aligned.bam -B
#### 2nd batch
macs2 callpeak -t $DAT/KO_K4_S27/picard/KO_K4_S27_mdup.withrg.csorted.cleaned.aligned.bam \
-c $DAT/Input_S23/picard/Input_S23_mdup.withrg.csorted.cleaned.aligned.bam \
-g mm \
--keep-dup all \
--outdir KO_K4_S27/macs2/narrowPeaks/ \
-n KO_K4_S27 \
-B \
--SPMR &

### macs2 calls H3K9me3 narrow peaks
macs2 callpeak -t $DAT/KO_K9_S28/picard/KO_K9_S28_mdup.withrg.csorted.cleaned.aligned.bam \
-c $DAT/Input_S23/picard/Input_S23_mdup.withrg.csorted.cleaned.aligned.bam \
-g mm \
--keep-dup all \
--outdir KO_K9_S28/macs2/narrowPeaks/ \
-n KO_K9_S28 \
-B \
--SPMR &
