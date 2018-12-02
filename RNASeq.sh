global:
  - indir: data/processed
  - outdir: data/analysis_LV
  - hisat_align: 'data/analysis_LV/{$sample}/samtools_sort/{$sample}_hisat2.sorted.bam'
  - sample_rule: (.*LV.*)$
  - by_sample_outdir: '1'
  - find_by_dir: '1'
  - analysis_dir: data/analysis_LV/
  - REFERENCE: '/scratch/Reference_Genomes/Public/Vertebrate_other/Danio_rerio/ENSEMBL_release_84/Danio_rerio.GRCz10.dna.toplevel'
  - ANNOTATION: '/scratch/Reference_Genomes/Public/Vertebrate_other/Danio_rerio/ENSEMBL_release_84/Danio_rerio.GRCz10.84.gtf'
  - HPC:
      - module: gencore gencore_dev gencore_rnaseq
      - partition: serial
      - commands_per_node: '1'
      - cpus_per_task: '1'
      - account: 'gencore'
rules:
  - hisat2:
      local:
        - INPUT: '{$self->indir}/trimmomatic'
        - OUTPUT: '{$self->outdir}/{$sample}_hisat2.sam'
        - HPC:
            - walltime: '12:00:00'
            - mem: 60GB
            - cpus_per_task: '12'
      process: |
        #TASK tags={$sample}
        hisat2 \
        -x {$self->REFERENCE} \
        --dta \
        -p 12 \
        -1 {$self->INPUT}/{$sample}_read1_trimmomatic_1PE.fastq.gz \
        -2 {$self->INPUT}/{$sample}_read2_trimmomatic_2PE.fastq.gz \
        -S {$self->OUTPUT}
  - samtools_view:
      local:
        - INPUT: '{$self->analysis_dir}/{$sample}/hisat2/{$sample}_hisat2.sam'
        - OUTPUT: '{$self->outdir}/{$sample}_hisat2.bam'
        - HPC:
            - deps: hisat2
            - walltime: '04:00:00'
            - mem: 60GB
      process: |
        #TASK tags={$sample}
        samtools view -Su {$self->INPUT} > {$self->OUTPUT}
  - samtools_sort:
      local:
        - INPUT: '{$self->analysis_dir}/{$sample}/samtools_view/{$sample}_hisat2.bam'
        - OUTPUT: '{$self->outdir}/{$sample}_hisat2.sorted.bam'
        - HPC:
            - deps: samtools_view
            - walltime: '04:00:00'
            - mem: 60GB
            - cpus_per_task: '12'
      process: |
        #TASK tags={$sample}
        samtools sort -@ 12 -o {$self->OUTPUT} {$self->INPUT}
  - stringtie_1:
      local:
        - INPUT: '{$self->hisat_align}'
        - OUTPUT: '{$self->outdir}/{$sample}_stringtie1.gtf'
        - HPC:
            - deps: samtools_sort
            - walltime: '08:00:00'
            - mem: 60GB
            - cpus_per_task: '12'
      process: |
        #TASK tags={$sample}
        stringtie \
        {$self->INPUT} \
        -G {$self->ANNOTATION} \
        -l {$sample} \
        -o {$self->OUTPUT} \
        -p 12
  - stringtie_merge:
      local:
        - override_process: '1'
        - outdir: 'data/analysis_LV'
        - OUTPUT: '{$self->outdir}/stringtie_merge.gtf'
        - HPC:
            - deps: stringtie_1
            - walltime: '08:00:00'
            - mem: 60GB
      process: |
        #TASK tags={$sample}
        ls -1 {$self->outdir}/*/stringtie_1/*.gtf > {$self->outdir}/transcripts_list.txt && \
        stringtie --merge \
        -G {$self->ANNOTATION} \
        -o {$self->OUTPUT} \
        {$self->outdir}/transcripts_list.txt
  - stringtie_2:
      local:
        - INPUT: '{$self->hisat_align}'
        - HPC:
            - deps: stringtie_merge
            - walltime: '08:00:00'
            - mem: 60GB
            - cpus_per_task: '12'
      process: |
        #TASK tags={$sample}
        stringtie \
        {$self->INPUT} \
        -eB \
        -G {$self->analysis_dir}/stringtie_merge.gtf  \
        -l {$sample} \
        -o {$self->outdir}/{$sample}_stringtie_2.gtf \
        -p 12
  - htseq_count:
      local:
        - INPUT: "{$self->hisat_align}"
        - OUTPUT: "{$self->outdir}/{$sample}_rawCounts.txt"
        - HPC:
            - deps: samtools_sort
            - walltime: '18:00:00'
            - mem: '35GB'
            - cpus_per_task: '1'
      process: |
        #TASK tags={$sample}
        htseq-count -f bam -s no -t exon \
        -i gene_id \
        {$self->INPUT} \
        {$self->ANNOTATION} > {$self->OUTPUT} && \
        sed -i '/^__.*/d'  {$self->OUTPUT}
