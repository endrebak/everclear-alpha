(defrule bwa-map
  "Map DNA sequences against a reference genome with BWA."
  {:externals  [:fastq :genome]
   :wildcards [:sample]
   :output    "bwa-map.bam"
   :shell     "bwa mem {{externals.genome}} {{externals.fastq}} |
               samtools view -Sb - > {{output.0}}"})

(defrule samtools-sort
  "Sort the bams."
  {:wildcards [:sample]
   :input     "bwa-map.bam"
   :output    ["bam/sorted.bam"]
   :publish   "{{sample}}/sorted.bam"
   :shell     "samtools sort -T {{wildcards.sample}} -O bam {{input.0}} > {{publish.0}}"})

(defrule samtools-index
         "Index read alignments for random access."
         {:wildcards [:sample :genome]
          :input     "bam/sorted.bam"
          :output    "bam/sorted.bam.bai"
          :publish   "{{sample}}/sorted.bam.bai"
          :shell     "samtools index {{input.0}}"})

(defrule bcftools-call
         "Aggregate mapped reads from all samples and jointly call genomic variants."
         {:input     {:sorted "bam/sorted.bam" :index "bam/sorted.bam.bai"}
          :output    "all.vcf"
          :wildcards [:genome]
          :externals  [:genome]
          :shell     "samtools mpileup -g -f {{externals.genome}} {{input.sorted}} | bcftools call -mv - > {{output.0}}"})


(defrule plot-quals
  {:input     "all.vcf"
   :wildcards [:genome]
   :expand    {:data [:sample]}
   :output    {:plot "quals.svg" :data "quals.tsv"}
   :script    "/Users/endrebakkenstovner/code/everclar_mirror/examples/snakemake/plot_quals.py"})
