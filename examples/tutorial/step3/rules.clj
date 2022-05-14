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

