(defrule bwa-map
  "Map DNA sequences against a reference genome with BWA."
  {:externals  [:fastq :genome]
   :output    "bwa-map.bam"
   :shell     "bwa mem {{externals.genome}} {{externals.fastq}} |
               samtools view -Sb - > {{output.0}}"})