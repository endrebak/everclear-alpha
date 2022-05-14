(ns rules
  (:require [clojure.set :as set]
            [everclear.state.state :refer [config]]
            [everclear.parse.create-user-ns :refer [defrule]]))

(def hooo "config")

(defn add-a [w]
      (zipmap (keys w) (map (fn [v] (str v "a")) (vals w))))

(defrule bwa-map
  "Map DNA sequences against a reference genome with BWA."
  {:wildcards [:sample :genome]
   :externals  [:genome :fastq]
   :output    "bwa-map.bam"
   :resources {:cpu 1}
   :params    {:rg "@RG\\tID:{{wildcards.sample}}\\tSM:{{wildcards.sample}}"}
   :shell     "bwa mem -R '{{params.rg}}' -t {{resources.cpu}} {{externals.genome}} {{externals.fastq}} | samtools view -Sb - > {{output.0}}"})

(defrule samtools-sort
  "Sort the bams."
  {:wildcards [:sample :genome]
   :input     "bwa-map.bam"
   :input-fns {0 add-a 1 (fn [] (println "hello sailor!"))}
   :output    ["bam/sorted.bam"]
   :publish   "{{genome}}/{{sample}}/sorted.bam"
   :shell     "samtools sort -T {{wildcards.sample}} -O bam {{input.0}} > {{publish.0}}"})

;; (defrule samtools-index
;;   "Index read alignments for random access."
;;   {:wildcards [:sample :genome]
;;    :input     "bam/sorted.bam"
;;    :output    "bam/sorted.bam.bai"
;;    ;; TODO must ensure output and publish has the same keys!
;;    :publish   "{{genome}}/{{sample}}/sorted.bam.bai"
;;    :shell     "samtools index {{input.0}}"})
;;
;; (defrule bcftools-call
;;   "Aggregate mapped reads from all samples and jointly call genomic variants."
;;   {:input     {:sorted "bam/sorted.bam" :index "bam/sorted.bam.bai"}
;;    :output    "all.vcf"
;;    :wildcards [:genome]
;;    :externals  [:genome]
;;    :shell     "samtools mpileup -g -f {{externals.genome}} {{input.sorted}} | bcftools call -mv - > {{output.0}}"})
;;
;; (defrule plot-quals
;;   {:input     "all.vcf"
;;    :wildcards [:genome]
;;    :output    {:plot "quals.svg" :data "quals.tsv"}
;;    :script    "/Users/endrebakkenstovner/code/everclar_mirror/examples/snakemake/plot_quals.py"})

;; (defrule test-alone
;;   {:output "test2.txt"
;;    :shell "echo 'hello!!' > {{output.0}}"})

;; (defrule test-no-wildcards
;;   {:input "quals.svg"
;;    :output "test.txt"
;;    :shell "echo 'hello!!' > {{output.0}}"})

;; (defrule test-another-no-wildcards
;;   {:input "test.txt"
;;    :output "hiya.txt"
;;    :shell "echo 'hiya!!' > {{output.0}}"})

;; (defrule test-all
;;   {:input "hiya.txt"
;;    :output "done.txt"
;;    :wildcards [:genome :sample]
;;    :shell "cat {{input.0}} > {{output.0}}"})


;; (defrule create-genome
;;   {:output "genome.txt"
;;    :wildcards [:genome]
;;    :expand {0 [:sample :genome]}
;;    :shell "echo genome > {{output.0}}"})

;; (defrule split
;;   {:input "genome.txt"
;;    :output "chromo.txt"
;;    :publish "{{genome}}/{{sample}}/chromo.txt"
;;    :wildcards [:genome :sample]
;;    :same-job [:sample]
;;    :shell "echo 'hi' > {{output}}"})


;; (defrule use-multiple
;;   {:input "chromo2.txt"
;;    :wildcards [:genome]
;;    :output "multi.txt"
;;    :shell "echo 'multi' > {{output.0}}"})
