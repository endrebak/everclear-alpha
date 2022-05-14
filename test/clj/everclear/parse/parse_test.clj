(ns everclear.parse.parse-test
  (:require
    [clojure.java.io :as io]
    [clojure.test :refer :all]
    [everclear.parse.parse :as parse]))

(def rulestring "(defrule samtools-index
         \"Index read alignments for random access.\"
         {:wildcards [:sample :genome]
          :input     \"bam/sorted.bam\"
          :output    \"bam/sorted.bam.bai\"
          :publish   \"{{genome}}/{{sample}}/sorted.bam.bai\"
          :shell     \"samtools index {{input.0}}\"})")

(with-open [rulefile (io/as-file (io/reader (char-array rulestring)))]
  (is (parse/->forms rulefile) nil))
