#!/usr/bin/env bb

(require '[clojure.tools.cli :refer [parse-opts]])
(require '[clojure.java.shell :refer [sh]])
(require '[clojure.java.io :refer [make-parents]])

(def cli-options
  ;; An option with a required argument
  [["-f" "--file FILE" "SAM/BAM file to summarize."]
   ["-o" "--output OUTPUT" "Where to place the summary"]
   ["-h" "--help"]])

(def args (parse-opts *command-line-args* cli-options))

(def file (get-in args [:options :file]))
(def outfile (get-in args [:options :output]))

(def result (sh "samtools" "flagstat" "-O" "tsv" file))
(def output (:out result))

(if outfile
  (do
    (make-parents output)
    (spit outfile output))
  (pprint output))
