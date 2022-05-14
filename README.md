# everclear

generated using Luminus version "4.06"

FIXME

## Prerequisites

You will need [Leiningen][1] 2.0 or above installed.

[1]: https://github.com/technomancy/leiningen

(let [pb (ProcessBuilder. ["python"])
      p (.start pb)
      out (.getOutputStream p)
      in (.getInputStream p)]
  (spit out "print(\"hello from lua\"\n)")
  (slurp in))

## Running

``` clojure
;; cider-jack-in
user>
(require '[everclear.dag.helpers :as helpers :refer [file-to-map]])

(require '[everclear.parse :as parse])
(require '[selmer.parser :refer [render]])
(require '[clojure.java.io :as io :refer [file]])
(require '[clojure.walk :as walk])
(require '[com.stuartsierra.dependency :as dep])

(prefer-method clojure.pprint/simple-dispatch clojure.lang.IPersistentMap clojure.lang.IDeref)
(everclear.core/start-with-map {:config-file "examples/snakemake/config.edn"})

(everclear.dag.helpers/get-jobs-by-rule :plot-quals @everclear.state.state/jobhash->jobid @everclear.state.state/jobhash->job)

(alias 'st 'everclear.state.state)
(def js st/jobhash->jobstate)
(def j st/jobhash->job)

(def job (second (second @st/jobhash->job)))

(swap! st/config assoc :config-file "/Users/endrebakkenstovner/code/epigenome_new/config.edn")
(swap! st/config assoc :config-file "examples/snakemake/config.edn")

(def p (everclear.dag.helpers/->Path "/Users/endrebakkenstovner/everclear/snakemake-flow/results/7b3a6ddd30/bam/sorted.bam.bai"))
(java.nio.file.Files/exists p (into-array [java.nio.file.LinkOption/NOFOLLOW_LINKS]))
(java.nio.file.Files/exists p (make-array java.nio.file.LinkOption 0))

(require '[everclear.parse :as parse])
(require '[selmer.parser :refer [render]])
(require '[clojure.java.io :as io :refer [file]])
(require '[clojure.walk :as walk])
(require '[com.stuartsierra.dependency :as dep])

(prefer-method clojure.pprint/simple-dispatch clojure.lang.IPersistentMap clojure.lang.IDeref)
(everclear.core/start-with-map {:config-file "examples/snakemake/config.edn"})


(def config-file "examples/snakemake/config.edn")
(def config (everclear.dag.helpers/file-to-map config-file))
(def rule-file (config :rule-file))
(def wildcards-file (config :wildcards-file))
(def externals-file (config :externals-file))

(def rules (everclear.parse/ruleinfo-keyed-by-rulename (config :rule-file)))
(def rulename->ruleinfo rules)
(def ruledata (vals rules))
(def wildcards (everclear.dag.helpers/file-to-map (config :wildcards-file)))
(def external (everclear.dag.helpers/file-to-map (config :externals-file)))
(def base-path (config :base-path))

(def gs (everclear.graph/->rulegraph+jobgraph+filegraph (vals rules) wildcards))
(def rg (first gs))
(def jg (second gs))
(def fg (nth gs 2))
(def rulegraph rg)
(def jobgraph jg)
(def filegraph fg)

(def jobid->job (everclear.jobs/job-keyed-by-jobid jobgraph filegraph rulename->ruleinfo external base-path))

(everclear.dag.core/create-graphs&jobs rule-file wildcards-file externals-file base-path)

(def job->hash (everclear.hash/finalhash-keyed-by-job jobgraph rulename->ruleinfo))
(def jobids (keys job->hash))
(def jobid->outmap (everclear.jobs/outmap-keyed-by-jobid job->hash rulename->ruleinfo base-path))
(def jobfile->outentry (everclear.jobs/outentry-keyed-by-jobfile jobid->outmap rulename->ruleinfo))
(def jobid->publishmap (everclear.jobs/publishmap-keyed-by-jobid jobids rulename->ruleinfo))
(def jobfile->publishentry (everclear.jobs/publishentry-keyed-by-jobfile jobid->publishmap rulename->ruleinfo))
(def rule+file->inalias (everclear.jobs/inalias-keyed-by-rulename-and-file (vals rulename->ruleinfo)))
(def jobid->jobpaths (everclear.jobs/jobpaths-keyed-by-jobid jobids filegraph jobfile->outentry jobfile->publishentry rule+file->inalias))

(def jobid->job (job-keyed-by-jobid jobgraph filegraph rulename->ruleinfo externals base-path))

(def jobid [:bcftools-call {:genome "hg38"}])
(everclear.jobs/modify-outalias-to-inalias :bcftools-call (dep/immediate-dependencies filegraph ) )

(def producer+file (everclear.graph/->producer+file ruledata))
(def file->producer (everclear.graph/producer-keyed-by-file producer+file))
(def producer+consumer (everclear.graph/->producer+consumer
 file->producer
 (everclear.graph/->file+consumer ruledata)))
(def rulegraph (everclear.graph/->rulegraph producer+consumer))
(def rulename->ruleinfo (zipmap (map :rulename ruledata) ruledata))
(def pjob+cjob (everclear.graph/->pjob+cjob rulename->ruleinfo rulegraph wildcards))
(def jobgraph (everclear.graph/->jobgraph pjob+cjob))
(def producer->file (zipmap (vals file->producer) (keys file->producer)))
(def pjob+jf+cjob (everclear.graph/->pjob+jobfile+cjob pjob+cjob producer->file))
(def filegraph (everclear.graph/->filegraph pjob+jf+cjob))










(clojure.pprint/pprint (:dependencies jg) (clojure.java.io/writer "jobgraph.clj"))
(clojure.pprint/pprint (:dependencies fg) (clojure.java.io/writer "filegraph.clj"))
(clojure.pprint/pprint (dep/topo-sort fg) (clojure.java.io/writer "toposorted.clj"))

To start a web server for the application, run:

    lein run

## License

Copyright © 2021 FIXME
# everclear

## TODO



#### Display joboutput

Also want to display info about

When a job starts what its input-files are
When a job ends, how long time it took, what files it produced

#### Hashing
If hash of input-files, code and settings-map are the same, the file does not need to be recreated.
If hash hinges on input-files, then final output-path cannot be created before job is set to run.
Alternative: keep the hash that tells how the files were made
Dispatch-job should compute the hash?


https://github.com/arachne-framework/valuehash

(require '[babashka.process :as p :refer [process]]
         '[clojure.java.io :as io])
(def bwa-mem (process ["bash" "-c" "bwa mem -R '@RG\\tID:A\\tSM:A'  /Users/endrebakkenstovner/everclear/snakemake-flow/data/genome.fa /Users/endrebakkenstovner/everclear/snakemake-flow/data/samples/A.fastq | samtools view -Sb - > /Users/endrebakkenstovner/everclear/snakemake-flow/bwa-map/genome/hg19/sample/A/bwa-map.bam"] {:err :inherit
                           :shutdown p/destroy}))
(with-open [rdr (io/reader (:out bwa-mem))]
  (binding [*in* rdr]
    (loop []
      (let [line (read-line)]
        (println :line line)
      (when (not (nil? line))
        (recur))))))


(def ls (process ["ls" "-latrh"] {:err :inherit :shutdown p/destroy}))
(with-open [rdr (io/reader (:out ls))]
  (binding [*in* rdr]
    (loop []
      (let [line (read-line)]
        (println :line line)
      (when (not (nil? line))
        (recur))))))


(def proc (:proc (babashka.process/process ["ls"])))
(def cf (.thenApply (.onExit proc) (reify java.util.function.Function (apply [this p] (.exitValue p)))))
(.get cf)

#{[[:samtools-index {:genome "hg38", :sample "A"}] "bam/sorted.bam.bai"]
  [[:samtools-index {:genome "hg38", :sample "B"}] "bam/sorted.bam.bai"]
  [[:samtools-index {:genome "hg38", :sample "C"}] "bam/sorted.bam.bai"]
  [[:samtools-sort {:genome "hg38", :sample "A"}] "bam/sorted.bam"]
  [[:samtools-sort {:genome "hg38", :sample "B"}] "bam/sorted.bam"]
  [[:samtools-sort {:genome "hg38", :sample "C"}] "bam/sorted.bam"]}


(def (atom a {}))

(defn fail [] (/ 1 0))

(add-watch a :watcher
  (fn [key atom old-state new-state]
    (prn "-- Atom Changed --")
    (prn "key" key)
    (prn "atom" atom)
    (fail)
    (prn "new-state" new-state)))


(def f (future (reset! a {:a 1})))


(def ruledata
  '({:doc "Index read alignments for random access.",
  :input {0 "bam/sorted.bam"},
  :output {0 "bam/sorted.bam.bai"},
  :publish {0 "{{genome}}/{{sample}}/sorted.bam.bai"},
  :resources {:cpu 1},
  :rulename :samtools-index,
  :shell "samtools index {{input.0}}",
  :wildcards [:sample :genome]}
  {:doc "Map DNA sequences against a reference genome with BWA.",
  :externals [:genome :fastq],
  :input nil,
  :output {0 "bwa-map.bam"}
  :params {:rg "@RG\\tID:{{wildcards.sample}}\\tSM:{{wildcards.sample}}"},
  :resources {:cpu 1},
  :rulename :bwa-map,
  :shell "bwa mem -R '{{params.rg}}' -t {{threads}} {{externals.genome}} {{externals.fastq}} | samtools view -Sb - > {{output.0}}",
  :threads 8,
  :wildcards [:sample :genome]}
  {:input {0 "hiya.txt"},
  :output {0 "done.txt"},
  :resources {:cpu 1},
  :rulename :test-all,
  :shell "cat {{input.0}} > {{output.0}}",
  :wildcards [:genome :sample]}
  {:input {0 "quals.svg"},
  :output {0 "test.txt"},
  :resources {:cpu 1},
  :rulename :test-no-wildcards,
  :shell "echo 'hello!!' > {{output.0}}"}
  {:input nil,
  :output {0 "single"},
  :resources {:cpu 1},
  :rulename :hiya,
  :shell "echo 'hiya' > {{output.0}}",
  :wildcards [:genome :sample]}
  {:input {0 "all.vcf"},
  :output {:data "quals.tsv", :plot "quals.svg"},
  :resources {:cpu 1},
  :rulename :plot-quals,
  :script "/Users/endrebakkenstovner/code/everclar_mirror/examples/snakemake/plot_quals.py",
  :wildcards [:genome]}
  {:doc "Aggregate mapped reads from all samples and jointly call genomic variants.",
  :externals [:genome],
  :input {:index "bam/sorted.bam.bai", :sorted "bam/sorted.bam"},
  :output {0 "all.vcf"},
  :resources {:cpu 1},
  :rulename :bcftools-call,
  :shell "samtools mpileup -g -f {{externals.genome}} {{input.sorted}} | bcftools call -mv - > {{output.0}}",
  :wildcards [:genome]}
  {:input {0 "test.txt"},
  :output {0 "hiya.txt"},
  :resources {:cpu 1},
  :rulename :test-another-no-wildcards,
  :shell "echo 'hiya!!' > {{output.0}}"}
  {:doc "Sort the bams.",
  :input {0 "bwa-map.bam"},
  :output {0 "bam/sorted.bam"},
  :publish {0 "{{genome}}/{{sample}}/sorted.bam"},
  :resources {:cpu 1},
  :rulename :samtools-sort,
  :shell "samtools sort -T {{wildcards.sample}} -O bam {{input.0}} > {{output.0}}",
  :wildcards [:sample :genome]}))

user>

#{{:child :bcftools-call, :file "bam/sorted.bam.bai", :parent :samtools-index}
{:child :samtools-index, :file "bam/sorted.bam", :parent :samtools-sort}
{:child :plot-quals, :file "all.vcf", :parent :bcftools-call}
{:child :test-no-wildcards, :file "quals.svg", :parent :plot-quals}
{:child :test-another-no-wildcards,
:file "test.txt",
:parent :test-no-wildcards}
{:child :samtools-sort, :file "bwa-map.bam", :parent :bwa-map}
{:child :test-all, :file "hiya.txt", :parent :test-another-no-wildcards}
{:child :bcftools-call, :file "bam/sorted.bam", :parent :samtools-sort}}

(def rs [{:doc "Sort the bams.",
:input {0 "bwa-map.bam"},
:output {0 "bam/sorted.bam"},
:publish {0 "{{genome}}/{{sample}}/sorted.bam"},
:resources {:cpu 1},
:rulename :samtools-sort,
:shell "samtools sort -T {{wildcards.sample}} -O bam {{input.0}} > {{output.0}}",
:wildcards [:sample :genome]}
{:doc "Map DNA sequences against a reference genome with BWA.",
:externals [:genome :fastq],
:input nil,
:output {0 "bwa-map.bam"}
:params {:rg "@RG\\tID:{{wildcards.sample}}\\tSM:{{wildcards.sample}}"},
:resources {:cpu 1},
:rulename :bwa-map,
:shell "bwa mem -R '{{params.rg}}' -t {{threads}} {{externals.genome}} {{externals.fastq}} | samtools view -Sb - > {{output.0}}",
:threads 8,
:wildcards [:sample :genome]}])

(def basename->wildcards->jobfile {"bwa-map.bam" {:genome "hg19" :sample "A"}})

;; TODO: add an externals-map

;; Need to merge-with: out-map, in-vec, externals
;; For in-vec, out-vec, input and output we need to conj the vectors inside the jobs
(apply merge-with merge (map :in-map (:split g)))

(def g (group-by :rulename (vals @st/jobhash->job)))
(apply merge-with concat (map :in-vec (:split g)))

{:out-map
 {0
  {{:genome "hg19", :sample "B"}
   "/Users/endrebakkenstovner/everclear/snakemake-flow/results/hg19/B/chromo.txt"}},
 :in-vec
 {0
  ("/Users/endrebakkenstovner/everclear/snakemake-flow/results/0a21024cb3/genome.txt")},
 :externals {},
 :rendered-shell "echo 'hi' > {output}",
 :params nil,
 :in-map
 {0
  {{:genome "hg19"}
   "/Users/endrebakkenstovner/everclear/snakemake-flow/results/0a21024cb3/genome.txt"}},
 :output
 {0
  "/Users/endrebakkenstovner/everclear/snakemake-flow/results/hg19/B/chromo.txt"},
 :same-job {0 [:sample]},
 :jobid [:split {:genome "hg19", :sample "B"}],
 :rulename :split,
 :parent-jobs #{[:create-genome {:genome "hg19"}]},
 :wildcards {:genome "hg19", :sample "B"},
 :out-vec
 {0
  ["/Users/endrebakkenstovner/everclear/snakemake-flow/results/hg19/B/chromo.txt"]},
 :input
 {0
  "/Users/endrebakkenstovner/everclear/snakemake-flow/results/0a21024cb3/genome.txt"},
 :shell "echo 'hi' > {output}",
 :resources {:cpu 1},
 :jobhash "20962c26f8",
 :publish {0 "{{genome}}/{{sample}}/chromo.txt"}}
