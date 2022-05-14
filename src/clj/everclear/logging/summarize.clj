(ns everclear.logging.summarize
  (:require [clojure.string :as str]
            [clojure.java.shell :as shell]
            [clojure.java.io :as io]))

(defn summarize-file [infile summarize-script outpath]
  (let [result (shell/sh "bb" summarize-script "-f" infile)]
    (io/make-parents outpath)
    (spit outpath result)
    result))

(defn write-summaries! [jobinfo summarizers]
  (apply merge-with merge
         (doall
           (for [[output-alias ext] (:output-extensions jobinfo)
                 :let [summarize-script (get summarizers ext)
                       current-summary (get-in jobinfo [:summary-map output-alias])]
                 :when (seq summarize-script) ;; TODO: change to default (binary/non-binary)
                 [xp-wcs outfile] (get-in jobinfo [:output-map output-alias])
                 :let [summary-path (get current-summary xp-wcs)]]
             {output-alias {xp-wcs (summarize-file outfile summarize-script summary-path)}}))))
