(ns everclear.state.jobstate
  (:require [everclear.dag.helpers :as helpers]))

(defn ->all-files [jobs]
  (set
    (filter (complement nil?)
            (flatten
              (for [job jobs
                    :let [infiles (vals (:in-vec job))
                          externals (vals (:externals-vec job))
                          outfiles (vals (:out-vec job))]] [infiles externals outfiles])))))

(defn jobid->jobstate
  ([jobid->jobinfo]
   (jobid->jobstate
     jobid->jobinfo
     (helpers/file-existence-keyed-by-file (->all-files (vals jobid->jobinfo)))))
  ([jobid->job file->file-exists?]
   (into {}
         (for [[jobid job] jobid->job
               :let [input (flatten (vals (:in-vec job)))
                     output (flatten (vals (:out-vec job)))
                     in-files-done? (every? file->file-exists? input)
                     out-files-done? (every? file->file-exists? output)]]
           [jobid
            (cond
              ;; this function is only to be run on startup, therefore none are :in-progress
              (and in-files-done? out-files-done?) :done
              (and in-files-done? (not out-files-done?)) :ready
              :else :not-ready)]))))
