(ns everclear.jobresult.jobresult
  (:require [everclear.dag.helpers :as helpers]
            [everclear.graph.graph :as graph]
            [everclear.logging.log :as log])
  (:use [hashp.core]))

(defn update-finished-jobs [jobid->jobstate jobid->jobresult jobgraph]
  ;; usage: (swap! jobstates update-jobstates jobids->jobresults hashgraph)
  (let [jobid->jobstate (merge jobid->jobstate jobid->jobresult)
        child->parents (graph/child->parents jobgraph)
        parent->children (graph/parent->children jobgraph)
        children-of-done-jobs (set (for [[jobid jobresult] jobid->jobresult
                                         :when (= jobresult :done)
                                         child (get parent->children jobid)]
                                     child))
        jobid->ready (into {}
                           (for [child children-of-done-jobs
                                 :let [parents (child->parents child)
                                       parent-jobstates (map jobid->jobstate parents)
                                       all-parents-done? (every? (partial = :done) parent-jobstates)]
                                 :when all-parents-done?]
                             [child :ready]))]
    (merge jobid->jobstate jobid->ready jobid->jobresult)))

(defn handle-success [jobid existence->outfile+publishfile]
  (tap> {:tap-id :job-finished :tap-data {:jobid jobid}})
  (helpers/symlink-missing-files! existence->outfile+publishfile)
  [jobid :done])

(defn handle-failure [jobid existence->outputfile+publishfile exit-code]
  (let [missing-files (map first (vals (select-keys existence->outputfile+publishfile [:both-missing :publish-missing])))
        existing-publishfiles (map second (:publish-missing existence->outputfile+publishfile))
        existing-outputfiles (map first (:publish-missing existence->outputfile+publishfile))]
    (tap> {:tap-id :job-failed :tap-data {:exit-code exit-code :missing-files missing-files}})
    (helpers/delete-files! (concat existing-outputfiles existing-publishfiles))
    [jobid :failed]))

(defn process-result [jobid exit-code jobinfo]
  (let [file->exists? (helpers/file-existence-keyed-by-file (helpers/files-and-symlinks-to-be-created jobinfo))
        existence->outfile+publishfile (helpers/->existence->outfile+publishfile jobinfo file->exists?)
        file-creation-successful? (every? (partial not= :both-missing) (vals existence->outfile+publishfile))]
    (if (and (= 0 exit-code) file-creation-successful?)
      (handle-success jobid existence->outfile+publishfile)
      (handle-failure jobid existence->outfile+publishfile exit-code))))

(defn process-results [jobid->completed jobid->jobinfo config]
  (let [jobid->jobresult
        (into {}
              (for [[jobid process] jobid->completed
                    :let [jobinfo (get jobid->jobinfo jobid)
                          jobresult (process-result jobid (.exitValue process) jobinfo)]]
                jobresult))]
    jobid->jobresult))
