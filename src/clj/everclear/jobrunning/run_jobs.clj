(ns everclear.jobrunning.run-jobs
  (:require
    [clojure.java.io :as io :refer [file]]
    [everclear.jobrunning.process :as process]
    [everclear.jobrunning.resources :as resources]
    [everclear.state.state :as app-state]
    [everclear.jobresult.jobresult :as jobresult]
    [everclear.dag.helpers :as helpers]
    [everclear.logging.log :as log]
    [clojure.set :as set]
    [everclear.state.jobstate :as jobstate])
  (:use [hashp.core]))

(defn prepare-dispatch! [jobinfo]
  (let [jobid (:jobid jobinfo)]
    (tap> {:tap-id :start-job :tap-data {:joid jobid :job jobinfo}})
    ;; In some cases (like when running samtools index on a publishpath)
    ;; only the publishpath (i.e. what is normally a symlink) is created.
    ;; To know that such jobs actually created a file, we need to remove
    ;; the symlinks first
    (helpers/delete-files! (:publish-vec jobinfo))
    (helpers/create-outfolders! jobinfo)))

(defn dispatch! [jobinfo]
  (prepare-dispatch! jobinfo)
  (let [running-process (process/->process jobinfo)]
    (future (process/read-process-and-write-logs running-process :err jobinfo))
    (future (process/read-process-and-write-logs running-process :out jobinfo))
    running-process))

(defn jobid->started-process! [jobs-to-run]
  (zipmap (map :jobid jobs-to-run)
          (mapv dispatch! jobs-to-run)))

(defn completed-procs+pending+procs [jobid->process]
  (let [{jobid+completed false jobid+pending true} (group-by #(-> % second .isAlive) jobid->process)]
    [(into {} jobid+completed) (into {} jobid+pending)]))

(defn jobs-to-run+resources-left [jobid->jobstate jobid->jobinfo job-sorter resources]
  (if-let [ready-jobs (seq (for [[k v] jobid->jobstate :when (= v :ready)] (jobid->jobinfo k)))]
        (resources/->runnable-jobs+resources-left resources (job-sorter ready-jobs))
        [nil resources]))

(defn update-resources [jobid->jobinfo jobid->completed-process resources-left]
  (let [completed-jobs (map jobid->jobinfo (keys jobid->completed-process))]
    (apply merge-with + resources-left (map :resources completed-jobs))))

;; TODO: should be able to kill jobs:
;;       just set all jobs not :done or :in-progress to :cancelled?
;;       AND kill all running processes
;; TODO: should be able to remove jobs to run:
;;       just set all jobs not :done or :in-progress to :cancelled?
(defn run-jobs! [jobid->jobinfo job-sorter jobgraph {:keys [resources summarizers]}]
  ;; should go outside run-jobs ? (tap-jobstart! @app-state/jobid->jobstate @app-state/to-run)
  (let [jobid->jobstate (jobstate/jobid->jobstate jobid->jobinfo)
        [jobs-to-run resources-left] (jobs-to-run+resources-left jobid->jobstate jobid->jobinfo job-sorter resources)]
    (loop [jobid->process (jobid->started-process! jobs-to-run)
           resources-left resources-left
           jobid->jobstate (merge jobid->jobstate (zipmap (map :jobid jobs-to-run)
                                                          (repeat :in-progress)))]
      (Thread/sleep 1000)
      (let [[jobid->completed-process jobid->pending-process] (completed-procs+pending+procs jobid->process)
            jobid->result (jobresult/process-results jobid->completed-process jobid->jobinfo resources)
            jobid->jobstate (jobresult/update-finished-jobs jobid->jobstate jobid->result jobgraph)
            resources-left (update-resources jobid->jobinfo jobid->completed-process resources-left)
            [jobs-to-run resources-left] (jobs-to-run+resources-left jobid->jobstate jobid->jobinfo job-sorter resources-left)
            jobid->pending-process (merge jobid->pending-process (jobid->started-process! jobs-to-run))]
        (log/write-logs jobid->jobinfo jobid->result summarizers)
        (when (and (seq jobid->pending-process) (not (:break @app-state/config)))
          (recur jobid->pending-process resources-left
                 (merge jobid->jobstate
                        (zipmap (map :jobid jobs-to-run)
                                (repeat :in-progress)))))))))
