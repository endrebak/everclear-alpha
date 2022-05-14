(ns everclear.jobrunning.resources)

(defn ->runnable-jobs+resources-left [resources-available prioritized-ready-jobs]
  ;; TODO: this function does not update the resources
  ;; TODO: some jobs might require more resources than the total available
  (loop [jobs prioritized-ready-jobs
         resources resources-available
         runnable-jobs []]
    (let [[job jobs-left] [(first jobs) (rest jobs)]
          resources-left (merge-with - resources (:resources job))
          job-runnable? (not-any? neg? (vals resources-left))
          runnable-jobs (if job-runnable? (conj runnable-jobs job) runnable-jobs)]
      (if-not (and (seq jobs-left) job-runnable?)
        [runnable-jobs resources]
        (recur jobs-left resources-left runnable-jobs)))))
