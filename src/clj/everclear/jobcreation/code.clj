(ns everclear.jobcreation.code)

(defn update-with-script [job]
  (let [script-template (slurp (:script job))
        script (render script-template job)]
    (assoc job :rendered-script script)))

(defn update-with-code-entry [job]
  (if (:shell job) ;; should be cond
    (assoc job :rendered-shell (render (:shell job) job))
    (update-with-script job)))
