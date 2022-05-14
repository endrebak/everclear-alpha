(ns everclear.checkers.checkers
  (:require [clojure.string :as str]))

(defn add-output-vars-used [rulename->ruleinfo]
  ;; this can probably be checked with selmer?
  (into {}
        (for [[rulename ruleinfo] rulename->ruleinfo
              :let [code (cond
                           (:shell ruleinfo) (:shell ruleinfo)
                           (:script ruleinfo) (slurp (:script ruleinfo))
                           :else (throw (Exception. (str rulename "Missing script or shell."))))
                    output-vars-used (map first (re-seq #"\{\{(output|publish)\.[a-zA-Z.]*\}\}" code))]]
          [rulename (assoc ruleinfo :output-vars-used output-vars-used)])))

(defn check-for-identical-outfiles [rulename->ruleinfo]
  (let [outfile->rulenames (apply
                             merge-with
                             concat
                             (for [[rulename ruleinfo] rulename->ruleinfo
                                   f (vals (:output ruleinfo))]
                               {f [[rulename]]}))
        duplicate-outfiles (for [[_k v] outfile->rulenames :when (> (count v) 1)] v)]
    (when (seq duplicate-outfiles)
      (throw (Exception.
               (str "These rules have the same outfile: "
                    (str/join ", " (first duplicate-outfiles))))))))

(defn ensure-publishpath-has-enough-wildcards [ruleinfo]
  ;; 1. for each publishentry
  ;; 2. check that the number of wildcards + expanded wildcards
  ;;    is equal to the number of {{}} in the publishpath
  )

(defn check-that-no-files-have-multiple-aliases [])

(defn check-input-fns [rulename->ruleinfo]
  (for [[rulename ruleinfo] rulename->ruleinfo
        [inalias input-fn] (:input-fns ruleinfo)]
    (assert
      (fn? input-fn)
      (throw (Exception. (str "input-fn for rule " rulename ", inputalias " inalias " is not a valid function: " input-fn))))))

