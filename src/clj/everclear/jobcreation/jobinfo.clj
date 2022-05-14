(ns everclear.jobcreation.jobinfo
  (:require
    [clojure.set :as set]
    [everclear.dag.helpers :as helpers]
    [everclear.jobcreation.wildcards :as wildcards]))

(defn ->infile->parents [f ruleinfo wildcards rulename->ruleinfo]
  (apply merge
         (for [[infile-alias [parent-rulename]] (:parents ruleinfo)
               :let [p-wc-ks (get-in rulename->ruleinfo [parent-rulename :wildcards])
                     p-wc-t (f wildcards p-wc-ks)]]
           {(get-in ruleinfo [:input infile-alias]) (vec (for [p-wcs (helpers/sort-wc-rows-by-ks p-wc-t)] [parent-rulename p-wcs]))})))

(defn apply-input-fn [input-fn [rulename wildcards] inalias]
  (if (nil? input-fn)
    [wildcards]
    (let [-fn (eval input-fn)
          wc-rows (-fn wildcards)
          wc-row-keys (set (map keys wc-rows))]
      (assert (= 1 (count wc-row-keys))
              (str "Different wildcards produced by input-fn " inalias " of rule " rulename "."
                   "\nThe wildcard rows look like: " wc-rows))
      wc-rows)))

(defn ->inalias->req-wcs [jobid ruleinfo]
  (apply merge
         (for [[inalias _infile] (:input ruleinfo)
               :let [input-fn (get-in ruleinfo [:input-fns inalias])]]
           {inalias (apply-input-fn input-fn jobid inalias)})))

(defn ->parent-jobid->children-jobids [jobid->protojob]
  (apply merge-with concat
         (for [[jobid job] jobid->protojob
               [inalias wildcard-rows] (->inalias->req-wcs jobid job)
               parent-wcs wildcard-rows
               :let [parent-rule (get-in job [:parents inalias])]]
           {[parent-rule parent-wcs] [jobid]})))

(defn outjobfiles [ruleinfo wcs]
  (apply merge-with concat
         (for [[outalias outfile] (:output ruleinfo)
               :let [expansions (get-in ruleinfo [:expand outalias] [{}])]
               xp-wc-row (map #(merge wcs %) expansions)]
           {outalias [outfile xp-wc-row]})))

(defn injobfiles [ruleinfo inalias->req-wcs]
  (apply merge-with concat
         (for [[inalias infile] (:input ruleinfo)
               req-wcs (get inalias->req-wcs inalias)]
           {inalias [[infile req-wcs]]})))

(defn jobparents [ruleinfo inalias->req-wcs]
  (apply merge-with concat
         (for [[inalias parent-rule] (:parents ruleinfo)
               req-wcs (get inalias->req-wcs inalias)]
           {inalias [[parent-rule req-wcs]]})))

(defn ->jobid->protojob [{:keys [wildcards rulename] :as ruleinfo} wc-rows]
  (apply
    merge
    (for [wcs wc-rows
          :let [jobid [rulename wcs]
                inalias->req-wcs (->inalias->req-wcs
                                   jobid
                                   ruleinfo)]]
      {[rulename wcs]
       (assoc ruleinfo
         :wildcards wcs
         :jobid [rulename wildcards]
         :output (outjobfiles ruleinfo wcs)
         :input (injobfiles ruleinfo inalias->req-wcs)
         :parents (jobparents ruleinfo inalias->req-wcs))})))

(defn ->wc-rows [rulename jobids wildcards wc-ks->wcs]
  (let [wc-rows (if (seq jobids)
                  (map #(select-keys % wildcards) (map second jobids))
                  (get wc-ks->wcs wildcards))
        missing-wildcards (set/difference (set wildcards) (-> wc-rows first keys set))
        _ (assert
            (not (seq missing-wildcards))
            (str "Rule " rulename " got too few wildcards requested. Missing: " missing-wildcards
                 ".\nSome input function in some child rule did not request enough wildcards."))]
    wc-rows))

(defn ->rulename->jobids [jobid->jobinfo]
  (let [inalias->parents (apply concat (map :parents (vals jobid->jobinfo)))
        rulename->jobids (group-by first (vals inalias->parents))]
    rulename->jobids))

(defn jobid->protojob [rules-topo-sorted rulename->ruleinfo wildcards-table]
  (let [wc-ks->wcs (wildcards/->wc-ks->wcs rulename->ruleinfo wildcards-table)]
    (loop [rules-reverse-topo-sorted (reverse rules-topo-sorted)
           rulename->jobids {}
           jobid->jobinfo {}]
      (let [rulename (first rules-reverse-topo-sorted)
            wildcards (get-in rulename->ruleinfo [rulename :wildcards])
            wc-rows (->wc-rows rulename (get rulename->jobids rulename) wildcards wc-ks->wcs)
            -jobid->jobinfo (->jobid->protojob (get rulename->ruleinfo rulename) wc-rows)
            -rulename->jobids (->rulename->jobids -jobid->jobinfo)
            jobid->jobinfo (merge jobid->jobinfo -jobid->jobinfo)
            remaining-rules (rest rules-reverse-topo-sorted)]
        (if (not (seq remaining-rules))
          jobid->jobinfo
          (recur
            remaining-rules
            (merge rulename->jobids -rulename->jobids)
            jobid->jobinfo))))))