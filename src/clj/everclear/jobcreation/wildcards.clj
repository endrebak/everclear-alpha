(ns everclear.jobcreation.wildcards)

(defn ->wc-ks->wcs [rulename->ruleinfo wildcards-table]
  (let [all-wc-ks (map :wildcards (vals rulename->ruleinfo))]
    (apply merge
           (for [wc-ks all-wc-ks
                 :let [wc-t (or
                              (-> (map #(select-keys % wc-ks) wildcards-table) set seq)
                              [{}])]]
             {wc-ks wc-t}))))

(defn submap?
  "Checks whether m contains all entries in sub."
  [^java.util.Map m ^java.util.Map sub]
  (.containsAll (.entrySet m) (.entrySet sub)))

(defn ->c-wcs+p-wc-ks->p-wc-t [wildcards-table c-wcs p-wc-ks]
  (for [wc-row wildcards-table
        :when (submap? wc-row c-wcs)]
    (select-keys wc-row p-wc-ks)))