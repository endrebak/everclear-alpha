(ns everclear.jobcreation.ruleinfo
  (:require [clojure.set :as set]
            [selmer.parser :refer [render]]
            [everclear.checkers.checkers :as checkers]
            [everclear.parse.parse :as parse]
            [clojure.string :as str]))

(defn add-expand [rulename->ruleinfo wildcards]
  (let [all-expand-ks (into #{}
                            (for [{:keys [expand]} (vals rulename->ruleinfo)
                                  expand-kw (vals expand)]
                              expand-kw))
        expand-ks->wildcards (into {}
                                    (for [expand-ks all-expand-ks]
                                      [expand-ks (map #(select-keys % expand-ks) wildcards)]))]
    (apply merge
           (for [[rulename ruleinfo] rulename->ruleinfo]
             (if-let [expand-map (:expand ruleinfo)]
               {rulename (assoc ruleinfo
                           :expand
                           (into {}
                             (for [[expand-alias expand-ks] expand-map]
                               {expand-alias (get expand-ks->wildcards expand-ks)})))}
               {rulename ruleinfo})))))

(defn ->node+file [ruledata filetype]
  (let [nodekey (condp = filetype
                  :input :child
                  :output :parent
                  (throw (Exception. "Only :input and :output valid values.")))]
    (into #{}
          (for [rule ruledata
                file (vals (get rule filetype))]
            {nodekey (rule :rulename) :file file}))))

(defn join-parent-rule-and-child [ruledata]
  "Find pairs of parent and child.

  Does the same thing as ->parent+child, but before we have the rulegraph.

  Connect them by a common output/input file, respectively."
  (let [parent+file (->node+file ruledata :output)
        child+file (->node+file ruledata :input)]
    (for [{:keys [parent child]} (clojure.set/join parent+file child+file)]
      [parent child])))

(defn ->file->producer [rulename->ruleinfo]
  (apply merge
    (for [{:keys [rulename output]} (vals rulename->ruleinfo)
          [_alias filename] output]
      {filename rulename})))

(defn ->file->outalias [rulename->ruleinfo]
  (apply merge
         (for [{:keys [output]} (vals rulename->ruleinfo)
               [alias filename] output]
           {filename alias})))

(defn add-parent-info [rulename->ruleinfo]
  (let [file->producer (->file->producer rulename->ruleinfo)
        file->outalias (->file->outalias rulename->ruleinfo)]
    (apply merge
           (for [[rulename ruleinfo] rulename->ruleinfo
                 :let [input-map (:input ruleinfo)
                       output-map (:output ruleinfo)]]
             {rulename
              (assoc ruleinfo
                :parents (zipmap (keys input-map)
                                 (map (comp vector file->producer) (vals input-map)))
                :output-basenames (:output ruleinfo)
                :output-extensions (zipmap (keys output-map)
                                           (map #(second (str/split % #"\." 2)) (vals output-map)))
                                           :infile->inalias (zipmap (vals input-map) (keys input-map))
                                           :inalias->parent-outalias (zipmap
                                                                       (keys input-map)
                                                                       (map file->outalias (vals input-map))))}))))

(defn update-with-script-template [ruleinfo]
  (assoc ruleinfo :script-template (slurp (:script ruleinfo))))

(defn update-with-code-entry [rulename->ruleinfo]
  (apply merge
         (for [[rulename ruleinfo] rulename->ruleinfo]
           {rulename
            (if (:shell ruleinfo)
              ruleinfo
              (update-with-script-template ruleinfo))})))

(defn rulename->ruleinfo [ruleforms wildcards base-path log-path]
  (let [rulename->ruleinfo (-> (parse/ruleinfo-keyed-by-rulename ruleforms)
                               (add-expand wildcards)
                               add-parent-info
                               update-with-code-entry)
        rulename->ruleinfo (zipmap (keys rulename->ruleinfo)
                                   (map #(assoc % :base-path base-path :log-path log-path) (vals rulename->ruleinfo)))]
    (checkers/check-for-identical-outfiles rulename->ruleinfo)
    (checkers/check-input-fns rulename->ruleinfo)
    rulename->ruleinfo))
