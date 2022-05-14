(ns everclear.jobcreation.paths
  (:require
    [clojure.java.io :as io]
    [everclear.dag.helpers :refer [make-absolute]]
    [everclear.jobcreation.wildcards :as wildcards]
    [selmer.parser :refer [render]]
    [cuerdas.core :as str]
    [everclear.dag.helpers :as helpers])
  (:use hashp.core))

(defn hashpath [hash xp-wcs outfile rulename]
  (let [xp-path (if (seq xp-wcs)
                  (apply io/file (map xp-wcs (sort (keys xp-wcs))))
                  "")]
    (.toString (io/file hash (name rulename) xp-path outfile))))

(defn wildcards->outpath [outjobfiles {:keys [rulename hash]}]
   (apply merge
          (for [[outfile expanded-wildcards] outjobfiles]
            {expanded-wildcards (hashpath hash expanded-wildcards outfile rulename)})))

(defn output-map [jobinfo]
  (apply merge
         (for [[outalias outjobfiles] (:output jobinfo)]
           {outalias
            (wildcards->outpath outjobfiles jobinfo)})))

(defn render-publishpath [base-path publishtemplate expanded-wildcards]
  (->> (render publishtemplate expanded-wildcards)
       (make-absolute base-path)))

(defn wildcards->publishpath [base-path publish-template outjobfiles]
  (apply merge
         (for [[expanded-wildcard-row _outfile] outjobfiles]
           {expanded-wildcard-row (render-publishpath base-path publish-template expanded-wildcard-row)})))

(defn publish-map [jobinfo]
  (apply merge
         (for [[outalias outjobfiles] (:output-map jobinfo)]
           {outalias
            (wildcards->publishpath (:base-path jobinfo) (get-in jobinfo [:publish outalias]) outjobfiles)})))

(defn wildcards->inpath [jobid->protojob inalias protojob parent-jobid]
  (let [parent-outalias (get-in protojob [:inalias->parent-outalias inalias])
        output-map (get-in jobid->protojob [parent-jobid :output-map parent-outalias])
        publish-map (get-in jobid->protojob [parent-jobid :publish-map parent-outalias])]
    (if (seq publish-map)
      publish-map
      output-map)))

(defn input-map [jobid->protojob protojob infile->inalias]
  (apply merge-with merge
         (for [[infile parent-jobids] (:parents protojob)
               :let [inalias (infile->inalias infile)]
               parent-jobid parent-jobids]
           {inalias (wildcards->inpath jobid->protojob inalias protojob parent-jobid)})))

(defn externals-map [external-aliases externals wildcards]
  (apply merge-with merge
         (for [external-alias external-aliases
               [wildcards-key external-file] (get externals external-alias)
               :when (wildcards/submap? wildcards wildcards-key)]
           {external-alias {wildcards-key external-file}})))

(defn wildcards-sorted-vec-from-map [alias->filemap]
  (into {}
        (for [[alias m] alias->filemap]
          (if (and (= 1 (count m)) (-> m first key empty?)) ;; the wildcards-map is empty
            [alias [(-> m first val)]]
            [alias (vec (helpers/sort-vals-by-wcmap-keys m))]))))

(defn add-vec-and-str [jobinfo alias->vs key-name]
  (merge jobinfo
         (apply merge-with merge
                (for [[alias v] alias->vs
                      :let [vec-name (-> key-name name (str/replace "put" "") (str "-vec") keyword)
                            s (str/join " " v)]]
                  {key-name {alias s} vec-name {alias v}}))))

(defn add-basepath [basepath m ext]
  (apply merge-with merge
         (for [[alias jobfiles] m
               [xp-wcs jobfile] jobfiles]
           {alias
            {xp-wcs (.toString (io/file basepath (str jobfile ext)))}})))

(defn add-output [jobinfo]
  (let [proto-output-map (output-map jobinfo)
        {:keys [base-path log-path]} jobinfo
        output-map (add-basepath base-path proto-output-map "")
        summary-map (add-basepath log-path proto-output-map ".log")]
    (-> jobinfo
        (assoc :output-map output-map)
        (assoc :protopath-map proto-output-map)
        (assoc :summary-map summary-map)
        (add-vec-and-str (wildcards-sorted-vec-from-map output-map) :output))))

(defn add-publish [jobinfo]
  (let [publish-map (publish-map jobinfo)]
    (-> jobinfo
        (assoc :publish-map publish-map)
        (add-vec-and-str (wildcards-sorted-vec-from-map publish-map) :publish))))

(defn add-externals [jobinfo externals]
  (let [externals-map (externals-map (:externals jobinfo) externals (:wildcards jobinfo))]
    (-> jobinfo
        (assoc :externals-map externals-map)
        (add-vec-and-str (wildcards-sorted-vec-from-map externals-map) :externals))))

(defn add-input [jobid->jobinfo]
  (zipmap (keys jobid->jobinfo)
          (for [jobinfo (vals jobid->jobinfo)]
            (if (seq (:input jobinfo))
              (let [infile->inalias (zipmap (vals (:input jobinfo)) (keys (:input jobinfo)))
                    input-map (input-map jobid->jobinfo jobinfo infile->inalias)]
                (-> jobinfo
                    (assoc :input-map input-map)
                    (add-vec-and-str (wildcards-sorted-vec-from-map input-map) :input)))
              jobinfo))))

(defn add-paths [jobid->jobinfo externals]
  (let [jobid->jobinfo
        (zipmap (keys jobid->jobinfo)
                (for [jobinfo (vals jobid->jobinfo)]
                  (cond-> jobinfo
                          (seq (:output jobinfo)) add-output
                          (seq (:publish jobinfo)) add-publish
                          (seq (:externals jobinfo)) (add-externals externals))))]
    (add-input jobid->jobinfo)))