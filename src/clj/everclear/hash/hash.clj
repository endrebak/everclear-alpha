(ns everclear.hash.hash
  (:require [clojure.walk :as walk]
            [nrepl.bencode :as bencode]))

(defn bencode [c]
  (-> (doto (java.io.ByteArrayOutputStream.) (bencode/write-bencode #p c)) .toString))

(defn sha256 [string]
  (let [digest (.digest (java.security.MessageDigest/getInstance "SHA-256") (.getBytes string "UTF-8"))]
    (apply str (map (partial format "%02x") digest))))

(defn nested-list [c]
  (walk/postwalk (fn [x] (if (map? x) (sort (into [] x)) x)) c))

(defn bencode-sha256 [c]
  (subs (-> c nested-list bencode sha256) 0 10))

(defn add-scripthash [rule]
  ;; In case the contents of the script change we need to update the hash of the rule
  (if-let [script (rule :script)] (assoc rule :script-hash (-> script slurp sha256)) rule))

;; IDEA: have a map of hashes in each jobid so you can tell what changed.

(defn jobid->jobhash [jobs-topo-sorted jobid->jobinfo]
  ;; TODO: should use seed-hash based on config and other stuff that affects results?
  ;; probably...
  ;; TODO: Should we hash the external files?
  ;; TODO: Should we hash the files created by rules without input. Eg. when downloading
  (let [jobid->hash (zipmap (keys jobid->jobinfo)
                            (map bencode-sha256 (vals jobid->jobinfo)))]
    (loop [jobs-topo-sorted jobs-topo-sorted
           jobid->job-and-dependencies-hash {}]
      (let [jobid (first jobs-topo-sorted)
            parent-ids (apply concat (-> jobid jobid->jobinfo :parents vals))
            dep-hashes (set (map jobid->job-and-dependencies-hash parent-ids))]
        (if-not jobid
          jobid->job-and-dependencies-hash
          (recur
           (rest jobs-topo-sorted)
           (assoc jobid->job-and-dependencies-hash jobid
                  (bencode-sha256
                   (conj dep-hashes (jobid->hash jobid))))))))))

(defn jobid->jobinfo [jobs-topo-sorted jobid->jobinfo]
  (let [jobid->jobhash (jobid->jobhash jobs-topo-sorted jobid->jobinfo)]
    (into {}
          (for [[jobid jobinfo] jobid->jobinfo
                :let [jobhash (jobid->jobhash jobid)]]
            [jobid (assoc jobinfo :hash jobhash)]))))
