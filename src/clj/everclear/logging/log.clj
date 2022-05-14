(ns everclear.logging.log
  (:require [clojure.java.io :as io]
            [everclear.logging.summarize :as summarize]))

(defn write-jobinfo [{:keys [log-path hash] :as jobinfo}]
  (let [jobinfo-file (io/file log-path hash "jobinfo.edn")]
    (io/make-parents jobinfo-file)
    (spit jobinfo-file jobinfo)))

(defn cp-graphs [{:keys [protopath-map output-map output-extensions]}]
  (doseq [[output-alias ext] output-extensions
          :when (#{"svg" "png" "jpg" "jpeg" "gif"} ext)
          [wildcards outfile] (get output-map output-alias)
          :let [protopath (get-in protopath-map [output-alias wildcards])
                resource-path (io/file "resources" "public" "img" protopath)]]
    (everclear.dag.helpers/cp! outfile resource-path)))

(defn read-logs [])

(defn write-logs [jobid->jobinfo jobid->result summarizers]
  (let [done-jobs (for [[k v] jobid->result :when (= v :done)]
                       (jobid->jobinfo k))]
    (doseq [{:keys [hash log-path] :as jobinfo} done-jobs]
      (do
        (summarize/write-summaries! jobinfo summarizers)
        (cp-graphs jobinfo)
        (spit (io/file log-path hash "jobinfo.edn") (with-out-str (clojure.pprint/pprint (into (sorted-map) jobinfo))))))))

