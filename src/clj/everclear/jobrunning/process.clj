(ns everclear.jobrunning.process
  (:require [clojure.java.io :as io]
            [selmer.parser :refer [render]]
            [clojure.string :as str]))

;; Writing the script to a file here instead of when the jobs are created to
;; avoid a lot of I/O upfront
(defn write-script-and-return-exec [{:keys [base-path hash script rendered-script]}]
  (let [script-output-path (->> (io/file script)
                                .getName
                                (io/file base-path hash)
                                .toString)]
    (spit script-output-path rendered-script)
    ["python" script-output-path]))

(defn code [job]
  (if (:shell job)
    ["bash" "-c" (:rendered-shell job)]
    (write-script-and-return-exec job)))

(defn ->process [job]
  (-> (code job)
      (ProcessBuilder.)
      (.start)))

(defn read-process [proc stream jobid]
  (try
    (with-open [rdr (io/reader (cond
                                 (= :out stream) (.getInputStream proc)
                                 (= :err stream) (.getErrorStream proc)
                                 :else (throw (Exception. (str "Only :stderr or :stdout valid streams. Was: " stream)))))]
      (binding [*in* rdr]
        (loop [lines []]
          (let [line (read-line)]
            (if (nil? line)
              lines
              (do
                ;; (tap> {:tap-id :process-output :tap-data {:jobid jobid :line line :stream stream}})
                (recur (conj lines line))))))))
    (catch Exception e
      (println "ooops\n" (pr-str e)))))

;; Do you also want to create the logs here?
;; Yes, want to do it in a future
(defn read-process-and-write-logs [proc stream {:keys [jobid hash log-path]}]
  (let [outfile (io/file log-path hash (str "std" (name stream) ".log"))
        lines (read-process proc stream jobid)]
    (io/make-parents outfile)
    (spit outfile (str/join "\n" lines))))