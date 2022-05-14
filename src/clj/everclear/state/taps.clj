(ns everclear.state.taps
  (:require
   [everclear.state.state :as app-state]
   [clojure.string :as str]
   [clojure.term.colors :refer [yellow red green]]))

(defn ->cmdline-str [[rulename job-wildcards] std-out-or-err line]
  (let [std-out-or-err-color (if (= std-out-or-err :err) red green)]
    (str (std-out-or-err-color (name rulename)) " ("
         (str/join " "
                   (for [[k v] job-wildcards]
                     (str/join " " [(green (name k)) v])))
         ")\n  " line)))

;; [Fri Jan 14 16:02:52 2022]
;; rule bwa_map:
;;     input: data/genome.fa, data/samples/A.fastq
;;     output: mapped_reads/A.bam
;;     jobid: 0
;;     resources: tmpdir=/var/folders/tz/9155ldys7xg36pvt1dxszwk80000gn/T

;; (assert false "Include whether job is stale or fresh in output.")
;; (assert false "Use emojis in output (stderr: ⚠️) etc)??")
;; (assert false "Have blank lines before after starting/finishing jobs")

(defn ->start-job-str [job]
  (let [now (.format
             (java.text.SimpleDateFormat. "EEE MMM dd HH:mm:ss yyyy")
             (new java.util.Date))]
    (str/join "\n"
              [(str "--- Starting job " (name (job :rulename)) " " (job :wildcards))
               (str "--- " now)
               (or
                (job :rendered-shell)
                (job :rendered-script))])))

(defn tap [{:keys [tap-id tap-data]}]
  (let [jobid (:jobid tap-data)]
    (try
      (condp = tap-id
        :process-output
        (let [cmdline (->cmdline-str jobid (:stream tap-data) (:line tap-data))]
          (println cmdline))
        :start-job
        (let [job (:job tap-data)]
          (println (->start-job-str job)))
        :job-finished
        (println (str "----- " jobid " finished"))
        :job-failed
        (println (str "----- Job FAILED: " jobid " with exit code " (tap-data :exit-code) " and missing files: " (str/join " " (vals (tap-data :missing-jobfiles)))))
        :no-ready-jobs
        (println (str "No ready jobs! Only jobs that are " (str/join ", " (seq (:jobstates tap-data)))))
        :no-to-run-jobs
        (println (str "No jobs scheduled to run!"))

        :no-to-run-and-ready-jobs
        (let [to-run (:to-run tap-data)
              jobhash->jobstate (:jobhash->jobstate tap-data)]
          (println (str "The jobs that are to-run have state: " (for [jobhash to-run] (jobhash->jobstate jobhash)))))
        ;; :jobs-in-progress
        ;; (println "Some ")
        :error
        (let [{:keys [cause trace]} tap-data]
          (println (str "ERROR: " cause trace)
                   (reset! app-state/error {:cause cause :trace (map str trace)})))
        (print (str "TAPTAPTAP:" tap-id tap-data)))
      (catch Throwable e
        (let [{:keys [_trace cause]} (Throwable->map e)]
          (println (str "Error in tap function: " cause)))))))
;;       (def jn (first j))
;;       (def jw (second j))
;;       (str (red (name jn)) "for " (str/join " " (for [m jw [k v] m] (str/join " "[(green (name k)) v])))))

;; (str (name jn) ": " (str/join " " (for [m jw [k v] m] (str/join " " [(name k) v]))))

;;; How do I connect to my running process from a second REPL?
;;; My workflow is:
;;; start cider
;;; start my luminus app
;;; see output from the app in my cider buffer.
;;;
;;; But when I try to connect to the nrepl port I am told the webserver uses with lein repl :connect 7000
;;; 2021-06-11 17:58:26,293 [nREPL-session-aa73c1f3-2614-4d77-b9bf-4a6e303c67f2] INFO  everclear.nrepl - starting nREPL server on port 7000
;;; I do not see any of the output from my app in that buffer.
;;;
;;; What gives?
