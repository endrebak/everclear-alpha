#!/usr/bin/env bb

(require '[babashka.curl :as curl])

(require '[clojure.tools.cli :refer [parse-opts]])
(require '[clojure.string :as str])
(require '[clojure.java.shell :as sh])
;; (require '[io.aviso.ansi :as ansi])

;; optional arguments:
;;   -h, --help            show this help message and exit

;; EXECUTION:
;;   target                Targets to build. May be rules or files. (default:
;;                         None)
;;   --dry-run, --dryrun, -n
;; --
;; GROUPING:
;;   --groups GROUPS [GROUPS ...]
;;                         Assign rules to groups (this overwrites any group
;;                         definitions from the workflow). (default: None)
;; --
;; REPORTS:
;;   --report [FILE]       Create an HTML report with results and statistics.
;;                         This can be either a .html file or a .zip file. In the
;;                         former case, all results are embedded into the .html
;; --
;; NOTEBOOKS:
;;   --draft-notebook TARGET
;;                         Draft a skeleton notebook for the rule used to
;;                         generate the given target file. This notebook can then
;; --
;; UTILITIES:
;;   --lint [{text,json}]  Perform linting on the given workflow. This will print
;;                         snakemake specific suggestions to improve code quality
;;                         (work in progress, more lints to be added in the
;; --
;; OUTPUT:
;;   --reason, -r          Print the reason for each executed rule. (default:
;;                         False)
;;   --gui [PORT]          Serve an HTML based user interface to the given

(def cli-options
  ;; An option with a required argument
  [["-P" "--port PORT" "Port number of everclear instance."
    :default 3000
    :parse-fn #(Integer/parseInt %)
    :validate [#(< 0 % 0x10000) "Must be a number between 0 and 65536"]]
   ["-S" "--start" "Start everclear."]
   ["-c" "--config CONFIG-FILE" "Which config file to use."]
   ["-C" "--switch-config" "Switch configuration file. Requires the -c/--config flag."]
   ["-U" "--url URL" "Url to the running everclear instance." :default "localhost"]
   ["-r" "--run" "Run selected DAG."]
   ["-f" "--from RULE" "Choose jobs from RULE in the DAG."]
   ["-u" "--until RULE" "Choose jobs until RULE in the DAG."]
   ["-w" "--wildcards WILDCARDS" "Choose jobs with these WILDCARDS."]
   ["-d" "--dry-run" "Do not execute command, but display what would be done."]
   ["-f" "--force" "Run selected rules even though their output exists."]
   ["-s" "--show-config" "Show the configuration options used."]
   ;; ["-cs" "--show-config" "Show the configuration options used."]
   ["-h" "--help"]])

(def args (parse-opts *command-line-args* cli-options))

(def flags (:options args))

(def COMMANDS #{:show-config :switch-config :run :dry-run :start})

(defn parse-flags [flags]
  (let [commands (remove nil? (map COMMANDS (keys flags)))]
    (case (count commands)
       0 (throw (Exception. (str "No commands found. Choose one of: " (str COMMANDS))))
       1 (let [command (first commands)
               options (dissoc flags command :url :port)]
               ;; _ (println command)
               ;; _ (println options)]
           {:options options :command command})
       (throw (Exception. (str "Multiple commands found: " (str/join " " commands)))))))

(when (:help flags)
  (println (:summary args))
  (System/exit 0))

(def url (:url flags))
(def port (:port flags))
(def url (str url ":" port "/json"))

(def parsed-flags (parse-flags flags))

(when (= :start (:command parsed-flags))
  (let [cmd ["lein" "run" "-c" (get-in parsed-flags [:options :config])]
        result (apply clojure.java.shell/sh cmd)]
    (println result)
    (System/exit 0)))

(defn send-command [url parsed-flags]
  (curl/post url {:body (str parsed-flags)
                  :headers {"Content-Type" "application/edn"
                            "Accept" "application/edn"}}))

(def result (send-command url parsed-flags))

(println (:body result))
