(ns everclear.dag.core
  (:require
    [everclear.parse.parse :as parse]
    [everclear.parse.create-user-ns :as create-user-ns]
    [everclear.jobcreation.ruleinfo :as ruleinfo]
    [everclear.jobcreation.jobinfo :as jobinfo]
    [everclear.graph.graph :as graph]
    [everclear.jobcreation.paths :as paths]
    [everclear.state.filewatch :as filewatch]
    [everclear.state.state :as app-state]
    [everclear.state.jobstate :as jobstate]
    [everclear.routes.ws :as ws]
    [everclear.dag.helpers :refer [file-to-map]]
    [everclear.hash.hash :as hash]
    [selmer.parser :refer [render]]))

(defn add-params [jobinfo]
  (merge jobinfo
         (if-let [params (seq (:params jobinfo))]
           {:params (zipmap (keys params)
                            (map #(render % jobinfo) (vals params)))})))

(defn add-code [jobinfo]
  (merge jobinfo
         (cond
           (seq (:shell jobinfo)) {:rendered-shell (render (jobinfo :shell) jobinfo)}
           (seq (:script jobinfo)) {:rendered-script (render (jobinfo :script-template) jobinfo)}
           :else (throw (Exception. "Jobinfo must contain :script or :shell.")))))

(defn parent-hash-map [jobinfo jobid->jobinfo]
  (let [infile->inalias (:infile->inalias jobinfo)]
    (apply merge-with #(into [] (apply concat %&))
           (for [[infile parent-jobids] (:parents jobinfo)
                 :let [inalias (infile->inalias infile)]
                 parent-jobid parent-jobids]
             {inalias [(get-in jobid->jobinfo [parent-jobid :hash])]}))))

(defn create-graphs&jobs [rule-file wildcards-file externals-file base-path log-path]
  (let [wildcards-table (file-to-map wildcards-file)
        forms+string (parse/read+string-forms-in-file rule-file)
        ruleforms #p (create-user-ns/eval-user-ns rule-file)
        var->val (create-user-ns/var->val forms+string)
        _ #p (create-user-ns/update-ruleforms ruleforms var->val)]
        ;; rulename->ruleinfo (ruleinfo/rulename->ruleinfo ruleforms wildcards-table base-path log-path)
        ;; rulegraph (graph/graph rulename->ruleinfo)
        ;; jobid->jobinfo (jobinfo/jobid->protojob (graph/topo-sort rulegraph) rulename->ruleinfo wildcards-table)
        ;; externals (file-to-map externals-file)
        ;; jobgraph (graph/graph jobid->jobinfo)]
        ;; jobid->jobinfo (-> (hash/jobid->jobinfo (graph/topo-sort jobgraph) jobid->jobinfo)
        ;;                    (paths/add-paths externals))]
        ;; jobid->jobinfo (zipmap (keys jobid->jobinfo)
        ;;                        (for [jobinfo (vals jobid->jobinfo)]
        ;;                          (assoc jobinfo
        ;;                            :parent-hashes (parent-hash-map jobinfo jobid->jobinfo))))
        ;; jobid->jobinfo (zipmap (keys jobid->jobinfo)
        ;;                        (map #(-> % add-params add-code) (vals jobid->jobinfo)))]
    {;; :rulegraph rulegraph
     ;; :jobgraph jobgraph
     ;; :jobid->jobinfo jobid->jobinfo
     ;; :rulename->ruleinfo rulename->ruleinfo

     }))

(defn update-graphs-and-jobs! [{:keys [rulegraph rulegraph-map rulename->ruleinfo jobid->jobinfo jobgraph] :as graphs&jobs} summarizers]
  (reset! app-state/error nil)
  (reset! app-state/rulename->ruleinfo rulename->ruleinfo)
  (reset! app-state/jobid->jobinfo jobid->jobinfo)
  (reset! app-state/rulegraph rulegraph)
  (reset! app-state/jobgraph jobgraph)
  (reset! app-state/rulegraph-map rulegraph-map)
  (ws/send-data-over-websockets :jobinfo/update jobid->jobinfo)
  (ws/send-data-over-websockets :rulegraph/update rulegraph-map)
  graphs&jobs)

(defn create-and-update! [rule-file wildcards-file externals-file config-file]
  (try
    (let [config (assoc (file-to-map config-file) :config-file config-file)
          {:keys [summarizers base-path log-path]} config]
      (reset! app-state/config config)
      (reset! app-state/resources (reset! app-state/resources (config :resources)))
      (->
       (create-graphs&jobs rule-file wildcards-file externals-file base-path log-path)
       (update-graphs-and-jobs! summarizers)))
    (catch Throwable e
      (let [{:keys [trace cause]} (Throwable->map e)]
        (tap> {:tap-id :error :tap-data {:cause cause :trace (map str trace)}})))))

(defn start-everclear [opts]
    (let [config-file (get-in opts [:options :config-file])
          config (file-to-map config-file)
          {:keys [wildcards-file rule-file externals-file]} config]

      (filewatch/initiate-watch!
       :jobinfo-watch
       [wildcards-file rule-file externals-file config-file]
       (fn [_ _] (create-and-update! rule-file wildcards-file externals-file config-file)))

      (create-and-update! rule-file wildcards-file externals-file config-file)))

;; The config-file is not part of the config iteself so when it is changed
;; (only possibly from within the program) we need to update the filewatcher so that it
;; points to the new config-file
(add-watch
 app-state/config
 :config-watcher
 (fn [_ _ old-config new-config]
   (if (nil? new-config)
     (throw (Exception. "Trying to set config to nil.")))
   (when (not= (old-config :config-file) (new-config :config-file))
     (try
       (let [config-file (new-config :config-file)
             {:keys [rule-file wildcards-file externals-file]} (file-to-map config-file)]
         (filewatch/update-watch!
          :jobinfo-watch
          [wildcards-file rule-file externals-file config-file]
          (fn [_ _] (create-and-update! rule-file wildcards-file externals-file config-file)))
         (create-and-update! rule-file wildcards-file externals-file config-file))
     (catch Throwable e
       (let [{:keys [trace cause]} (Throwable->map e)]
         (tap> {:tap-id :error :tap-data {:cause cause :trace (map str trace)}})))))))
