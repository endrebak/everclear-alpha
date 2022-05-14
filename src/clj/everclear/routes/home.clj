(ns everclear.routes.home
  (:require
    [everclear.layout :as layout]
    [everclear.graph.graph :as graph]
    [everclear.state.state :as app-state]
    [clojure.java.io :as io]
    [everclear.middleware :as middleware]
    [everclear.routes.ws :as ws]
    [muuntaja.middleware :as muuntaja]
    [ring.util.response]
    [ring.util.request :refer [body-string]]
    [clojure.pprint]
    [clojure.data.json :as json]
    [ring.util.http-response :as response]
    [clojure.string :as str]
    [everclear.dag.helpers :as helpers]))

;; ;; TODO: make nicer function, refactor
;; ;; should fetch all relevant flags and options and select a command to run
;; ;;
;; (assert false "Handle options here!")

(defn handle-options [flags]
  (let [all-tasks #{:run-jobs}])
  ;; check wildcards
  ;; check job
  ;; check that no more than one command is given
  )

(defn read-wildcards [wc-str]
  (let [m (clojure.edn/read-string wc-str)]
    (zipmap (keys m) (map str (vals m)))))

(defn rulename->wildcard-rows [jobs]
  (for [[rulename groups] (group-by first jobs)
        :let [wildcards (map second groups)]]
    [rulename (helpers/sort-wc-rows-by-ks wildcards)]))

(defn dry-run [jobs]
  (let [srtd (helpers/topo-sort-jobids-by-rulename
               (rulename->wildcard-rows jobs)
               @everclear.state.state/rulegraph)]
    (clojure.string/join "\n\n"
                         (for [[rn wcs] srtd
                               :let [wc-str (str/join " " wcs)]]
                           (str rn "\n" wc-str)))))

(defn jobs [{:keys [until wildcards]}]
        (let [jg @everclear.state.state/jobgraph
              until (keyword until)
              wildcards (read-wildcards wildcards)]
          (graph/all-ancestors jg (graph/jobs jg until wildcards))))

(defn json-handler [request]
  (let [options (get-in request [:body-params :options])
        command (get-in request [:body-params :command])
        jobgraph (jobs options)]
    (case command
      :run (let [jobid->jobinfo (select-keys @app-state/jobid->jobinfo jobgraph)]
             (response/ok (everclear.jobrunning.run-jobs/run-jobs!
             jobid->jobinfo identity @app-state/jobgraph @app-state/config)))
      :show-config (response/ok @everclear.state.state/config)
      ;; Not working reliably
      ;; :switch-config (response/ok
      ;;                  (reset! everclear.state.state/config
      ;;                          (assoc (helpers/file-to-map (:config options))
      ;;                            :config-file (:config options))))
      :dry-run (response/ok (dry-run jobgraph))
          (response/ok (str "No such command: " command)))))

(defn wrap-nocache [handler]
  (fn [request]
    (-> request
        handler
        (assoc-in [:headers "Pragma"] "no-cache"))))

(defn wrap-formats [handler]
  (-> handler
      (muuntaja/wrap-format)))

(defn home-page [request]
  (layout/render request "home.html"))

(defn home-routes []
  [""
   ["/" {:get home-page :middleware [middleware/wrap-csrf]}]
   ["/json" {:post (fn [req] (json-handler req)) :middleware [wrap-formats]}]
   ["/docs" {:get (fn [_]
                    (-> (response/ok (-> "docs/docs.md" io/resource slurp))
                        (response/header "Content-Type" "text/plain; charset=utf-8")))}]
   ["/ws" {:get  (:ajax-get-or-ws-handshake-fn ws/connection)
           :post (:ajax-post-fn ws/connection)
           :middleware [middleware/wrap-csrf]}] ])
