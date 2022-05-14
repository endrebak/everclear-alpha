(ns everclear.core
  (:require
   [clojure.java.io :as io]
   [everclear.state.filewatch]
   [everclear.handler :as handler]
   [everclear.nrepl :as nrepl]
   [everclear.dag.core :as dag]
   [everclear.state.taps :as tap]
   [luminus.http-server :as http]
   [selmer.util]
   [everclear.config :refer [env]]
   [clojure.tools.cli :refer [parse-opts]]
   [clojure.tools.logging :as log]
   [mount.core :as mount])
  (:gen-class))

;; log uncaught exceptions in threads
(Thread/setDefaultUncaughtExceptionHandler
  (reify Thread$UncaughtExceptionHandler
    (uncaughtException [_ thread ex]
      (log/error {:what :uncaught-exception
                  :exception ex
                  :where (str "Uncaught exception on" (.getName thread))}))))

(def cli-options
  ;; An option with a required argument
  [["-c" "--config-file FILE" "Config file to use."
    :validate [#(and
                 (not (nil? %))
                 (.exists (io/file %)))]]
   ["-p" "--port PORT" "Port number"
    :parse-fn #(Integer/parseInt %)]
   ;; A boolean option defaulting to nil
   ["-h" "--help"]])

(mount/defstate ^{:on-reload :noop} http-server
  :start
  (http/start
    (-> env
        (assoc  :handler (handler/app))
        (update :port #(or (-> env :options :port) %))
        (select-keys [:handler :host :port])))
  :stop
  (http/stop http-server))

(mount/defstate ^{:on-reload :noop} repl-server
  :start
  (when (env :nrepl-port)
    (nrepl/start {:bind (env :nrepl-bind)
                  :port (env :nrepl-port)}))
  :stop
  (when repl-server
    (nrepl/stop repl-server)))

(defn stop-app []
  (doseq [component (:stopped (mount/stop))]
    (log/info component "stopped"))
  (shutdown-agents))

(defn start-app [args]
  (println "in start app")
  (Thread/setDefaultUncaughtExceptionHandler
   (reify Thread$UncaughtExceptionHandler
     (uncaughtException [_ thread ex]
       (println (str ex " Uncaught exception on " (.getName thread))))))
  (doseq [component (-> args
                        mount/start-with-args
                        :started)]
    (log/info component "started"))
  (.addShutdownHook (Runtime/getRuntime) (Thread. stop-app)))

(defn start-with-map [opts]
  (let [opts (if-not (:options opts)
               {:options opts}
               opts)]
    (start-app opts)
    (dag/start-everclear opts)
    (do
      ;; TODO: find better place for selmer-settings
      (selmer.util/set-missing-value-formatter!
       (fn [tag context-map]
         (let [missing-path (mapv keyword (clojure.string/split (:tag-value tag) #"\."))
               [missing-key missing-in-map] (last (for [i (range 1 (inc (count missing-path)))
                                                        :let [missing-path-subsequence (take i missing-path)
                                                              m (get-in context-map missing-path-subsequence)]
                                                        :when m]
                                                    [(last missing-path-subsequence) m]))]
           (throw (Exception.
                   (str "Missing key " missing-key " in tag " (tag :tag-value) " of map " missing-in-map " for rule " (context-map :rulename)))))))
      (Thread/sleep 1000)
      (println "****************adding tap****************")
      (add-tap #'tap/tap)
      ;; (try
      ;;   (future (run/run-jobs!))
      ;;   (catch Throwable t
      ;;     (let [{:keys [trace cause]} (Throwable->map t)]
      ;;       (tap> {:tap-id :error :tap-data {:cause cause :trace (map str trace)}}))))
      )))


(defn -main [& args]
  (let [opts (parse-opts args cli-options)]
    (start-with-map opts)))
