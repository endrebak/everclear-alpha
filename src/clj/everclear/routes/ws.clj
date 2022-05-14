(ns everclear.routes.ws
  (:require
    [clojure.pprint :as pp]
    [everclear.jobrunning.run-jobs :as run]
    [everclear.state.state :as state]
    [mount.core :refer [defstate]]
    [taoensso.sente :as sente]
    [taoensso.sente.server-adapters.http-kit :refer (get-sch-adapter)]
    [everclear.state.state :as app-state]
    [everclear.graph.graph :as graph]))

(defn client-id [ring-req]
  (get-in ring-req [:params :client-id]))

(defstate connection
  :start (sente/make-channel-socket!
          (get-sch-adapter)
          {:packer :edn
           :user-id-fn client-id}))

(defn send-data-over-websockets [event-id data]
  (println (str "Sending to " event-id))
  (doseq [uid (:any @(:connected-uids connection))]
    ((:send-fn connection) uid [event-id data])))

(defn handle-message! [{:keys [id client-id ?data] :as message}]
  (case id
    :request/jobinfo (do
                       (println (str "in case :request/jobinfo :" id))
                       (swap! app-state/jobid->jobinfo identity))
    ;; :request/ruleinfo (do
    ;;                    (println (str "in case :request/jobinfo :" id))
    ;;                    (swap! state/ruleinfo identity))
    :request/rulegraph (send-data-over-websockets
                         :rulegraph/update
                         (graph/graph-data-js @everclear.state.state/rulegraph name))
    :commands/start-job (everclear.jobrunning.run-jobs/run-jobs!
                          @app-state/jobid->jobinfo identity @app-state/jobgraph @app-state/config)
    :chsk/ws-ping))

(defn stop-router! [stop-fn]
  (when stop-fn (stop-fn)))

(defn start-router! []
  (sente/start-chsk-router! (:ch-recv connection) #'handle-message!))

(defstate router
  :start (start-router!)
  :stop (stop-router! router))
