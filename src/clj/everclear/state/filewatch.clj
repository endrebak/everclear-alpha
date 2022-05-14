(ns everclear.state.filewatch
  (:require
   [hawk.core :as hawk]
   [everclear.state.state :as state]
   [everclear.routes.ws :as ws]))

(defn initiate-watch! [n fs fn]
  (swap! state/watches
         assoc n (hawk/watch! {:watcher :polling :sensitivity :high} [{:paths fs :handler fn}])))

(defn end-watch! [n]
  (hawk/stop! (@state/watches n))
  (swap! state/watches dissoc n))

(defn update-watch! [n fs fn]
  (end-watch! n)
  (initiate-watch! n fs fn))

(add-watch
 state/error
 :error-watcher
 (fn [_ _ _ new-state]
   (doseq [uid (:any @(:connected-uids ws/connection))]
     ((:send-fn ws/connection) uid [:jobinfo/error new-state]))))
