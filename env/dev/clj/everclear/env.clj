(ns everclear.env
  (:require
    [selmer.parser :as parser]
    [clojure.tools.logging :as log]
    [everclear.dev-middleware :refer [wrap-dev]]))

(def defaults
  {:init
   (fn []
     (parser/cache-off!)
     (log/info "\n-=[everclear started successfully using the development profile]=-"))
   :stop
   (fn []
     (log/info "\n-=[everclear has shut down successfully]=-"))
   :middleware wrap-dev})

;; (alter-var-root #'clojure.core/get-in
;;                 (fn [v]
;;                   (fn [m ks]
;;                     (let [v (v m ks ::not-found)]
;;                       (if (= v ::not-found)
;;                         (throw (ex-info "" {}))
;;                         v)))))
