(ns everclear.ws
  (:require
   [re-frame.core :as rf]
   [taoensso.sente  :as sente :refer (cb-success?)]))

(let [connection (sente/make-channel-socket! "/ws" js/csrfToken {:type :auto :wrap-recv-evs? false})]
  (def ch-chsk (:ch-recv connection)) ; ChannelSocket's receive channel
  (def send-message! (:send-fn connection)))

(defmulti -event-msg-handler
  "Multimethod to handle Sente `event-msg`s"
  :id ; Dispatch on event-id
  )

(defmethod -event-msg-handler :jobinfo/update
  [{:as ev-msg :keys [?data]}]
  (js/console.log (str "jobinfo"))
  (rf/dispatch [:save-jobinfo ?data]))

(defmethod -event-msg-handler :rulegraph/update
  [{:as ev-msg :keys [?data]}]
  (js/console.log (str "rulegraph received" ?data))
  (rf/dispatch [:save-rulegraph-map ?data]))

(defmethod -event-msg-handler :ruleinfo/update
  [{:as ev-msg :keys [?data]}]
  ;; (js/console.log (str "ruleinfo " (with-out-str (cljs.pprint/pprint (map (juxt :in-rule :out-rule ?data))))))
  (rf/dispatch [:save-ruleinfo ?data]))

(defmethod -event-msg-handler :jobinfo/error
  [{:as ev-msg :keys [?data]}]
  ;; (js/console.log (str "correctamundo " ?data))
  (rf/dispatch [:save-jobinfo-error ?data]))

(defmethod -event-msg-handler :reply/hi
  [{:as ev-msg :keys [?data]}]
  (js/console.log (str ":replya/hia " ?data)))

(defmethod -event-msg-handler
  :default ; Default/fallback case (no other matching handler)
  [{:as ev-msg :keys [event]}]
  (println "Unhandled event: %s" event))

(defn state-handler [{:keys [?data]}]
  (.log js/console (str "state changed: " ?data)))

(defn handshake-handler [{:keys [?data]}]
  (.log js/console (str "connection established: " ?data)))

(defn default-event-handler [ev-msg]
  ;; (js/console.log "In default-event-handler")
  (.log js/console (str "Unhandled event: " (:event ev-msg))))

(defn event-msg-handler [& [{:keys [message state handshake]
                             :or {state state-handler
                                  handshake handshake-handler}}]]
  #(-event-msg-handler %))

(def router (atom nil))

(defn stop-router! []
  (when-let [stop-f @router] (stop-f)))

(defn start-router! [message-handler]
  (stop-router!)
  (reset! router (sente/start-chsk-router!
                   ch-chsk
                   (event-msg-handler
                     {:message   message-handler
                      :state     state-handler
                      :handshake handshake-handler}))))
