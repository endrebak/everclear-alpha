(ns everclear.events
  (:require
    [re-frame.core :as rf]
    [ajax.core :as ajax]
    [reitit.frontend.easy :as rfe]
    [reitit.frontend.controllers :as rfc]))

;;dispatchers

(rf/reg-sub
 :rulegraph-map
 (fn [db _]
   ;; (js/console.log (str "rulegraph triggered: " (-> db :rulegraph)))
   (-> db :rulegraph-map)))

(rf/reg-event-db
 :save-rulegraph-map
 (fn [db [_ rulegraph-map]]
   ;; (js/console.log (str "rulegraph updated: " rulegraph))
   (assoc db :rulegraph-map rulegraph-map)))

(rf/reg-sub
 :ruleinfo
 (fn [db _]
   ;; (js/console.log (str "ruleinfo sub triggered"))
   (-> db :ruleinfo)))

(rf/reg-event-db
 :save-ruleinfo
 (fn [db [_ ruleinfo]]
   ;; (js/console.log "ruleinfo updated: ")
   (assoc db :ruleinfo ruleinfo)))

(rf/reg-sub
 :jobinfo
 (fn [db _]
   (js/console.log (str "jobinfo sub triggered"))
   ;; (js/console.log (clj->js (-> db :jobinfo)))
   (-> db :jobinfo)))


(rf/reg-event-db
 :save-jobinfo
 (fn [db [_ jobinfo]]
   (js/console.log "jobinfo updated: " jobinfo)
   (assoc db :jobinfo jobinfo)))

(rf/reg-event-db
 :save-jobinfo-error
 (fn [db [_ error]]
   ;; (js/console.log (str "error updated: " error))
   (assoc db :jobinfo-error error)))

(rf/reg-sub
 :jobinfo-error
 (fn [db _]
   ;; (js/console.log (str "jobinfo-error sub triggered"))
   ;; (js/console.log (clj->js (-> db :jobinfo)))
   (-> db :jobinfo-error)))



(rf/reg-event-db
  :common/navigate
  (fn [db [_ match]]
    (let [old-match (:common/route db)
          new-match (assoc match :controllers
                                 (rfc/apply-controllers (:controllers old-match) match))]
      (assoc db :common/route new-match))))

(rf/reg-fx
  :common/navigate-fx!
  (fn [[k & [params query]]]
    (rfe/push-state k params query)))

(rf/reg-event-fx
  :common/navigate!
  (fn [_ [_ url-key params query]]
    {:common/navigate-fx! [url-key params query]}))

(rf/reg-event-db
  :set-docs
  (fn [db [_ docs]]
    (assoc db :docs docs)))

(rf/reg-event-fx
  :fetch-docs
  (fn [_ _]
    {:http-xhrio {:method          :get
                  :uri             "/docs"
                  :response-format (ajax/raw-response-format)
                  :on-success       [:set-docs]}}))

(rf/reg-event-db
  :common/set-error
  (fn [db [_ error]]
    (assoc db :common/error error)))

(rf/reg-event-fx
  :page/init-home
  (fn [_ _]
    {:dispatch [:fetch-docs]}))

;;subscriptions

(rf/reg-sub
  :common/route
  (fn [db _]
    (-> db :common/route)))

(rf/reg-sub
  :common/page-id
  :<- [:common/route]
  (fn [route _]
    (-> route :data :name)))

(rf/reg-sub
  :common/page
  :<- [:common/route]
  (fn [route _]
    (-> route :data :view)))

(rf/reg-sub
  :docs
  (fn [db _]
    (:docs db)))

(rf/reg-sub
  :common/error
  (fn [db _]
    (:common/error db)))
