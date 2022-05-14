(ns everclear.core
  (:require
    [day8.re-frame.http-fx]
    [reagent.dom :as rdom]
    [reagent.core :as r]
    [re-frame.core :as rf]
    [goog.events :as events]
    [goog.history.EventType :as HistoryEventType]
    [markdown.core :refer [md->html]]
    [everclear.ajax :as ajax]
    [everclear.ws :as ws]
    [everclear.dag.draw :as dag]
    [everclear.events]
    [reitit.core :as reitit]
    [reitit.frontend.easy :as rfe]
    [clojure.string :as str])
  (:import goog.History))

(defn nav-link [uri title page]
  [:a.navbar-item
   {:href   uri
    :class (when (= page @(rf/subscribe [:common/page-id])) :is-active)}
   title])

(defn navbar []
  (r/with-let [expanded? (r/atom false)]
              [:nav.navbar.is-info>div.container
               [:div.navbar-brand
                [:a.navbar-item {:href "/" :style {:font-weight :bold}} "everclear"]
                [:span.navbar-burger.burger
                 {:data-target :nav-menu
                  :on-click #(swap! expanded? not)
                  :class (when @expanded? :is-active)}
                 [:span][:span][:span]]]
               [:div#nav-menu.navbar-menu
                {:class (when @expanded? :is-active)}
                [:div.navbar-start
                 [nav-link "#/" "Home" :home]
                 [nav-link "#/about" "About" :about]]]]))

(defn about-page []
  [:section.section>div.container>div.content
   [:img {:src "/img/warning_clojure.png"}]])

(defn format-files [header filemap]
  (if-not (seq filemap)
    nil
    [:div
     [:h4 header]
     (for [[alias jobfiles] filemap]
       ^{:key jobfiles}
       [:div
        [:h6 alias]
        (for [[xp-wcs file] jobfiles]
          ^{:key file}
          [:p file])])]))

(defn hiya []
  (when-let [jobinfo (rf/subscribe [:jobinfo])]
    [:div
     (for [[jobid jobinfo] @jobinfo
           :let [[rulename wildcards] jobid]]
       ^{:key jobid}
       [:div
        [:h2 (name rulename)]
        [:h4 (str/join " "
                       (apply concat
                              (sort-by first
                                       (for [[k v] wildcards]
                                         ^{:key v}
                                         [(name k) v]))))]
        [:div
         [:div (format-files "Externals" (:externals-map jobinfo))]
         [:div (format-files "Input" (:input-map jobinfo))]
         [:div (format-files "Output" (:output-map jobinfo))]
        [:pre (or (jobinfo :rendered-script) (jobinfo :rendered-shell))]
        [:br]]])]))


(defn error-div []
  (let [{:keys [trace cause]} @(rf/subscribe [:jobinfo-error])]
    (if (seq trace)
      (let [more-relevant-errors (filter #(->> % (re-matches #".*\[(user|clojure).*") not) trace)]
        [:div
         {:style {:position "fixed"
                  :z-index "10000"
                  :left "0px"
                  :bottom "0px"
                  :right "0px"
                  :background-color :red
                  :font-size 24
                  :color :white}}
         [:div {:style
                ;; Does not actually center the error div
                {:margin "auto"}}
          [:div {:style {:font-size 32}} cause]
          (for [err more-relevant-errors]
            ^{:key err}
            [:div {:style {:font-size 24}}
             err
             [:br]])]])
      [:div#error ""])))


(defn home-page []
  (js/console.log "requesting data from server")
  (ws/send-message! [:request/rulegraph])
  ;; (ws/send-message! [:request/ruleinfo])
  (let [];; jobinfo @(rf/subscribe [:jobinfo])]
    [:section.section>div.container>div.content
     [:input.button.is-primary
      {:type     :submit
       :on-click #(do
                    (ws/send-message! [:commands/start-job :whichever] 0))
       :value    "Run all jobs"}]
     [:div#wrapper]
     [:img {:src "/img/warning_clojure.png"}]
     [:br]
     [:f> dag/graph-page]
     [:f> error-div]
     [hiya]
     ;; (when-let [docs @(rf/subscribe [:docs])]
     ;;   [:div {:dangerouslySetInnerHTML {:__html (md->html docs)}}])
     ]))

(defn page []
  (if-let [page @(rf/subscribe [:common/page])]
    [:div
     [navbar]
     [page]]))

(defn navigate! [match _]
  (rf/dispatch [:common/navigate match]))

(def router
  (reitit/router
    [["/" {:name        :home
           :view        #'home-page
           :controllers [{:start (fn [_] (rf/dispatch [:page/init-home]))}]}]
     ["/about" {:name :about
                :view #'about-page}]]))

(defn start-router! []
  (rfe/start!
    router
    navigate!
    {}))


(defn response-handler [messages fields errors]
  (fn [{[_ message] :?data}]
    (if-let [response-errors (:errors message)]
      (reset! errors response-errors)
      (do
        (reset! errors nil)
        (reset! fields nil)
        (swap! messages conj message)))))

(def session (r/atom {}))

;; -------------------------
;; Initialize app
(defn ^:dev/after-load mount-components []
  (rf/clear-subscription-cache!)
  (rdom/render [#'page] (.getElementById js/document "app")))

(defn init! []
  (ws/start-router!
   (response-handler
    (r/cursor session [:messages])
    (r/cursor session [:fields])
    (r/cursor session [:errors])))
  (start-router!)
  (ajax/load-interceptors!)
  (mount-components))
