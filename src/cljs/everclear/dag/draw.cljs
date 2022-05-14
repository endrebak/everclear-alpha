(ns everclear.dag.draw
  (:require
   [cljs.core.async :refer [go]]
   [cljs.core.async.interop :refer-macros [<p!]]
   ["d3" :as d3]
   ["d3-dag" :as d3-dag]
   [re-frame.core :as rf]
   [clojure.set :as set]))

(defn p [arg]
  (js/console.log arg))

(defn draw-dag [rulejson]
  (let [_ (-> (.select d3 "#dag-graph") (.remove))
        reader (d3-dag/dagStratify)
        dag (reader (clj->js rulejson))
        node-radius 55
        padding 1.5
        base (* node-radius padding 2)
        node-size (fn [] #js [base base])
        layout (->
                (d3-dag/sugiyama)
                (.nodeSize #js [60 60])
                (.layering (d3-dag/layeringSimplex))
                (.decross (d3-dag/decrossOpt))
                (.coord (d3-dag/coordQuad))
                (.nodeSize node-size))
        new-dag ^js (layout dag)
        width (* 1 (.-width new-dag))
        height (* 1 (.-height new-dag))
        svg-selection (->
                       (d3/select "#wrapper")
                       (.append "svg")
                       (.attr "id" "dag-graph")
                       (.attr "width" width)
                       (.attr "height" height))
        defs (.append svg-selection "defs")
        steps (.size dag)
        interp d3/interpolateRainbow
        color-map (into {} (for [[i node] (.. dag idescendants entries)]
                             [node.data.id (interp (/ i steps))]))
        line (-> (.line d3)
                 (.curve (.-curveCatmullRom d3))
                 (.x (fn [d] (.-x d)))
                 (.y (fn [d] (.-y d))))
        _ (-> svg-selection
              (.append "g")
              (.selectAll "path")
              (.data (.links dag))
              (.enter)
              (.append "path")
              (.attr "d" (fn [obj]
                           (line (.-points obj))))
              (.attr "fill" "none")
              (.attr "stroke-width" 3)
              (.attr "stroke"
                     (fn [obj]
                       (let [source (.-source obj)
                             target (.-target obj)
                             source-data-id (-> source .-data .-id)
                             target-data-id (-> target .-data .-id)
                             grad-id (js/encodeURIComponent (str source-data-id "--" target-data-id))
                             grad (->
                                   (.append defs "linearGradient")
                                   (.attr "id" grad-id)
                                   (.attr "gradientUnits" "userSpaceOnUse")
                                   (.attr "x1" (.-x source))
                                   (.attr "x2" (.-x source))
                                   (.attr "y1" (.-y source))
                                   (.attr "y2" (.-y source)))
                             _ (-> grad
                                   (.append "stop")
                                   (.attr "offset" "0%")
                                   (.attr "stop-color" (color-map source-data-id)))
                             _ (-> grad
                                   (.append "stop")
                                   (.attr "offset" "100%")
                                   (.attr "stop-color" (color-map target-data-id)))
                             url (str "url(#" grad-id ")")]
                         url))))
        nodes (-> svg-selection
                  (.append "g")
                  (.selectAll "g")
                  (.data (.descendants dag))
                  (.enter)
                  (.append "g")
                  (.attr "transform" (fn [obj] (str "translate(" (+ (.-x obj) 10) ", " (.-y obj) ")"))))

        ;; _ (-> nodes
        ;;       (.append "rect")
        ;;       (.attr "width" 200)
        ;;       (.attr "height" 50)
        ;;       (.attr "fill" (fn [n] (color-map (-> n .-data .-id)))))

        _ (-> nodes
              (.append "text")
              (.text (fn [d] (-> d .-data .-id)))
              (.attr "font-weight" "bold")
              (.attr "font-family" "sans-serif")
              (.attr "text-anchor" "right")
              (.attr "alignment-baseline" "middle")
              (.attr "fill" "black"))]))

(defn load-dag [rulejson]
  (if (seq rulejson)
    (do
      (draw-dag rulejson)
      [:h3 "load-dag"])
    [:h3 "(No graph is drawn because the rule-data is empty)"]))

(defn graph-page []
  (when-let [rulegraph-map @(rf/subscribe [:rulegraph-map])]
      (when (seq rulegraph-map)
        [:div
         [:h1 "Hiya from graph-page"]
         [load-dag rulegraph-map]])))
