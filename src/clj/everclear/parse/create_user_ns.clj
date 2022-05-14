(ns everclear.parse.create-user-ns
  (:require [everclear.parse.parse :as parse]))

(defn fetch-rules [namespace]
  (for [[_k v] (ns-publics namespace)
        :when (:rule (meta v))]
    (deref v)))

(defn var->val [form+string]
  (apply merge
         (for [[form form-s] form+string
               :let [n (second form)]
               :when (#{'defn 'def} (first form))]
           {n (clojure.edn/read-string form-s)})))

(defn update-ruleforms [ruleforms var->val]
  (doall (for [ruleform ruleforms]
           (if-let [input-fns (:input-fns ruleform)]
             (assoc ruleform
               :input-code (doall
                             (into {}
                                   (for [[alias func] input-fns
                                         :let [func-var-without-namespace (-> func class .getName clojure.repl/demunge symbol name symbol)]]
                                     [alias (get var->val func-var-without-namespace (get-in ruleform [:input-code alias]))]))))
               ruleform))))

(defmacro defrule
  ([rulename body] (defrule rulename body nil))
  ([rulename doc body]
   (let [input-fns (:input-fns body)
         _ input-fns
         body (-> (if doc (assoc body :doc doc) body)
                  (dissoc :input-fns)
                  (assoc :rulename (keyword rulename)))]
     (if input-fns
       `(def ~(with-meta rulename {:rule true})
          (assoc ~body
            :input-fns ~input-fns
            :input-code '~input-fns))
       `(def ~(with-meta rulename {:rule true}) ~body)))))

(defn eval-user-ns [file]
  ;; TODO (assert (= ns-name 'rules) "User-defined namespace must be called rules.")
  (remove-ns 'rules)
  (load-file file)
  (fetch-rules 'rules))

;; (defn eval-user-ns [defs]
;;   (let [prev-ns (-> *ns* str symbol)]
;;     (println defs)
;;     (println prev-ns)
;;     (ns user-ns)
;;     ;; Make the config available in the name-space
;;     (eval '(def config @everclear.state.state/config))
;;     (doseq [form defs] (eval form))
;;     (in-ns prev-ns)))

(defn eval-rules-in-user-ns [ruleforms ptr-list]
  (let [prev-ns (-> *ns* str symbol)]
    (ns user-ns)
    (doseq [form ruleforms]
      (swap! ptr-list conj (eval form)))
    (in-ns prev-ns)))