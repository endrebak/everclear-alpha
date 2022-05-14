(ns everclear.parse.parse
  (:import [java.io PushbackReader]
           (clojure.lang LineNumberingPushbackReader))
  (:require
   [clojure.java.io :as io]))

(defn read-forms-in-file [file]
  (let [rdr (-> file io/file io/reader PushbackReader.)]
    (loop [forms []]
      (let [form (try (read rdr) (catch Exception e nil))]
        (if form
          (recur (conj forms form))
          forms)))))

(defn read+string-forms-in-file [file]
  (let [rdr (-> file io/file io/reader LineNumberingPushbackReader.)]
    (loop [forms []]
      (let [form (try (read+string rdr) (catch Exception e nil))]
        (if form
          (recur (conj forms form))
          forms)))))

(defn ->forms [rule-file]
  ;; split into two: defrule and the rest
  ;; the first becomes maps, the second evaled in a namespace
  (let [forms (read-forms-in-file rule-file)
        {rules true defs false} (group-by #(= 'defrule (first %)) forms)]
    {:rules rules :defs defs}))

(defn files->namedlist [files]
  (cond
    (or (map? files) (nil? files)) files
    (vector? files) (into {} (map-indexed vector files))
    (string? files) {0 files}
    :else (throw (Exception. ":input/:output must be map, string or vector"))))

(defn add-doc
  [[_ name & body]]
  (let [name (keyword name)]
    (if (= 2 (count body))
      (assoc (second body) :doc (first body) :rulename name)
      (assoc (first body) :rulename name))))

(defn input-and-output-to-namedlist [rule]
  (let [output (files->namedlist (rule :output))
        input (files->namedlist (rule :input))
        files
        {:input input
         :output output}]
    (merge rule
           (if-let [publish (rule :publish)]
             (assoc files :publish (files->namedlist publish))
             files))))

(defn add-default-resources [rulemap]
  (merge-with merge {:resources {:cpu 1}} rulemap))

(defn ->ruleinfo [ruleforms]
  (let [rules ruleforms]
    (for [rule rules
          :let [rulemap (-> rule
                            input-and-output-to-namedlist
                            add-default-resources)]]
      rulemap)))

(defn ruleinfo-keyed-by-rulename [ruleforms]
  (let [rulemaps (->ruleinfo ruleforms)]
    (into {} (for [rulemap rulemaps
                   :let [name (:rulename rulemap)]]
               [name rulemap]))))
