(ns everclear.state.state
  (:require [everclear.dag.helpers :as helpers]
            [clojure.java.io]
            [clojure.set :as set])
  (:use [hashp.core]))

;; should have different atoms for regular and hawk watchers?
(def watches
  (atom {}))

;; One reason to keep these in state: front-end might request them

(def jobinfo
  (atom {}))

(def resources (atom {}))

(def rulename->ruleinfo (atom {}))

(def rulegraph-map (atom {}))

(def jobid->jobinfo (atom {}))

(def error
  (atom nil))

(def jobgraph (atom nil))

(def rulegraph (atom nil))

;; how to send the exception further up?
(defn config-validator [{:keys [config-file] :as config}]
  (when (seq config) ;; do not validate when atom is empty, i.e. created
    (when-not (and (seq config-file) (.exists (clojure.java.io/file config-file)))
      (throw (ex-info (str "Config file does not exist: " config-file) {:invalid-state config}))))
    true)

(def config (atom {}))

;; Before we wanted to have a loop to avoid starting jobs in new threads.
;; But we learned the reason is that we used futures, not that we used watchers
;; Perhaps better design with having a watcher? Then if a job finishes can
;; just call blah blah
;;
;; Problem is that functions that change state should not be called in atoms
