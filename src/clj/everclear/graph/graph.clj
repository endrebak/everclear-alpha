(ns everclear.graph.graph
  (:require [loom.graph]
            [loom.alg]
            [loom.derived]
            [clojure.set :as set]
            [everclear.jobcreation.wildcards :as wildcards]))

(defn unconnected-nodes [id->info parent+child]
  (set/difference (set (keys id->info))
                  (set (apply concat parent+child))))

(defn parent+child [id->info]
  (for [[child-id child-info] id->info
        parent-id (apply concat (-> child-info :parents vals))]
    [parent-id child-id]))

(defn graph [id->info]
  (let [parent+child (parent+child id->info)
        proto-graph (loom.graph/add-nodes*
                      (loom.graph/digraph)
                      (unconnected-nodes id->info parent+child))]
  (loom.graph/add-edges* proto-graph parent+child)))

(defn topo-sort [graph] (loom.alg/topsort graph))

(defn child->parents [graph] (:in graph))

(defn parent->children [graph] (:adj graph))

(defn dag? [graph] (loom.alg/dag? graph))

(defn nodes [graph] (loom.graph/nodes graph))

(defn lonely-nodes [graph] (loom.alg/loners graph))


(defn predecessors [gr node] (loom.graph/predecessors gr node))

(defn successors [gr node] (loom.graph/successors gr node))

(defn jobs-by-rulename [jobs rulename]
  (for [[job-rulename wildcards] jobs
        :when (= rulename job-rulename)]
    [job-rulename wildcards]))

(defn jobs-by-wildcards [jobs wildcards]
  (let [ks (keys wildcards)]
    (for [[job-rulename job-wildcards] jobs
          :when (wildcards/submap? (select-keys job-wildcards ks) wildcards)]
      [job-rulename job-wildcards])))

(defn- recurse-from-nodes-using-f [nodes node->nodes]
  (loop [unvisited (set (mapcat node->nodes nodes))
         visited (set nodes)]
    (let [new-unvisited (set
                          (apply concat
                                 (map node->nodes unvisited)))
          new-visited (set/union visited unvisited)]
      (if-not (seq unvisited)
        visited
        (recur new-unvisited new-visited)))))

;; We need to get all the dependencies of a job

(defn jobs
  ([jobgraph rulename] (jobs jobgraph rulename nil))
  ([jobgraph rulename wildcards]
   (cond-> (nodes jobgraph)
           rulename (jobs-by-rulename rulename)
           wildcards (jobs-by-wildcards wildcards))))

(defn all-ancestors
  ([graph children]
   (recurse-from-nodes-using-f children (child->parents graph))))

(defn all-successors
  ([graph from] (recurse-from-nodes-using-f from (parent->children graph))))

(defn change-nodes [f g]
  (loom.derived/mapped-by f g))

(defn non-orphans [g]
 (-> g child->parents keys set))

(defn parentless [g]
  (let [all-nodes (:nodeset g)]
    (set/difference all-nodes (non-orphans g))))

(defn graph-data-js
  ([g] (graph-data-js g identity))
  ([g f]
   (let [g (change-nodes f g)
         parentless (for [n (parentless g)] {:id n})
         child+parents (into #{} (for [[c ps] (child->parents g)]
                                   {:id c :parentIds ps}))]
     (concat child+parents parentless))))