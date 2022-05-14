(ns everclear.graph.jobgraph
  (:require [everclear.graph.graph :as g]
            [everclear.jobcreation.wildcards :as wildcards]))

(defn keep-jobs-by-rulename [jobgraph rulename]
    (for [[job-rulename wildcards] (g/nodes jobgraph)
          :when (= rulename job-rulename)]
      [job-rulename wildcards]))

(defn keep-jobs-by-wildcards [nodes wildcards]
  (for [[job-rulename job-wildcards] nodes
        :when (wildcards/submap? job-wildcards wildcards)]
    [job-rulename job-wildcards]))

(defn until
  ([jobgraph rulename] (until jobgraph rulename {}))
  ([jobgraph rulename wildcards]
   (-> (keep-jobs-by-rulename jobgraph rulename)
       (keep-jobs-by-wildcards wildcards))))
