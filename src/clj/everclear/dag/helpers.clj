(ns everclear.dag.helpers
  (:require
    [loom.alg]
    [clojure.edn :as edn]
    [clojure.java.shell :refer [sh]]
    [clojure.string :as str]
    [clojure.set :as set]
    [clojure.data.csv :as csv]
    [clojure.java.io :as io])
  (:use [hashp.core])
  (:import java.nio.file.Paths))

(defn make-absolute [basepath path]
  (if (.isAbsolute (io/file path))
    path
    (.toString
      (Paths/get
        basepath
        (into-array [path])))))

(defn read-delimited [f sep]
  (with-open [reader (io/reader f)]
    (doall
      (csv/read-csv reader :separator sep))))

(defn delim-file->maps [f sep]
  (let [delim-data (read-delimited f sep)]
    (map zipmap
         (->> (first delim-data) ;; First row is the header
              (map keyword) ;; Drop if you want string keys instead
              repeat)
         (rest delim-data))))

(defn file-to-map [f]
  ;; TODO: add support for more extensions
  (if-not (.exists (io/file f))
    (throw (Exception. (str "File does not exist: " f)))
    (let [ext (last (clojure.string/split f #"\."))]
      (case ext
        "edn" (-> (slurp f) edn/read-string)
        "clj" (-> (slurp f) edn/read-string)
        "tsv" (delim-file->maps f \tab)
        "csv" (delim-file->maps f \,)
        (throw (Exception. (str "Invalid extension for " f ". Must be edn, clj, tsv, csv.")))))))

(defn invert-many-to-one
  "Returns a one-to-many mapping"
  ([m] (invert-many-to-one #{} m))
  ([to m]
   (persistent!
    (reduce (fn [m [k v]]
              (assoc! m v (conj (get m v to) k)))
            (transient {}) m))))

(defn ->Path [p]
  (java.nio.file.Paths/get (java.net.URI. (str "file://" p))))

(defn delete-files! [fs]
  (map io/delete-file fs))

(defn file-existence-keyed-by-file [all-files]
    (into {} (for [f all-files
                   :when (not (nil? f))]
               [f (.exists (io/file f))])))

(defn create-folders! [fs]
  (doseq [f fs] (io/make-parents f)))

(defn files-to-be-created [jobinfo]
  (remove nil? (flatten (vals (:out-vec jobinfo)))))

(defn files-and-symlinks-to-be-created [job]
  (remove nil? (flatten (concat (vals (job :publish-vec)) (vals (job :out-vec))))))

(defn create-outfolders! [jobinfo]
  (create-folders! (files-and-symlinks-to-be-created jobinfo)))

(defn sort-wcmap-by-keys [wcs->val ks]
  (if (seq ks)
    ;; Don't use `first` and `second` with map entries - use `key` and `val` instead.
    (sort-by #((apply juxt ks) (key %)) wcs->val)
    {}))

(defn sort-vals-by-wcmap-keys [wcs->val]
  (map val
       (sort-wcmap-by-keys wcs->val
                           (sort (-> wcs->val first key keys)))))

(defn sort-wc-rows-by-ks [wc-rows]
  (if-let [ks (-> wc-rows first keys seq)]
    (sort-by (apply juxt ks) wc-rows)
    [{}]))

(defn symlink! [srcfile symlink]
  (java.nio.file.Files/createSymbolicLink
    (java.nio.file.Paths/get symlink (into-array String []))
    (java.nio.file.Paths/get srcfile (into-array String []))
    (into-array java.nio.file.attribute.FileAttribute [])))

(defn mv! [from to]
  (let [to-file (io/file to)]
    (io/make-parents to-file)
    (.renameTo (io/file from) to-file)))

(defn cp! [from to]
  (let [to-file (io/file to)]
    (io/make-parents to-file)
    (io/copy (io/file from) to-file)))

(defn is-binary? [f]
  (let [stdout (:out (clojure.java.shell/sh "file" "-I" "-b" f))
        charset (second (re-find #"charset=(.*)" stdout))]
    (= "binary" charset)))

(defn file-existence [outfile-exists publishfile-exists]
  (condp = [outfile-exists publishfile-exists]
    [true nil] :ok
    [true true] :ok
    [true false] :publish-missing
    [false true] :output-missing
    [false nil] :both-missing
    [false false] :both-missing))

(defn ->existence->outfile+publishfile [jobinfo file->exists?]
  (apply merge-with concat
         (for [[outalias outfiles] (:out-vec jobinfo)
               :let [publishfiles (get-in jobinfo [:publish-vec outalias] (repeat nil))]
               [outfile publishfile] (map vector outfiles publishfiles)]
           {(file-existence (file->exists? outfile) (file->exists? publishfile)) [[outfile publishfile]]})))

(defn outfile+publishfile [jobhash->jobinfo]
  (for [[_jobhash jobinfo] jobhash->jobinfo
        [output-alias publishfiles] (:publish jobinfo)
        [idx publish-file] (map vector (range) publishfiles)
        :let [output-file (get-in jobinfo [:output output-alias idx])]]
    [output-file publish-file]))

(defn symlink-missing-files! [existence->outfile+publishfile]
  (doseq [[outputfile publishfile] (:publish-missing existence->outfile+publishfile)]
    (symlink! outputfile publishfile))
  (doseq [[outputfile publishfile] (:output-missing existence->outfile+publishfile)]
    (do
      (mv! publishfile outputfile)
      (symlink! outputfile publishfile))))

(defn get-jobs-by-rule [rule jobhash->jobid jobhash->job]
  (let [correct-jobhashes (for [[k [rulename _wildcards]] jobhash->jobid
                                :when (= rulename rule)]
                            k)]
    (vals (select-keys jobhash->job correct-jobhashes))))

(defn topo-sort [rulegraph]
  (let [topo-sorted-rules (loom.alg/topsort rulegraph)]
    (into {}
          (map
            (fn [[a b]] [b a])
            (map-indexed vector topo-sorted-rules)))))

(defn topo-sort-mapseq-by-rulename [jobids rulegraph]
  (let [rulename+idx (topo-sort rulegraph)]
    (sort-by #(-> % :rulename rulename+idx) jobids)))

(defn topo-sort-jobids-by-rulename [mapseq rulegraph]
  (let [rulename+idx (topo-sort rulegraph)]
    (sort-by #(-> % first rulename+idx) mapseq)))

(defn outfile-maps [jobid->jobinfo]
  (for [[jobid job] jobid->jobinfo
        [alias ext] (:output-extensions job)
        :let [[rulename wildcards] jobid
              outfile-map (get-in job [:output-map alias])
              summary-map (get-in job [:summary-map alias])
              protopath-map (get-in job [:protopath-map alias])]]
    {:extension ext :jobid [rulename wildcards]
     :protopath-map protopath-map
     :alias alias :rulename rulename :wildcards wildcards
     :outfile-map outfile-map :summary-map summary-map}))