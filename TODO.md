# Have sample sheets and other wildcards. Do not want to keep it all in same table?

- No: a file is a very long thing, not suited to be a wildcard (takes up space)

# Pick out all {{}} from code/defrule/params and check that they all are present

If request "sample A", genome "hg19" must ensure that it is present in sample sheets etc

# Add source code highlighting

Pygments python

Clojure clojure

# Efficiency improvements

- avoid writing existing files
- check if an atom is changed before updating it

# Switch config does not work

The fn that watches is not updated when the config is switched?

# Should the outmap be a sorted map?

Would make everything easier.

Think: are all the filemaps sorted-maps? Not necessarily after subsetting and so on.

# Separate CLI into subcommands

# Store jobresults in atom?

Or otherwise, how can I remember which jobs failed?

# Use mnemonics instead of hashes

ranged-rabbit instead of 5cebb13e6f

# How to start job from command line?

Until?

Need to find all corresponding jobs to the rule in the jobgraph.

Then need to subset on the wildcards

# Perhaps it is best to only start job-runner when it is needed? I think so.

# Delete input for

# Investigate why selmer is not throwing error when field is missing

# Update d3-dag to most recent version:

Use plain javascript instead and get help to write it from Erik Brinkman

https://shadow-cljs.github.io/docs/UsersGuide.html#classpath-js
