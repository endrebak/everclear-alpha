This tutorial is a modified version of the Snakemake documentation. Its purpose is to explain how Everclear works. Because we copy their example, readers can juxtapose our guide with the Snakemake tutorial, which might be familiar.

Remember that we wrote Everclear to make medium-sized or larger workflows easier to write and maintain, so for this specific example, Everclear requires more setup. First, however, it serves to explain the basic concepts of Everclear. Second, it shows off the interactive Everclear software and how it helps development.

Let us begin.

Step 1: Mapping reads.

(defrule bwa-map
"Map DNA sequences against a reference genome with BWA."
{:externals  [:genome :fastq]
:output    "bwa-map.bam"
:shell     "bwa mem {{externals.genome}} {{externals.fastq}} |
samtools view -Sb - > {{output}}"})

An Everclear rule has a name, here bwa-map, and several entries called directives that specify the code the rule should run, the files it should produce, and the files it needs to execute.

The external directive specifies the external files the rule needs to run. The name external makes clear that these files are not produced by any rules in the Everclear workflow. We see that the entries in the external directive are keywords that map to specific input files. There are two different external files in this example, named genome and fastq.

The output directive is different; it specifies the output file's basename.

The shell directive contains the command to execute. Here one needs to refer to the input and output files one wants to use. We refer to the external input with {{externals}} and the output with {{output}}.

We store the mapping between keywords and external files in an additional file. Let us call it externals.edn:

{:genome "/Users/endrebakkenstovner/everclear/snakemake-flow/data/genome.fa"
:fastq "/Users/endrebakkenstovner/everclear/snakemake-flow/data/samples/A.fastq"}

We must also create a configuration that contains all the information Everclear needs to run the workflow.

Our configuration file looks like the following:

{:rule-file "examples/tutorial/step1/rules.clj"
:externals-file "examples/tutorial/step1/externals.edn"
:wildcards-file "examples/tutorial/step1/wildcards.edn"
:base-path "/Users/endrebakkenstovner/everclear/tutorial/step1/results"
:log-path "/Users/endrebakkenstovner/everclear/tutorial/step1/logs"}

The rule file contains the rule definitions. In our case, it only contains one rule, namely bwa-map. The externals file contains, as we saw, the mapping of keywords to files. The new entries are base-path and log-path. The base path describes where to output the files produced by the pipeline. The log path is where to store the logs.

Now use "everclear --start --config examples/tutorial/step1/config.edn &" to start Everclear. The Everclear software needs a few seconds to start, but this only has to be done once. Since the Everclear program is continually running, we start it in the background by adding "&" to the end of the command. Adding & gives us back control of the terminal.

Everclear will now serve a webpage. You can request it in the browser by writing localhost:3000 in the address bar. Using the browser dashboard is optional, but it does aid development by providing information not possible to display in the text-only terminal.

If you look at the webpage, you should see the job at the bottom.

We could run this job from the browser by clicking the "run all jobs" button. However, we will run it from the command line with "everclear --run".

If you now look at the output folder (given in the config file), you will see a subfolder named by the hash of the job, for example, 16589cf56b. In it, you will find the file bwa-map/bwa-map.bam.

Unless we delete the file, Everclear will never try to create this exact file again. However, if we change the rule ever so slightly, the job will be given a new hash, which does not exist on the file system.

Let us switch the order of the external keys, so fastq comes before genome. Save the file. Now we see that the job has gotten a new hash. If we call "everclear --run" again, Everclear runs the second job. This hashing of the job info ensures that files can never become stale and the different outputs of the pipeline out of sync.

Also, did you notice that when you saved the rule file after changing it, Everclear automatically updated the job? This reloading happens because Everclear watches the files referred to in the config file and recreates the DAG and job information upon any change. Neat, huh?
