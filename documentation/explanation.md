Explanation of how Everclear works, by first explaining Snakemake and then comparing Everclear to it

A high-level overview of how Snakemake workflows are defined

Snakemake is a WMS inspired by GNU Make. Like GNU Make, the Snakemake language allows you to define a series of computational steps, called rules, where each step might take one or more inputs and produces one or more outputs. The inputs and outputs are paths to files. Inputs may either be files that already exist or files that some other rule can create as their output. So the users implicitly define the sequence of steps to be performed by ensuring step i's output path matches the input path of step i+1.

For example, one might have a simple workflow with two rules: download_genome, which downloads a genome from the internet, and split_genome, which splits the genome into one file per chromosome. The rule download_genome has no inputs and one output: the downloaded genome fasta. On the other hand, the rule split_genome has the just mentioned genome fasta as its input, and the output is multiple chromosome files.

The above description contains a simplification: rules are templates and can instantiate multiple jobs. So a rule is a function. And a job is a call to that function with a specific input.

How does Snakemake make it possible to instantiate rules into jobs? By adding one or more named variables, called wildcards, to the rule. To make the rules download_genome and split_genome work on different genomes, we must add a genome wildcard. This wildcard must be present in every output path (otherwise, the different instantiations of the jobs would overwrite each other's output). So for the rule download_genome, the output path could be "downloads/{genome}.fa". Each instantiation of the rule with different wildcards, like hg19 or mm10, would produce a different output file, like downloads/hg19.fa or downloads/mm10.fa.

What determines what values the wildcard genome will take? Proximately, other rules in the workflow decide it. For example, if some job, like split_genome, has genomes/hg19.fa as its input, it will ask the rule download_genome to output the fasta before starting its execution. However, the question is the same: what determines which values the wildcards in split_genome will take? Ultimately, the writer or invoker of the workflow decides the wildcard value.

So, in essence, the sequence of steps a Snakemake workflow should perform is determined by matching a set of strings against each other. Each string represents either the input or output of a rule in the workflow and can have one or more wildcards.

How does Snakemake execute its workflows?

A Snakemake workflow is executed by requesting a file. (This request is called a target in Make lingo). An example target might be "sorted/hg19/A/chr1.bam". If the file exists, Snakemake does nothing. If the file does not exist, Snakemake looks for a rule with an output template that matches the file. For example, the rule samtools_sort has the output string "sorted/{genome}/{sample}/{chromosome}.bam". This path matches the file "sorted/hg19/A/chr1.bam". Therefore, the job with wildcards genome hg19, sample A, and chromosome chr1 can create the file.

If the inputs to samtools_sort exist, for example, "bwa-map/hg19/A/chromosome.bam", Snakemake will run the job. If not, Snakemake will look for a rule with an output path that matches "bwa-map/hg19/A/chromosome.bam".

So Snakemake looks for a file, and if it does not find it, it looks for matching outputs in a rule.

Snakemake command line

The user starts Snakemake by requesting a file. Then Snakemake starts executing the steps described above until it finds the directed acyclic graph of jobs that can produce the requested file.

To check whether the pipeline is valid, one needs to invoke the Snakemake command-line tool. One must request every target the workflow can produce.

The Snakemake software and execution environment

Snakemake is a Python command-line script. It can run one command for each invocation.

The Snakemake workflow definition language

The language is a Python DSL (domain-specific language), which means that it uses the Python runtime and can embed and execute arbitrary Python code.

Logging in Snakemake

Snakemake supports a logging directive to write logs. It works like the output directive, but with one difference: Snakemake does not delete log files after a job fails. Keeping log files might help debug the reason a job failed.

Snakemake does not write any logs for the user, but users can add code to write logs to the rule.

How Everclear workflows are defined

Like in Snakemake, the building blocks of Everclear workflows are rules, and the links of input and output paths determine the sequence of steps to be performed. Furthermore, Everclear rules can contain wildcards so that a single rule can spawn multiple jobs.

However, Everclear unbraids wildcards and paths and makes them two separate fields. This Snakemake rule

rule samtools_sort:
input:
"mapped_reads/{sample}.bam"
output:
"sorted_reads/{sample}.bam"
shell:
"samtools sort -T sorted_reads/{wildcards.sample} "
"-O bam {input} > {output}"

would look like this in Everclear:

(defrule samtools-sort
"Sort the bams."
{:wildcards [:sample :genome]
:input     "bwa-map.bam"
:output    "bam/sorted.bam"
:publish   "{{genome}}/{{sample}}/sorted.bam"
;; TODO must ensure output and publish has the same keys!
:shell     "samtools sort -T {{wildcards.sample}} -O bam {{input.0}} > {{publish.0}}"})

Notice how the inputs and outputs defined in the Everclear rule can be short, like "bam/sorted.bam", or just a basename, like "bwa-map.bam". So, to find the DAG in Everclear, only simple string matching is needed; one does not have to pattern match against wildcards.

This leaves the question, how are the correct wildcards determined? Instead of requesting a path and matching that against a wildcard string, we use the sample sheet. A sample sheet is a structured spreadsheet, like .tsv, that includes the information needed to set up and analyze an experiment. Their first use-case was sequencing experiments, but they are well suited to describe the relationship between variables in any experiment.

For this little example, the below sample sheet suffices to illustrate:

genome sample p-value
hg19 A 0.01
hg19 B 0.01
hg19 A 0.05
hg19 B 0.05

Here we see three variables; only the two "genome" and "sample" are relevant for the rule samtools-sort. So, for this rule, we remove the column p-value and then find the unique rows left; these are hg19 A and hg19 B. So, the sample sheet and rule definition together determine that this job should instantiate two jobs, hg19 A and hg19 B.

Other differences between Snakemake and Everclear

Snakemake does not distinguish regular inputs from those the user should supply themselves. Everclear, unlike Snakemake, uses a specific external directive for inputs the workflow depends upon but cannot produce itself.

Everclear's logging capabilities.

Executing an Everclear workflow

Everclear computes the DAG and whether all files expected as output from any rule in the DAG exist. If some rules are missing their output, it marks those rules as todo.

Everclear command line

Whereas Snakemake receives a file request like sorted/hg19/A/chr1.bam and builds the DAG of jobs needed to produce it, Everclear can take a rule as an input, i.e., samtools-sort. When Everclear receives a rule as a target, it will start a job for each of its possible wildcard combinations. The jobs to run can be restricted by restricting the wildcards with an additional flag. In summary, in Snakemake, one needs to introduce multiple file targets to start multiple target jobs. However, in Everclear, a file is never the target, but rules are. Moreover, by default, Everclear starts all jobs possible for that rule.

Everclear logging

Everclear includes comprehensive logging by default. It stores the code run by the job and all information in the job data structure. Furthermore, it stores information about all the files produced by every job, like file size and, if the file is not binary, the head and the tail, and the number of lines in the file. It also includes a user-extendable library of file loggers.

This library maps file extensions to a script that produces logs for that file type. For example, the script that writes logs for bam files calls samtools flagstat -O tsv to produce a spreadsheet of alignment information. The user may also add scripts for unknown file types.

Unlike Snakemake, Everclear does not overwrite logs since it uses the job hash to give each run of the same job a different output path.

Advantages of the Everclear model

Source and target rules are knowable from the workflow definition alone.

Since the sequence of simple, wildcard-less strings define the DAG in Everclear, it is possible to know which rules are sources and which rules are targets from the workflow definition itself.

For example, take the simple Everclear pseudocode workflow definition with two rules, bwa-map and samtools-sort.

bwa-map
input: "sample.fastq"
output: "bwa-map.bam"

samtools-sort
input: "bwa-map.bam"
output: "sorted.bam"

Here, we can tell that bwa-map comes first and samtools-sort second. For Snakemake, a simplified workflow might look like

bwa-map
input: "{sample}.fastq"
output: "{sample}/{genome}/bwa-map.bam"

samtools-sort
input: "{sample}/{genome}/bwa-map.bam"
output: "{sample}/{genome}/sorted.bam"

While Snakemake could, in theory, derive the sequence of steps in this simple case, it does not do so. Furthermore, it cannot find the DAG in all cases. For comprehensive workflows with many rules, it would be computationally infeasible. It might also be theoretically impossible since inputs to Snakemake rules can be functions that produce a path based on the wildcards it receives.

The advantages of knowing targets from the workflow definition alone are significant in terms of understandability. For one, a workflow might have hundreds of rules and be hard to understand. However, if a workflow can know and present the user with all its targets, it becomes self-documenting. Indeed, this knowability makes it possible for Everclear to display the implied DAG of its rules quickly and continually.

In Snakemake, one can denote targets with rules that have no output. So, for example, one could create a target to produce all the sorted bam files like so:

rule sorted-bams:
input:
expand("{sample}/{genome}/sorted.bam",
sample=["A", "B"], genome="hg19")

However, this is additional work for the user, which is easy to forget, requires maintenance if the output path of the rule samtools-sort changes, and leads to an unnecessarily large code-base.

No DAG inconsistencies are possible.

In Snakemake, some inconsistencies might pop up that are not visible until the wildcards are declared, like recursive definitions and two rules being able to produce the same file.

In addition, Snakemake cannot find or warn about certain classes of hard-to-debug DAG errors. For example, the user might have changed the output of a rule A from {A}/{B}/{C}/{D}/{E}.txt to {A}/{B}/{C}/{D}/{E}/{F}.txt. In this example, the user forgot to change the input path of rule B, which depended on the output of rule A. That is, rule B's input is still {A}/{B}/{C}/{D}/{E}.txt. The two rules are not connected anymore, but if the user runs the pipeline in dry-run mode to check its integrity, Snakemake might say everything is well. Snakemake will consider the workflow intact if the user had previously run the pipeline and produced {A}/{B}/{C}/{D}/{E}.txt. The file then exists, so rule B seems satisfiable. This error in judgment happens because Snakemake does not distinguish between internal and external inputs, as discussed previously.

Less complexity in path names

Snakemake braids together paths and wildcards. This braiding leads to several problems. We just mentioned one: finding the DAG from the workflow definition alone becomes hard.

Furthermore, every file the workflow should produce requires the user to specify the exact location the file should have. Many of these will be intermediate results the user cares little about; therefore, requiring the user to specify the intended full path of every file is an unnecessary burden. Furthermore, it precludes some very desirable features, like having the paths include metadata.

Everclear, on the other hand, produces the output paths for the user. Take the example Everclear rule below:

(defrule samtools-sort
"Sort the bams."
{:wildcards [:sample :genome]
:input     "bwa-map.bam"
:output    "bam/sorted.bam"
:publish   "{{genome}}/{{sample}}/sorted.bam"
:shell     "samtools sort -T {{wildcards.sample}} -O bam {{input.0}} > {{publish.0}}"})

Here, the output is "bam/sorted.bam", and the wildcards are genome and sample. With our sample sheet, these wildcards give rise to two different jobs: hg19 A and hg19 B. Everclear will produce two different output paths for these two different jobs: c335dfb384/samtools-sort/hg19/A/bam/sorted.bam and 3e889ac412/samtools-sort/hg19/B/bam/sorted.bam.

Here we see an additional point: the path contains metadata. In this case, the metadata is a deterministic hash of each job. Therefore, identical jobs (those with the same rule metadata and wildcard values) will produce the same hash. By including the hash in the path, Everclear can distinguish two runs of the same job (for example, samtools-sort with the wildcards genome "hg19" and sample "A") that has a different rule definition (for example, different shell commands). However, in Snakemake, each invocation of the same rule, differently defined, would overwrite each other's output.

Of course, in some cases, the user might want to define the exact location of an output path. Therefore, Everclear includes a publish directive that looks much like regular Snakemake paths, like
"{{genome}}/{{sample}}/sorted.bam".

Everclear still uses the regular hash-including path in the workflow but copies the produced files to the publish location afterward.

What jobs the workflow will run can be derived from the workflow definition and sample sheet alone.

Since a sample sheet contains the wildcards, the jobs to produce from each rule are known ahead of time. This knowledge has advantages similar to knowing the DAG from the workflow definition alone.

Additionally, it makes it easier to request files (produce targets) from the command line. Instead of writing out the full path of every file you wish to produce, the user can list the rules for producing the desired files. If the user only desires to run some jobs produced by a rule, the user can restrict which by additional wildcard flags.

So, in Snakemake, creating every file from a target looks like

$ snakemake -f sorted/hg19/A/samtools_sort.bam sorted/hg19/B/samtools_sort.bam

while in Everclear, it looks like

$ everclear samtools-sort

The advantages do not look significant for this simple example, but rules with more wildcards have longer path names, which are harder to remember and type out. Moreover, the user must type out every target desired. In addition, if the user gives an incorrect path, like sorted/A/hg19/samtools_sort.bam, Snakemake will either try to produce a file where the genome wildcard is A and the sample wildcard is hg19 or fail with an error message that no rule could create the desired output.

Advantages of the Everclear software
No wait time to start a workflow

With Snakemake, the user must start the Snakemake software every time they want to execute a command. Then Snakemake must try to find a rule that can produce the desired target, walk backward until it sees whether the constraints are satisfiable, investigate the state of the hard drive to see which files exist, and then run the pipeline.

Since Everclear is continually running and always knows the DAG, it dispatches the jobs needed to produce a target without delay.

Mention inventing on principle by Bret Victor when discussing the server model. Link: https://jamesclear.com/great-speeches/inventing-on-principle-by-bret-victor

Everclear continually evaluates the workflow definition during development and gives immediate feedback.

In Snakemake, to know whether the workflow is valid (for example, cycle-free, no syntax errors, the requirements to produce any target is satisfiable), the user must run Snakemake in --dry-run mode to get feedback. So Snakemake must be started anew for each validity check. Starting the software is time-consuming for the reasons outlined in the previous section. Moreover, to check the whole pipeline, the user must enumerate every target possible to produce in the --dry-run invocation, which is onerous. Furthermore, as explained previously, Snakemake cannot find every DAG inconsistency, including some hard-to-debug errors.

All this is in contrast to Everclear, which quickly, robustly, and continually checks the integrity of the whole DAG.

Everclear makes writing real-time, user-defined dashboards easy.

Everclear keeps all data in atoms, like a database. Whatever is requested is sent to the front-end and kept in a re-frame DB. The Everclear front-end displays these data in a Single Page Application (SPA). Users can add additional panels to this application to produce views and dashboards with data that gives real-time information about the results and the state of the workflow run.

Everclear logging capabilities create logs that allow users to know exactly how the workflow produced a file or arrived at a result.

Workflow management systems know everything about the code run and in what order. It also knows the files it produces. This complete knowledge allows the WMS, in theory, to write comprehensive logs that allow the user, reviewers, and those reading the supplementaries of a paper to understand precisely how the scientists arrived at a result. However, as far as we can see, this complete knowledge is not utilized in any current WMSes.

This automatic log writing dramatically lessens the burden of those writing a paper as they can use it as-is as supplementary material. It is also a boon for reviewers and those wanting to understand the pipeline because it gives them a precise understanding that no journal paper text could give. Human descriptions of a method, no matter how precisely crafted, are necessarily ambiguous. It also helps the workflow creator and project collaborators sanity check the results. Since Everclear logs pertinent details of every file produced, readers can comfortably read through an ordered log and easily spot unexpected results or states of the data without sitting in a terminal and investigating files.

Since Everclear stores the logs as computer-readable data, multiple log views can be produced, each with the details relevant to the reader. The logs can, for example, be interactive web pages that are searchable. Furthermore, while the amount of log data available is massive, one can restrict the view of the log to certain properties, specific wildcards, and subsets of the DAG. By choosing to view all data, but only for a specific set of wildcards, the user gets a complete view of the workflow run, but only for one strand of the execution. For example, in an analysis of the Epigenome Roadmap data, one can restrict the wildcards to a specific tissue type. While the data output for each strand will be different, the flow of execution and code run will be the same. This log would give a complete understanding of the logic of the pipeline.

Since each job has a unique hash, one would never get an inconsistent log with details from different pipeline states. For example, in Snakemake, since each log has a fixed path, the log file output may be stale or from different pipeline versions.

In Everclear, this will not happen: each job has a unique hash, and each job links to its parents by their hash. In this way, Everclear can store logs from all different pipeline versions. This complete log obviates the need to track commit hashes - an error-prone method - to know what pipeline version produced the results.

Since Everclear stores the code and job info, the log files themselves are enough to rerun that specific version of the pipeline.

With other workflow management systems, one can technically read the open-source code to understand the steps performed. However, this is not failsafe for a reason outlined above: one has no guarantee that the code published is the one claimed to produce the results. Furthermore, there are considerable difficulties involved in understanding the workflow; without knowing the final targets, one would not know where to start to work out the rules involved. Also, the rule code is a string template instantiated differently for each job, which adds a layer of complexity to understanding the code.

Everclear's logging system addresses all these problems. Since the log contains the executed code, Everclear knows precisely the code executed for each job. Moreover, since the log also contains the hashes of the parents, one can reproduce the exact sequence of steps used to produce any result. Lastly, since the log contains pertinent info about each file produced, with previews (such as the head and tail), the workflow becomes much easier to understand. This quote from Fred Brooks underlines that seeing the data makes the code and its intent straightforward to understand; and perhaps even obviates the need:

'Show me your flowcharts [code], and conceal your tables [data], and I shall continue to be mystified; show me your tables [data], and I won't usually need your flowcharts [code]: they'll be obvious.'