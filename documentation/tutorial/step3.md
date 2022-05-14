We need to sort the bam files to index them later. We do so with the command samtools-sort.

Therefore, we need a new rule: samtools-sort

(defrule samtools-sort
"Sort the bams."
{:wildcards [:sample]
:input     "bwa-map.bam"
:output    ["bam/sorted.bam"]
:publish   "{{sample}}/sorted.bam"
:shell     "samtools sort -T {{wildcards.sample}} -O bam {{input.0}} > {{publish.0}}"})

There are two new entries here, namely input and publish. We will first look at the input directive. Notice that it is "bwa-map.bam", the output of the rule we wrote in the last step. Everclear connects these inputs and outputs to produce a graph where the first step is bwa-map and the second is samtools-sort.

Since the rule wildcard is sample, and our wildcards file contains the entries A and B, the rule will also spawn two jobs, A and B. Furthermore, Everclear will connect the job bwa-map {:sample A} to samtools-sort {:sample A} (and similarly for {:sample B}).

The output is a vector ([]) with the file "bam/sorted.bam". We do not need to write the square brackets if we only have one output file; in that case, we could rather use "bam/sorted.bam". We added the vector notation to explain why we had to write {{output.0}} in the previous rule. The index zero denoted that we wanted the first entry of the output vector.

However, note in the shell directive that we use publish, not output, in this rule. The publish-directive allows the user to place the output files in a specific place on the hard drive. For example, in the publish-path, "{{sample}}/sorted.bam", Everclear fills in each double square bracket with the wildcard values for that job. Since the publish-path is relative, Everclear uses the absolute base-path entry in the config file as the parent.

Remember that Everclear creates a path for each output file, which is not user-configurable. This programmatically created output path allows Everclear to run multiple versions of the same job, for example :bwa-map {:genome "hg38} and keep the results separate.

Everclear will also programmatically create an output path even though the rule contains a publish directive. However, after outputting the files, Everclear will copy them to the publish path. We will see later why this is needed.

Finally, note that we need to use {{publish.0}} in the shell-directive to refer to the first (and only) entry of publish.

Now run these jobs.