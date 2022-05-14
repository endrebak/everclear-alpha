This step will aggregate the mapped reads from all the samples and call genomic variants on them. We will use both the samtools and bcftools utilities.

(defrule bcftools-call
"Aggregate mapped reads from all samples and jointly call genomic variants."
{:input     {:sorted "bam/sorted.bam" :index "bam/sorted.bam.bai"}
:output    "all.vcf"
:wildcards []
:externals  [:genome]
:shell     "samtools mpileup -g -f {{externals.genome}} {{input.sorted}} | bcftools call -mv - > {{output.0}}"})

Above, we see three new things: we have multiple input files, use a map to give each input file an alias, and remove the wildcards.

To begin with: having no wildcards means that the rule should only spawn one job. Furthermore, using a map to alias the input files means that we can refer to the input files with input.sorted and input.index (note that we could have used a vector instead, like ["bam/sorted.bam"  "bam/sorted.bam.bai"] and referred to each with input.0 and input.1, but that is less descriptive and more error-prone). The tricky thing to understand here is: what does Everclear do when the rules producing the input (samtools-sorted and samtools-index) have wildcards different from the current rule?

The interpretation is straightforward. Bcftools-call has no wildcards, but samtools-sort and index have one, making bcftools-call the more general rule. So every job produced by samtools-sort and index will be an input to bcftools-call. We confirm this by inspecting the first part of the shell directive in the browser. The template command

samtools mpileup -g -f {{externals.genome}} {{input.sorted}}

becomes

samtools mpileup -g -f /Users/endrebakkenstovner/everclear/snakemake-flow/data/genome.fa /Users/endrebakkenstovner/everclear/tutorial/step5/results/A/sorted.bam /Users/endrebakkenstovner/everclear/tutorial/step5/results/B/sorted.bam

That is, {{input.sorted}} contains two files.

We implicitly defined the rule to accept multiple inputs by using fewer wildcards in the consumer/child (bcftools-call) than the producers/parents (samtools-sort and samtools-index). This technique is powerful as it means we do not have to define the relationship between inputs and outputs with code (like the expand function in Snakemake). Instead, the configuration data that describes the workflow, namely the wildcards table and rule definitions, determine the inputs. The meaning and advantages of this might be a bit abstract now, but we will learn how this can be used to good effect later.