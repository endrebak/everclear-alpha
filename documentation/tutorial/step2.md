Step 2: Generalizing the read mapping rule.

The first rule we wrote will only work to produce one single file, the alignment of "data/samples/A.fastq". If we want to make the rule take any sample, not just A, we need to provide a wildcard for sample. We do so by adding a wildcards directive to the rule:

:wildcards [:sample]

This directive tells Everclear that it should spawn jobs for every sample wildcard value. For example, A and B. How does Everclear know what wildcards to accept? We tell it so in a wildcards file, a structured file where the columns represent one wildcard and the rows are the different wildcard values.

In our case, we will add this simple wildcards file:

sample
A
B

We must also update our config to point to the file:

:wildcards-file "examples/tutorial/step2/wildcards.tsv"

There is one step left remaining, namely updating our externals file so that each sample wildcard value points to the correct input file:

{:genome {{} "/Users/endrebakkenstovner/everclear/snakemake-flow/data/genome.fa"}
:fastq {{:sample "A"} "/Users/endrebakkenstovner/everclear/snakemake-flow/data/samples/A.fastq"
{:sample "B"} "/Users/endrebakkenstovner/everclear/snakemake-flow/data/samples/B.fastq"}}

Here each map entry for fastq contains the sample wildcard value. The genome entry contains no wildcards {} since we want to align each sample to the same genome.

We now see two jobs in the Everclear dashboard, ":bwa-map {:sample A}" and ":bwa-map {:sample B}". We can see the same by running:

$ everclear.clj --dry-run                                                                                                  
:bwa-map
{:sample "A"} {:sample "B"}

The line :bwa-map represents the rule, and each map in the line immediately below represents a job.

Now we can run every job by

everclear --run

However, instead, let us select only a single job. We do so by restricting the jobs to run by their wildcards.

everclear.clj --dry-run --wildcards {:sample "B"}                                                     
:bwa-map
{:sample "B"}