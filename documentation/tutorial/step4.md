Step 4: Indexing read alignments

We now need to index the read alignments. Then index allows us to access reads by their mapping location quickly. We will create the index with the samtools index command.

This lesson will not teach us much new syntax but rather more about the workings of Everclear and the publish directive. (Should this lesson be an advanced tutorial?).

(defrule samtools-index
"Index read alignments for random access."
{:wildcards [:sample]
:input     "bam/sorted.bam"
:output    "bam/sorted.bam.bai"
:publish   "{{sample}}/sorted.bam.bai"
:shell     "samtools index {{input.0}}"})

Here we refer to the sorted bam file we created in the previous step. Remember that we used a publish path to specify where to put the output exactly. So which path does Everclear consider the input path, the programmatically created hashed path (<base-path>/c55f6a006f/samtools-sort/A/bam/sorted.bam), or the publish-path (<base-path>/A/sorted.bam). If you inspect the jobs in the browser Everclear dashboard, you will see that Everclear uses the publish path. So the command is:

samtools index /Users/endrebakkenstovner/everclear/tutorial/step4/results/A/sorted.bam

Why does Everclear choose to use the publish path? Samtools index <file.bam> creates an index at <file.bam>.bai. So if Everclear used the hashed-path, like                                    "<base-path>/c55f6a006f/samtools-sort/A/bam/sorted.bam", the index would have ended up at "<base-path>/c55f6a006f/samtools-sort/A/bam/sorted.bam.bai". However, we would want to put it somewhere like "<base-path>/<new-hash>/samtools-index/A/bam/sorted.bam.bai" where <new-hash> is the hash of the samtools-index job. However, software that uses indexed bam files expect the indexes to be located at the exact same place, with only the .bai extension differing. Therefore, we must use the publish directive, even if we do not care about the positioning of the sorted bam files and their indexes.