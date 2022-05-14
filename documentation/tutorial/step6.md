In the final step of the beginners tutorial, we will learn how to use a custom Python script instead of the shell commands we have used so far. Then, we will use it to create a plot of our data, a histogram of the variant call quality scores in "all.vcf". Furthermore, we will showcase the expand keyword and how to use it to output multiple files.

(defrule plot-quals
{:input     "all.vcf"
:wildcards [:genome]
:output    {:plot "quals.svg" :data "quals.tsv"}
:script    "/Users/endrebakkenstovner/code/everclar_mirror/examples/snakemake/plot_quals.py"})