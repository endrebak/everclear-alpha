# Expand

For expand do not have the expansion as a wildcard, but rather as an expand entry:

{
:output {:chroms "chromosome.txt"}
:wildcards [:genome]
:expand [:chromosome]
}

How does the child rule know which one of the output files from the parent rule to use?
Must look up the corresponding file in the in-map based on the child wildcards.
