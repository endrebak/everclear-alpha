import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pysam import VariantFile

quals = [record.qual for record in VariantFile("{{input.0}}")]
plt.hist(quals)

with open("{{out-vec.data.0}}", "w+") as h:
    h.write("hi!")

plt.savefig("{{out-vec.plot.0}}")
