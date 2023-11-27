# Eclipsebio Customer Toolkit

Thank you for purchasing one of Eclipsebio's products! This repository contains code that can be used to interact with or further process data from Eclipsebio. The following sections list the provided tools and how to use them.

If you need assistance, you can contact our tech support at <techsupport@eclipsebio.com>.

## 1) eSHAPE
[eSHAPE](https://eclipsebio.com/technology/eshape/) is our technology for determining the secondary structure of an RNA. 

### 1.1) extract_shape_reactivities.py
This script is used to extract the per-base reactivities for a given region of interest from an eSHAPE dataset. It requires that Python with pandas and NumPy installed.

The code can be used with:

```
python extract_shape_reactivities.py --bed [BED] --bedgraph [BEDGRAPH] --fasta [FASTA] --scaling [SCALING] --output_folder [OUTPUT]
```

Where the following are provided as inputs
* BED: is the region of interest that you want to obtain reactivity information from in a [BED](https://genome.ucsc.edu/FAQ/FAQformat.html#format1) format. A BED-6 format with strand information is required.
* BEDGRAPH: reactivity information provided by Eclipsebio in a [bedGraph](https://genome.ucsc.edu/goldenPath/help/bedgraph.html) format
* FASTA (optional): if provided, the nucleotide sequence in the region of interest will also be generated
* SCALING (optional): method for normalizing the data. If not provided, raw reactivities are generated. Normalization can be performed using IQR ("IQR") or min-max scaling ("Min-Max")
* OUTPUT (optional): output folder for generated files, defaults to the current directory

The output files end with .shape, .fa, and .map and these files can be used as inputs to a variety of RNA folding algorithms.
