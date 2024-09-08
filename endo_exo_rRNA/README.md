# Measure endogenous and exogenous rRNA

Segment GC from NPM1 images and measure the endogenous and exogenous rRNA
intensity, colocalization, and NOP56 distribution.

Image intensity is background corrected by subtracting the 5 percentile intensity
of each channel, for each image.  GC objects are segmented using a global,
minimum cross-entropy threshold and retaining objects between 50 and 300
pixels in diameter.  The size and intensities of the GC are measured.  For
intensity distribution, the GC object is converted to a binary image to allow
for post-hoc combination of adjacent distribution bins.  20 bins are measured
in cellprofiler and the outer 4 are used to estimate the rim enrichment.
Prior to measuring colocalization, the GC object is dialated with a 10 pixel
disk.

The `endo_exo_rRNA.ipynb` notebook shows how to parse the information from
the resulting csvs into a combined pandas dataframe that can be analyzed and
further plotted.
