# Rim-based RDF of different processing factors

Measure the distribution of processing factors around de novo nucleoli.  Since
some nucleoli do not have NPM1 intensity, the NOP56 is used to segment nucleoli.
This requires some modification as the NOP56 intensity is more rough than NPM1.
The workflow starts by background subtracting the 5 percentile intensity value
of each channel, for each image.  The NOP56 image is filtered with a 5 pixel
Gaussian blur, from which nucleoli are segmented with a global, 2-class Otsu
threshold for objects between 10 and 150 pixels in diameter.  The nucleoli
are dilated with a 10 pixel disk to act as the support for RDF measurements.
To measure the total intensity directly outside the nucleolus, the dilated objects
are masked by the original nucleoli to produce a rim region.  The intensity of
each channel and sizes are measured in the nucleoli and rim.  A rim-based
RDF measurement is taken using the nucleoli as sources and dilated nucleoli as the
support with a 10 pixel radius.

The `rim_rdf.ipynb` parses the resulting csv to produce data with the intensities
and a long-form RDF dataframe.  Subsequent cells show how to normalize, average,
and estimate variance of the intensities as a function of distance from the object
boundary.
