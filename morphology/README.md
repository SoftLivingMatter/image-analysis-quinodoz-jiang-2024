# Measure morphological features

This pipeline measures several quantitative features of nucleoli where there is
no marker for nuclei.  Instead, nearby nucleoli are grouped together and assumed
to be part of one cell.  First, nucleoli are segmented from NPM1 images using
a global, minimum cross-entropy threshold for objects with diameters between
30 and 400 pixels in diameter.  Nucleoli within 100 pixels are grouped as cells.
To sample some of the nucleoplasm for later measures, the merged nucleoli are
dilated with a 50 pixel disk.

To find FCs, the RPA194 image is masked with the dilated GC object and then
enhanced for speckles with 20 pixel feature sizes.  The FCs are segmented from
the enhanced image with an adaptive, three class Otsu threshold where the center
class is set to background, to remove out of plane FCs.  FCs are retained if they
have diameters between 7 and 50 pixels.  FCs which fall inside of the GC are
found by masking the FCs with the GC mask and keeping those that have 50% overlap
with the GC.

DFC objects are found by masking the NOP56 image with the dilated GC objects and
thresholding with and adaptive, three-class Otsu where the middle class is set
to background.  Preliminary DFC objects are filtered to have an area greater than
20 pixels.

To measure the intensity distribution of NOP56 and RPA194, the initial GC objects
were converted to a binary image to normalize the area of each bin.  This was
used during rim enrichment calculation to combine variable numbers of bins.  In
the module, the GC was segmented into 20 bins.

For colocalization measures, the FC, DFC, and GC objects were merged to act as
the objects to measure correlation and overlap.  This limits the effect of background
pixels.  The size, shape and intensity of several objects were measured and
then related to the DilatedGC to act as the primary parent object.

The `morphology.ipynb` notebook shows how to import the measurements and how
to facet several measurements for one example experiment.
