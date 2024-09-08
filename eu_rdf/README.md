# Measure RDF of EU images

This pipeline shows how to use the MeasureRdf module to measure point-based
radial distribution functions.  At a minimum, the point source objects and
an object containing them are required.  Here FCs are the point sources and
are contained in GC objects.  Measurements are stored in the GC objects,
which can be related to the cells which contain them to produce averages
per cell.

First, GC objects are found from the NPM1 image using a global, minimum
cross entropy threshold, retaining objects with pixel diameters between
30 and 400 pixels.  GC objects within 100 pixels are merged and considered
to be from the same cell.  The GC object is used to mask the RPA194 signal
which is further enhanced for speckles with a feature size of 20 pixels.
FC objects are found from the enhanced image with an adaptive, 3-class Otsu
threshold where the middle class is assigned to the background.  This selects
primarily FCs which are in plane with the image.  Adding more FCs will improve
deconvolution, but also dilute the most relevant signal and should be assessed
on a per experiment basis.

The `rdf.ipynb` has code for parsing the csvs, averaging rdfs per cell,
per image and producing plots and color bars.
