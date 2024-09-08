# Measure cytoplasmic intensity of exogenous plasmid rRNA

The pipeline uses cellpose to find and segment cytoplasm and nuclei in order to
measure the cytoplasmic intensity of endogenous and exogenous rRNA.  First,
cellpose is run on the endogenous 28S intensity with an expected diameter of
300 pixels on the inverted image with the nuclei model.  Next, cells are detected
with the cyto2 model on the non-inverted endogenous 28S image and an expected
diameter of 500 pixels.  At this point, there are masks for nuclei and cells.
A series of masking operations ensures that each cell has a nucleus, each nucleus
is within a cell, and the cytoplasm is estimated from the cell without the
nucleus.  The endogenous and exogenous intensities are measured in each region
and all objects are related to the cell object.

The `plasmid_cyto.ipynb` notebook shows how to read in the csvs and plot plasmid
intensity within different regions.
