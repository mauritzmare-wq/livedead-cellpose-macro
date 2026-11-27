# Live/dead image analysis macro (Fiji)

`livedead_cellpose.ijm` batch-analyses confocal calcein / EthD-1 images of A375 cells in GelMA constructs.

The macro prepares the images, segments the cells with the published Cellpose-SAM model, and measures the masks in Fiji. It adds no new segmentation method.

**Requirements:** Fiji and Cellpose-SAM (cpsam), run from the command line on a GPU. Versions: Fiji (ImageJ 1.54p, Java 21.0.7 64-bit), Cellpose 4.2.1.1.

**Folder layout:** `ROOT/TYPE A/DAY 1/Set 1/Snap-7204/Snap-7204.tif` (one calcein, one EthD and one brightfield image per set). Images need the burned-in 200 µm scale bar, which the macro uses to set the scale.

**Use:** open the macro in Fiji, run it, and choose the ROOT folder. Results are written to each set's `Analysis` folder and joined into `ALL_*.csv` files in ROOT. Parameters are listed at the top of the macro (Cellpose-SAM flow threshold 0.3, cellprob threshold 1.0, minimum size 100 px, 1 px mask shrink applied in Fiji). Cellpose-SAM takes no cell diameter; the `cyto3` and 24 px entries in the macro only name folders and the printed Cellpose 3 command.

**Outputs:** per-cell measurements (`ALL_cells_measurements.csv`), per-image summaries with a Note on each processing step (`ALL_per_image_summary.csv`), and per-group calcein and EthD measurements (`ALL_groups_*.csv`). Check the `*_check.jpg` images before trusting results.

Cite Cellpose (Stringer et al. 2021, Nature Methods 18:100-106) and Fiji (Schindelin et al. 2012, Nature Methods 9:676-682) if you use it.
