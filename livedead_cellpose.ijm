// ==================================================================================
// Live/dead (calcein + EthD) batch analysis - SET based, safe to run overnight
//
// Folder layout (names must start with DAY and SET; the level between them is optional):
//   ROOT \ TYPE A \ DAY 1 \ [Scaffold 1 \] Set 1 \ Snap-8554 \ Snap-8554.tif
//                                                  Snap-8555 \ Snap-8555.tif
//                                                  Snap-8556 \ Snap-8556.tif
// Choose ROOT when asked (e.g. D:\Processed\3d Culture). Only files named Snap-<number>.tif are used.
//
// For every Set folder the macro:
//   1. classifies each tif as calcein / EthD / brightfield from its colour (filenames are not used),
//   2. records whether the set is complete (one of each) in set_report.csv,
//   3. measures the scale from the burned-in scale bar of each image,
//   4. EthD: removes hot pixels, subtracts calcein bleed-through using THAT set's own green image
//      (fit red = a + k*green on this field only, ignoring outliers),
//   5. rolling-ball background subtraction, Cellpose, per-cell measurements,
//   6. writes everything into  <Set folder>\Analysis\  (CSVs, label images, ROI zips, _check.jpg),
//   7. at the very end builds ALL_*.csv in ROOT by joining all the per-set files.
// If the run is interrupted, start it again: a set is skipped only if Analysis\set_done.txt says it was
// finished with the SAME number of images. To force a full redo set skipDoneSets = false.
// ==================================================================================

// ---------- SETTINGS ----------
var skipDoneSets = true;     // false = redo every set

// Scale
var autoScale = true;        // measure the burned-in scale bar on every image
var scaleBarUm = 200;        // length of that bar in um
var pxPerUm = 1.465;         // fallback if no bar is found (test image: 293 px = 200 um)
var reportInPixels = false;  // true = report in pixels (then set clusterEqDiam / minArea in px)
var maskW = 0.20;            // bottom-right rectangle (fraction of width/height) painted over to hide the bar and label
var maskH = 0.07;

// Channel identification
var colourRatio = 1.5;       // green (red) must be this many times brighter than the other channels
var blankThr = 0.3;          // all channel means below this = completely dark image, assumed EthD

// Hot pixels (single very bright pixels). Removes pixels this much brighter than their 5x5 median. 0 = off.
var hotPixelsEthD = 50;
var hotPixelsCalcein = 0;

// Calcein bleed-through into the red image, estimated per set from its own green image
var doBleedCorrection = true;
var minBleedCorr = 0.3;      // if red and green correlate less than this the correction is skipped (and noted)
var maxBleedK = 0.10;        // if the fitted leak (red per unit green) is above this, the red image is mostly calcein:
                             // EthD is then NOT analysed and recorded as NA (not as 0 dead cells). Typical good sets: about 0.02

// Where Cellpose runs. We used "external" (command line, GPU).
//   "external": Cellpose runs OUTSIDE Fiji, e.g. from the command line with your GPU. You run this macro twice:
//               step 1 prepares the images and saves them in <ROOT>\ToSegment, you run one Cellpose command (the macro prints it),
//               step 2 reads Cellpose's masks back in and does all the measuring. The macro asks which step when it starts.
//   "fiji":     Fiji's own Cellpose command (CPU unless a GPU environment is installed in Fiji). One run does everything.
var cellposeMode = "external";
var stage = 0;               // set by the question at the start (external mode only)
var cellposeModel = "cyto3"; // only used in the printed Cellpose 3 command

// Cellpose settings. The data in the article were segmented with Cellpose-SAM (cpsam), which takes no cell
// diameter and no model name: cellposeModel and cellDiameter only matter for the Cellpose 3 command the macro prints.
var cellposeCmd   = "Cellpose...";
var cellDiameter  = 24;      // px; with Cellpose-SAM it only names the d24 folder
var cellprobThr   = 1.0;     // higher = tighter masks (0 = Cellpose default)
var flowThr       = 0.3;     // lower = rejects odd shapes (0.4 = Cellpose default)
var minSize       = 100;     // px area; smallest object kept (also passed to Cellpose in external mode)
var maskShrinkPx  = 1;       // px each cell mask is shrunk by before measuring (removes the blurry rim; 0 = off). The saved label image is not changed.
var cellposeReturnRois = "false";

// EthD stains nuclei, so dead cells appear as small, compact bright spots on a noisy background. Cellpose is not suited to
// that (it labels the noise), so EthD uses a spot detector instead: a band-pass filter tuned to the spot size, a threshold set
// from the image's own robust noise level, and a size range. Set ethdMethod = "cellpose" to use Cellpose for EthD as before.
var ethdMethod = "spots";
var spotSigma = 2.0;         // px; fine scale of the band-pass, about 1/3 of a nucleus diameter (nucleus ~7-15 px)
var spotK = 5;               // spot must stand this many robust noise-SDs above the background
var spotMinContrast = 4;     // and be at least this much brighter than its surroundings (units of the subtracted image)
var spotMinArea = 10;        // px^2, smaller is ignored
var spotMaxArea = 600;       // px^2, larger is ignored (not a nucleus)

// Background subtraction for CALCEIN: "masked" = object-aware (recommended: keeps dim cores inside bright spheroids),
// "rolling" = rolling ball with the per-day radius below (can wipe out dim cores).
var bgMethodCalcein = "masked";
var bgMargin = 6;            // objects = smoothed image brighter than its median by this much (raw units)
var bgMaskGrow = 8;          // px the object mask is grown so blurred edges are not used as background
var bgSigma = 40;            // px, smoothness of the background estimate

// Per-day background radius and Cellpose diameter (px). The values below were set by eye for our images; re-measure them for other data.
// Rule: rolling-ball radius >= radius of the largest object that day (about 0.75 x its diameter, never below 25).
var bgSize        = 25;      // used for days not in the table
var dayNames      = newArray("DAY 1", "DAY 3", "DAY 6", "DAY 9", "DAY 12", "DAY 15");
var dayDiameter   = newArray(24,       24,      24,      24,      24,       24);   // single-cell diameter on every day: clusters are handled by the grouping below
var dayBgSize     = newArray(25,       25,      40,      80,      150,      220);

// SPHEROID GROUPING (calcein images). The cell-number and diameter thresholds below were chosen by looking at our images;
// they are not taken from the literature. Cellpose finds single cells everywhere; cells that touch or lie within linkDistUm of
// each other are then merged into one group. A group of at least minCellsSpheroid cells = "spheroid", 2 to
// minCellsSpheroid-1 cells = "small cluster", 1 cell = "single cell". Spheroid-level area, diameter and calcein
// intensity go to groups_measurements.csv, and the masks are saved as *_spheroid_mask.tif and *_singlecells_mask.tif.
// 1 = group on that day (same order as dayNames above). Days 1 and 3 are single cells only.
var dayGroup = newArray(0, 0, 1, 1, 1, 1);
// How the group (spheroid) outline is made:
//   "blob":  the whole bright blob in the smoothed calcein image (so a dimmer core inside the same blob is included),
//            kept only if it contains in-focus cells.
//   "cells": the joined-up masks of the individual cells only (a dim core with no detected cells can be left out).
var groupMethod = "blob";
var blobThr = 10;            // blob = smoothed calcein (after background subtraction) above this value; lower it if dim cores are cut out
var blobSigma = 2;           // px, smoothing before thresholding
var blobRel = 0.35;          // edge of a spheroid = where the smoothed calcein falls to this fraction of the brightest value within blobMaxR px
                             // (a half-maximum style edge: it keeps the glow around very bright spheroids out of the outline and separates
                             // neighbouring spheroids). 0 = off (fixed threshold blobThr only). Lower = looser outline, higher = tighter.
var blobMaxR = 30;           // px, how far to look for the local brightest value
var blobGrow = 0;            // px added around the blob outline
var linkDistUm = 5;          // ("cells" method) gap (um) up to which neighbouring cells are joined into one group
var groupGrow = 2;           // px added around the joined cells so the group mask reaches the real edge
var minCellsSpheroid = 10;       // a group needs at least this many cells ...
var ringPx = 15;             // px width of the ring around each group used as local EthD background
var spheroidMinMean = 30;        // a group called a spheroid only by its size must also be at least this bright (calcein mean after background subtraction), so large dim out-of-focus discs are not spheroids
var satLevel = 250;               // raw calcein value counted as saturated (8-bit max is 255); SatFraction in groups_measurements.csv = share of a group at or above it
var minSpheroidDiamUm = 100;     // ... or at least this diameter (for spheroids that Cellpose returns as one big object)
var ethdUseSingleCellSettings = true;   // EthD images always use cellDiameter and bgSizeEthD
var bgSizeEthD    = 50;

// "No cells" check: needs at least minPixAbove pixels brighter than minPeak (after subtraction). One hot pixel is not enough.
var minPeakCalcein = 30;
var minPeakEthD    = 8;      // tune with a field that has REAL dead cells; a field with no dead cells should come out as "no signal"
var minPixAbove    = 150;

// Noise floor for SEGMENTATION ONLY (subtracted from a copy before Cellpose; measurements are unaffected)
var segFloorCalcein = 0;
var segFloorEthD    = 0;

// FOCUS FILTER (optional). Out-of-focus cells are dim, smooth discs; in-focus cells have sharp edges.
// The macro finds bright objects, scores each by how sharp its edge is (band-pass energy divided by local brightness),
// and blacks out objects scoring below focusThr in the copy Cellpose sees. The measurements still use the full image.
// Check the *_focus_check.jpg (green = kept, red = removed) on several images before trusting it.
// Large objects (spheroids) can be scored too low: leave it off for days dominated by spheroids, or lower focusThr.
var focusFilter = true;
var focusThr    = 0.10;      // objects with a sharpness score below this are removed (higher = stricter)
var focusS1 = 1.0;           // band-pass: fine Gaussian (px)
var focusS2 = 2.5;           // band-pass: coarse Gaussian (px)
var focusWin = 10;           // averaging scale of the sharpness map (px), about one cell radius
var focusMinSig = 3;         // added to local brightness so dark noise is not amplified
var focusObjThr = 4;         // objects = areas brighter than this (after background subtraction and smoothing)
var focusKeepDiamUm = 80;    // objects at least this wide (um) are never removed as out of focus (large spheroids have smooth edges and would score low). 0 = off
var focusMinObj = 60;        // ignore objects smaller than this (px^2)
var focusVerbose = false;    // true = print every object's score to the Log (for choosing focusThr)

// Cells are never deleted; the CSV has Keep = 1/0 so you can filter later. 0 = filter off.
var minCellMeanCalcein = 0;
var minCellMeanEthD    = 0;
var minArea = 0;
var maxAspectEthD = 0;       // >0: EthD objects with Feret/MinFeret above this get Keep = 0 (e.g. 2.5)
var clusterEqDiam = 30;      // equivalent diameter (um) below this = single cell, above = cluster/spheroid
// ------------------------------

// Only rootDir and registry are global, and both are set once here and then only read.
// Everything else is passed between functions as arguments or return values.
var rootDir, registry;

rootDir = getDirectory("Choose ROOT folder (e.g. D:\\Processed\\3d Culture)");
registry = rootDir + "sets_registry.txt";
if (cellposeMode == "external") stage = chooseStage();
if (!(cellposeMode == "external" && stage == 1) && File.exists(registry)) ok = File.delete(registry);
print("Run started: " + rootDir);
setBatchMode(true);
setOption("ScaleConversions", false);
setOption("BlackBackground", true);
run("Set Measurements...", "area mean integrated feret's redirect=None decimal=4");
nLoose = walk(rootDir);
setBatchMode(false);
if (cellposeMode == "external" && stage == 1) {
    showCellposeInstructions();
} else {
    buildMasters();
    if (nLoose > 0) print("NOTE: " + nLoose + " Snap tif(s) were not inside any Set folder and were skipped.");
    nSetsHandled = 0;
    if (File.exists(registry)) {
        regArr = split(File.openAsString(registry), "\n");
        nSetsHandled = regArr.length;
    }
    print("DONE. " + nSetsHandled + " set folders handled. Combined files are in " + rootDir);
}
if (isOpen("Log")) { selectWindow("Log"); saveAs("Text", rootDir + "run_log.txt"); }

// ================================ GPU / external Cellpose helpers ================================
function chooseStage() {
    Dialog.create("Cellpose on the GPU");
    Dialog.addMessage("Step 1: prepare the images for Cellpose.\nStep 2: measure, after you have run Cellpose.");
    Dialog.addChoice("Which step?", newArray("1: prepare images for Cellpose", "2: measure (after Cellpose has run)"));
    Dialog.show();
    c = Dialog.getChoice();
    if (startsWith(c, "1")) return 1;
    return 2;
}

// Prints (and saves) the Cellpose command(s) to run on the GPU after step 1
function showCellposeInstructions() {
    base = rootDir + "ToSegment" + File.separator;
    txt = "";
    if (File.exists(base)) {
        lst = getFileList(base);
        for (i = 0; i < lst.length; i++) {
            if (endsWith(lst[i], "/") && startsWith(lst[i], "d")) {
                dname = substring(lst[i], 0, lengthOf(lst[i]) - 1);
                diam = substring(dname, 1);
                folder = base + dname;
                txt = txt + "Folder " + folder + "  (cell diameter " + diam + " px)\n" +
                      "  Cellpose 3 (cyto3, nuclei or your own model):\n" +
                      "    python -m cellpose --dir \"" + folder + "\" --pretrained_model " + cellposeModel + " --diameter " + diam +
                      " --flow_threshold " + flowThr + " --cellprob_threshold " + cellprobThr + " --min_size " + minSize + " --use_gpu --save_tif --no_npy\n" +
                      "  Cellpose 4 / Cellpose-SAM:\n" +
                      "    python -m cellpose --dir \"" + folder + "\" --pretrained_model cpsam" +
                      " --flow_threshold " + flowThr + " --cellprob_threshold " + cellprobThr + " --min_size " + minSize + " --use_gpu --save_tif --no_npy\n\n";
            }
        }
    }
    if (txt == "") txt = "No images were prepared (nothing found in " + base + ").\n";
    txt = "STEP 1 FINISHED.\nNow open the terminal / Anaconda prompt of the environment where the Cellpose GUI works with your GPU\n" +
          "and run the command for your Cellpose version (if an option is not recognised, delete that option):\n\n" + txt +
          "When Cellpose has finished (each image gets a *_cp_masks.tif next to it), run this macro again and choose step 2.\n";
    print(txt);
    if (File.exists(base)) File.saveString(txt, base + "RUN_CELLPOSE.txt");
}

// ================================ folder walking ================================
// Only files named exactly Snap-<digits>.tif are used (so exports like Snap-7242_c1-2.tif are ignored)
function isSnap(f) {
    return matches(f, "(?i)Snap-[0-9]+\\.tiff?");
}

// Returns the number of Snap tifs found outside any Set folder
function walk(dir) {
    loose = 0;
    list = getFileList(dir);
    for (i = 0; i < list.length; i++) {
        if (endsWith(list[i], "/")) {
            nm = substring(list[i], 0, lengthOf(list[i]) - 1);
            if (startsWith(toUpperCase(nm), "SET")) processSet(dir + list[i]);
            else if (toUpperCase(nm) != "ANALYSIS" && toUpperCase(nm) != "CELLPOSE_OUTPUT") loose = loose + walk(dir + list[i]);
        } else if (isSnap(list[i])) {
            loose = loose + 1;
        }
    }
    return loose;
}

// Returns all Snap tif paths below dir, each followed by "|"
function collectTifs(dir) {
    acc = "";
    list = getFileList(dir);
    for (i = 0; i < list.length; i++) {
        if (endsWith(list[i], "/")) {
            nm = substring(list[i], 0, lengthOf(list[i]) - 1);
            if (toUpperCase(nm) != "ANALYSIS" && toUpperCase(nm) != "CELLPOSE_OUTPUT") acc = acc + collectTifs(dir + list[i]);
        } else if (isSnap(list[i])) {
            acc = acc + dir + list[i] + "|";
        }
    }
    return acc;
}

// ================================ one set folder ================================
function processSet(setDir) {
    fparts = split(replace(setDir, "\\", "/"), "/");
    n = fparts.length;
    sSet = fparts[n - 1];
    sGroup = ""; sDay = ""; sScaffold = "";
    dayIdx = -1;
    for (q = 0; q < n - 1; q++) {
        if (startsWith(toUpperCase(fparts[q]), "DAY")) dayIdx = q;
    }
    if (dayIdx < 0) {
        print("WARNING: Set folder is not inside a DAY folder, skipped: " + setDir);
        return;
    }
    sDay = fparts[dayIdx];
    if (dayIdx > 0) sGroup = fparts[dayIdx - 1];
    if (n - 2 > dayIdx) sScaffold = fparts[n - 2];
    setLabel = sGroup + " | " + sDay + " | " + sScaffold + " | " + sSet;

    sDir = setDir + "Analysis" + File.separator;
    tifStr = collectTifs(setDir);
    nT = 0;
    if (lengthOf(tifStr) > 0) { tifArr = split(tifStr, "|"); nT = tifArr.length; }

    stage1 = (cellposeMode == "external" && stage == 1);

    // skip only sets that were finished with the same, non-zero number of images
    doneF = sDir + "set_done.txt";
    if (skipDoneSets && nT > 0 && File.exists(doneF) && !stage1) {
        doneTxt = replace(replace(File.openAsString(doneF), "\n", ""), "\r", "");
        // a set whose earlier run had an ERROR (for example Cellpose masks that did not exist yet) is never skipped
        prevTxt = "";
        if (File.exists(sDir + "per_image_summary.csv")) prevTxt = File.openAsString(sDir + "per_image_summary.csv");
        if (doneTxt == "done " + nT && indexOf(prevTxt, "ERROR") < 0) {
            print("Already done, skipped: " + setLabel + "  (" + nT + " images)");
            File.append(sDir, registry);
            return;
        }
    }

    print("=== SET: " + setLabel + " ===  found " + nT + " Snap tif(s)");
    if (nT == 0) {
        lst = getFileList(setDir);
        shown = "";
        for (q = 0; q < lst.length && q < 12; q++) shown = shown + lst[q] + "  ";
        print("  nothing found. This folder contains: " + shown);
    }

    File.makeDirectory(sDir);
    if (!stage1) {
    cellsF = sDir + "cells_measurements.csv";
    imgF   = sDir + "per_image_summary.csv";
    setF   = sDir + "set_report.csv";
    groupsF = sDir + "groups_measurements.csv";
    if (File.exists(cellsF)) ok = File.delete(cellsF);
    if (File.exists(imgF))   ok = File.delete(imgF);
    if (File.exists(setF))   ok = File.delete(setF);
    if (File.exists(groupsF)) ok = File.delete(groupsF);
    groupsEF = sDir + "groups_ethd.csv";
    if (File.exists(groupsEF)) ok = File.delete(groupsEF);
    if (File.exists(cellsF) || File.exists(imgF) || File.exists(setF) || File.exists(groupsF) || File.exists(groupsEF)) {
        print("WARNING: result files in " + sDir + " could not be replaced (are they open in Excel or another program?). SET SKIPPED - close them and run step 2 again.");
        return;
    }
    File.append("Group,Day,Scaffold,Set,SetComplete,Channel,File,CellID,Area,Mean,IntDen,Feret,MinFeret,EqDiameter,Keep", cellsF);
    File.append("Group,Day,Scaffold,Set,SetComplete,Channel,File,NumCells,NumKept,NumSingle,NumCluster,KeptArea,KeptIntDen,ImageArea,BleedK,BleedOffset,BleedR,Note", imgF);
    File.append("Group,Day,Scaffold,Set,SetComplete,Channel,File,GroupID,Class,NumCells,Area,EqDiameter,Feret,MinFeret,Mean,IntDen,SatFraction", groupsF);
    File.append("Group,Day,Scaffold,Set,SetComplete,CalceinFile,EthDFile,GroupID,Class,NumCells,Area,CalceinMean,CalceinIntDen,EthDInsideMean,EthDRingMean,EthDNetMean,EthDNetIntDen,EthDNetToCalceinMean", groupsEF);
    File.append("Group,Day,Scaffold,Set,NumTifs,NumCalcein,NumEthD,NumBrightfield,NumUnclassified,SetComplete,Issues", setF);
    }

    chan = newArray(maxOf(nT, 1));
    notes = newArray(maxOf(nT, 1));
    nCal = 0; nEth = 0; nBf = 0; nUnk = 0;
    greenPath = "";
    for (t = 0; t < nT; t++) {
        res = split(classifyFile(tifArr[t]), "|");
        chan[t] = res[0];
        notes[t] = "";
        if (res.length > 1) notes[t] = res[1];
        if (chan[t] == "calcein") { nCal = nCal + 1; greenPath = tifArr[t]; }
        else if (chan[t] == "EthD") nEth = nEth + 1;
        else if (chan[t] == "brightfield") nBf = nBf + 1;
        else nUnk = nUnk + 1;
    }
    issues = "";
    if (nT == 0) issues = issues + "no tif files; ";
    if (nCal == 0) issues = issues + "missing calcein; ";
    if (nCal > 1)  issues = issues + nCal + " calcein images; ";
    if (nEth == 0) issues = issues + "missing EthD; ";
    if (nEth > 1)  issues = issues + nEth + " EthD images; ";
    if (nBf == 0)  issues = issues + "missing brightfield; ";
    if (nBf > 1)   issues = issues + nBf + " brightfield images; ";
    if (nUnk > 0)  issues = issues + nUnk + " unclassified images; ";
    sComplete = "yes";
    if (issues != "") sComplete = "NO";
    if (!stage1) File.append(sGroup + "," + sDay + "," + sScaffold + "," + sSet + "," + nT + "," + nCal + "," + nEth + "," + nBf + "," + nUnk + "," + sComplete + "," + issues, setF);
    if (sComplete == "NO") print("  INCOMPLETE SET: " + issues);

    if (nCal != 1) greenPath = "";   // no unique green image to pair with
    prefix = sGroup + "," + sDay + "," + sScaffold + "," + sSet + "," + sComplete;
    // calcein first: the EthD step measures inside the spheroid masks found on the calcein image
    for (t = 0; t < nT; t++) {
        if (chan[t] == "calcein") analyseImage(tifArr[t], "calcein", "", notes[t], prefix, sDir, sDay);
    }
    if (!stage1) {
        for (t = 0; t < nT; t++) {
            if (chan[t] == "EthD") analyseImage(tifArr[t], "EthD", greenPath, notes[t], prefix, sDir, sDay);
        }
        File.saveString("done " + nT, doneF);
        File.append(sDir, registry);
    }
    run("Close All");
    call("java.lang.System.gc");
}

// ================================ channel classification ================================
// Returns "channel|note"
function classifyFile(path) {
    clsNote = "";
    open(path);
    t = getTitle();
    ch = "unclassified";
    if (bitDepth() == 24) {
        w = getWidth(); h = getHeight();
        run("Split Channels");
        selectWindow(t + " (red)");
        makeRectangle(0, 0, round(w * 0.8), round(h * 0.93));   // leaves out the scale bar corner
        getStatistics(ar, mr);
        selectWindow(t + " (green)");
        makeRectangle(0, 0, round(w * 0.8), round(h * 0.93));
        getStatistics(ag, mg);
        selectWindow(t + " (blue)");
        makeRectangle(0, 0, round(w * 0.8), round(h * 0.93));
        getStatistics(ab, mb);
        mx = maxOf(mr, maxOf(mg, mb));
        if (mx < blankThr) { ch = "EthD"; clsNote = "completely dark image; assumed EthD"; }
        else if (mg > colourRatio * mr && mg > colourRatio * mb) ch = "calcein";
        else if (mr > colourRatio * mg && mr > colourRatio * mb) ch = "EthD";
        else ch = "brightfield";
    } else {
        getLut(lr, lg, lb);
        if (lg[255] > 200 && lr[255] < 100 && lb[255] < 100) ch = "calcein";
        else if (lr[255] > 200 && lg[255] < 100 && lb[255] < 100) ch = "EthD";
    }
    run("Close All");
    return ch + "|" + clsNote;
}

// ================================ scale bar ================================
// Longest run of pure-white pixels in the bottom-right corner = scale bar length in px
function detectBarPx() {
    w = getWidth(); h = getHeight();
    x0 = round(w * 0.75); y0 = round(h * 0.90);
    best = 0;
    for (y = y0; y < h; y += 2) {
        runLen = 0;
        for (x = x0; x < w; x++) {
            v = getPixel(x, y) & 0xffffff;
            if (v == 0xffffff) { runLen++; if (runLen > best) best = runLen; } else runLen = 0;
        }
    }
    return best;
}

// re-applies a pixel size (um per pixel) to the active image
function setPixelSize(pw) {
    if (reportInPixels) run("Set Scale...", "distance=0 known=0 unit=pixel");
    else run("Set Scale...", "distance=1 known=" + d2s(pw, 8) + " unit=um");
}

// ================================ load one channel ================================
// Opens a tif, sets the scale, keeps the wanted colour channel, paints over the scale bar,
// removes hot pixels, converts to 32-bit and renames it outTitle. Returns "imageID|scaleNote".
function loadChannel(path, channel, outTitle) {
    open(path);
    orig = getTitle();
    scaleNote = "";
    if (reportInPixels) {
        run("Set Scale...", "distance=0 known=0 unit=pixel");
        scaleNote = "pixels";
    } else {
        barPx = 0;
        if (autoScale && bitDepth() == 24) barPx = detectBarPx();
        if (barPx >= 50) {
            run("Set Scale...", "distance=" + barPx + " known=" + scaleBarUm + " unit=um");
            scaleNote = "scale bar " + barPx + " px = " + scaleBarUm + " um";
        } else {
            run("Set Scale...", "distance=" + pxPerUm + " known=1 unit=um");
            scaleNote = "NO scale bar found: used " + pxPerUm + " px/um";
            print("  WARNING: no scale bar found in " + orig);
        }
    }

    if (bitDepth() == 24) {
        run("Split Channels");
        rN = orig + " (red)"; gN = orig + " (green)"; bN = orig + " (blue)";
        if (channel == "EthD") keepName = rN; else keepName = gN;
        names = newArray(rN, gN, bN);
        for (c = 0; c < 3; c++) {
            if (names[c] != keepName) { selectWindow(names[c]); close(); }
        }
        selectWindow(keepName);
    }
    rename(outTitle);

    // paint over the scale bar corner with the image's background level
    mw = round(getWidth() * maskW);
    mh = round(getHeight() * maskH);
    if (mw > 0 && mh > 0) {
        bgLevel = getValue("Median");
        setColor(bgLevel);
        fillRect(getWidth() - mw, getHeight() - mh, mw, mh);
        run("Select None");
    }

    hp = hotPixelsCalcein;
    if (channel == "EthD") hp = hotPixelsEthD;
    if (hp > 0) run("Remove Outliers...", "radius=2 threshold=" + hp + " which=Bright");

    run("32-bit");
    return "" + getImageID() + "|" + scaleNote;
}

// ================================ bleed-through fit ================================
// Fits  red = a + k * green  on 16x16-pixel block averages of THIS field, rejecting outliers
// (real EthD signal) so they do not pull k upwards. Rejection is two-sided so k cannot drift downwards.
// Returns "k|a|r"
function fitBleed(gID, rID) {
    selectImage(gID);
    run("Scale...", "x=0.0625 y=0.0625 interpolation=Bilinear average create title=gs");
    gsID = getImageID();
    selectImage(rID);
    run("Scale...", "x=0.0625 y=0.0625 interpolation=Bilinear average create title=rs");
    rsID = getImageID();

    selectImage(gsID);
    w = getWidth(); h = getHeight();
    n = w * h;
    xs = newArray(n); ys = newArray(n); use = newArray(n);
    i = 0;
    for (yy = 0; yy < h; yy++) {
        for (xx = 0; xx < w; xx++) { xs[i] = getPixel(xx, yy); use[i] = 1; i++; }
    }
    selectImage(rsID);
    i = 0;
    for (yy = 0; yy < h; yy++) {
        for (xx = 0; xx < w; xx++) { ys[i] = getPixel(xx, yy); i++; }
    }
    selectImage(gsID); close();
    selectImage(rsID); close();

    k = 0; a = 0; r = 0;
    for (iter = 0; iter < 4; iter++) {
        sx = 0; sy = 0; sxx = 0; sxy = 0; syy = 0; m = 0;
        for (i = 0; i < n; i++) {
            if (use[i] == 1) {
                sx = sx + xs[i]; sy = sy + ys[i];
                sxx = sxx + xs[i] * xs[i]; sxy = sxy + xs[i] * ys[i]; syy = syy + ys[i] * ys[i];
                m++;
            }
        }
        den = m * sxx - sx * sx;
        if (m < 100 || den <= 0) iter = 99;
        else {
            k = (m * sxy - sx * sy) / den;
            a = (sy - k * sx) / m;
            ss = 0;
            for (i = 0; i < n; i++) {
                if (use[i] == 1) { res = ys[i] - (a + k * xs[i]); ss = ss + res * res; }
            }
            sd = sqrt(ss / m);
            if (sd < 0.1) sd = 0.1;
            for (i = 0; i < n; i++) {
                if (use[i] == 1 && abs(ys[i] - (a + k * xs[i])) > 3 * sd) use[i] = 0;
            }
            dy = m * syy - sy * sy;
            if (dy > 0) r = (m * sxy - sx * sy) / sqrt(den * dy);
        }
    }
    return d2s(k, 8) + "|" + d2s(a, 8) + "|" + d2s(r, 8);
}

// ================================ focus filter ================================
// Blacks out out-of-focus objects in the segmentation copy (segID). Returns "nObjects|nRemoved".
// Saves a check picture to checkPath: green outline = kept, red outline = removed.
function removeOutOfFocus(subID, segID, checkPath) {
    // sharpness map = sqrt( local mean of (G1 - G2)^2 ) / (local mean brightness + focusMinSig)
    selectImage(subID);
    run("Duplicate...", "title=fa");
    run("Gaussian Blur...", "sigma=" + focusS1);
    faID = getImageID();
    selectImage(subID);
    run("Duplicate...", "title=fb");
    run("Gaussian Blur...", "sigma=" + focusS2);
    fbID = getImageID();
    imageCalculator("Subtract create 32-bit", faID, fbID);
    enID = getImageID();
    run("Square");
    run("Gaussian Blur...", "sigma=" + focusWin);
    run("Square Root");
    selectImage(faID); close();
    selectImage(fbID); close();
    selectImage(subID);
    run("Duplicate...", "title=lv");
    run("Gaussian Blur...", "sigma=" + focusWin);
    run("Add...", "value=" + focusMinSig);
    lvID = getImageID();
    imageCalculator("Divide create 32-bit", enID, lvID);
    shID = getImageID();
    selectImage(enID); close();
    selectImage(lvID); close();

    // objects = bright areas of the lightly smoothed image
    selectImage(subID);
    run("Duplicate...", "title=ob");
    run("Gaussian Blur...", "sigma=2");
    obID = getImageID();
    setThreshold(focusObjThr, 1000000000);
    roiManager("reset");
    run("Analyze Particles...", "size=" + focusMinObj + "-Infinity add");
    resetThreshold();
    nObj = roiManager("count");

    // check picture
    selectImage(subID);
    run("Duplicate...", "title=fchk");
    run("Enhance Contrast", "saturated=0.35");
    chkID = getImageID();

    nRem = 0;
    for (i = 0; i < nObj; i++) {
        selectImage(shID);
        roiManager("select", i);
        List.setMeasurements;
        objDiam = 2 * sqrt(List.getValue("Area") / PI);
        run("Enlarge...", "enlarge=4");
        getStatistics(ar, mn, mi, mx);
        score = mx;
        run("Select None");
        if (focusVerbose) print("  focus object " + (i + 1) + ": score " + d2s(score, 3));
        selectImage(chkID);
        roiManager("select", i);
        if (score < focusThr && !(focusKeepDiamUm > 0 && objDiam >= focusKeepDiamUm)) {
            Overlay.addSelection("red", 2);
            nRem = nRem + 1;
            selectImage(segID);
            roiManager("select", i);
            run("Enlarge...", "enlarge=3");
            setColor(0);
            fill();
            run("Select None");
        } else {
            Overlay.addSelection("green", 2);
        }
        selectImage(chkID);
        run("Select None");
    }
    selectImage(chkID);
    run("Flatten");
    saveAs("Jpeg", checkPath);
    selectImage(shID); close();
    selectImage(obID); close();
    roiManager("reset");
    return "" + nObj + "|" + nRem;
}

// ================================ group cells into spheroids ================================
// Uses the cell ROIs in the ROI Manager. Returns "nSpheroids|nSmallClusters|nIsolatedCells|spheroidArea".
function groupSpheroids(subID, nCells, prefix, channel, fname, sDir, base, rawID) {
    selectImage(subID);
    getPixelSize(pu, pw, ph);
    w = getWidth(); h = getHeight();
    cx = newArray(nCells); cy = newArray(nCells);

    // binary picture of all cells (also remembers each cell's centre)
    newImage("cm", "8-bit black", w, h, 1);
    cmID = getImageID();
    for (r = 0; r < nCells; r++) {
        selectImage(subID);
        roiManager("select", r);
        List.setMeasurements;
        cx[r] = List.getValue("X") / pw;
        cy[r] = List.getValue("Y") / ph;
        if (groupMethod == "cells") {
            selectImage(cmID);
            roiManager("select", r);
            setColor(255);
            fill();
            run("Select None");
        }
    }

    if (groupMethod == "blob") {
        // the whole bright blob in the smoothed calcein image (includes a dimmer core), holes filled
        selectImage(cmID); close();
        selectImage(subID);
        run("Select None");
        run("Duplicate...", "title=bm");
        run("Gaussian Blur...", "sigma=" + blobSigma);
        bmID = getImageID();
        if (blobRel > 0) {
            // local brightest value (grayscale dilation, smoothed), times blobRel = the edge level for each object
            run("Duplicate...", "title=lm");
            run("Maximum...", "radius=" + blobMaxR);
            run("Gaussian Blur...", "sigma=" + (blobMaxR / 3));
            run("Multiply...", "value=" + blobRel);
            lmID = getImageID();
            imageCalculator("Subtract create 32-bit", bmID, lmID);
            dID = getImageID();
            setThreshold(0, 1000000000);
            run("Convert to Mask");
            m1ID = getImageID();
            selectImage(bmID);
            setThreshold(blobThr, 1000000000);
            run("Convert to Mask");
            imageCalculator("AND create", m1ID, bmID);
            cmID = getImageID();
            selectImage(m1ID); close();
            selectImage(bmID); close();
            selectImage(lmID); close();
            selectImage(cmID);
            run("Fill Holes");
            run("Minimum...", "radius=2");     // opening: removes thin spurs
            run("Maximum...", "radius=2");
        } else {
            setThreshold(blobThr, 1000000000);
            run("Convert to Mask");
            run("Fill Holes");
            cmID = getImageID();
        }
        if (blobGrow > 0) { selectImage(cmID); run("Maximum...", "radius=" + blobGrow); }
    } else {
        // join neighbouring cells: closing with radius = half the link distance, fill holes, then grow a little
        selectImage(cmID);
        rLink = round(linkDistUm / pw / 2);
        if (rLink > 0) {
            run("Maximum...", "radius=" + rLink);
            run("Fill Holes");
            run("Minimum...", "radius=" + rLink);
        } else {
            run("Fill Holes");
        }
        if (groupGrow > 0) run("Maximum...", "radius=" + groupGrow);
    }
    selectImage(cmID);
    setThreshold(255, 255);
    roiManager("reset");
    run("Analyze Particles...", "size=0-Infinity add");
    resetThreshold();
    nG = roiManager("count");

    newImage("sphmask", "8-bit black", w, h, 1);
    sphID = getImageID();
    newImage("isomask", "8-bit black", w, h, 1);
    isoID = getImageID();
    selectImage(subID);
    run("Select None");          // a leftover selection would make the copy only as big as that selection
    run("Duplicate...", "title=gchk");
    run("Enhance Contrast", "saturated=0.35");
    gchkID = getImageID();
    grpF = sDir + "groups_measurements.csv";
    infoF = sDir + base + "_groups_info.csv";
    if (File.exists(infoF)) ok = File.delete(infoF);
    File.append("GroupID,Class,NumCells,Area,Mean,IntDen", infoF);

    nSph = 0; nClu = 0; nIso = 0; sphArea = 0;
    for (g = 0; g < nG; g++) {
        selectImage(subID);
        roiManager("select", g);
        List.setMeasurements;
        area = List.getValue("Area");
        mean = List.getValue("Mean");
        intden = List.getValue("IntDen");
        feretV = List.getValue("Feret");
        minF = List.getValue("MinFeret");
        nIn = 0;
        for (j = 0; j < nCells; j++) {
            if (selectionContains(cx[j], cy[j])) nIn = nIn + 1;
        }
        run("Select None");
        // share of the group's pixels that are saturated in the raw calcein image (saturated intensities are unreliable)
        selectImage(rawID);
        setThreshold(satLevel, 1000000000);
        roiManager("select", g);
        List.setMeasurements("limit");
        satFrac = List.getValue("%Area") / 100;
        resetThreshold();
        run("Select None");
        cls = "single cell";
        colour = "cyan";
        if (nIn >= minCellsSpheroid || (minSpheroidDiamUm > 0 && 2 * sqrt(area / PI) >= minSpheroidDiamUm && mean >= spheroidMinMean)) { cls = "spheroid"; colour = "yellow"; nSph = nSph + 1; sphArea = sphArea + area; }
        else if (nIn >= 2) { cls = "small cluster"; colour = "orange"; nClu = nClu + 1; }
        else if (nIn == 0) { cls = "no cell"; colour = "gray"; }
        else nIso = nIso + 1;
        eqD = 2 * sqrt(area / PI);
        if (cls != "no cell") {
            File.append(prefix + "," + channel + "," + fname + "," + (g + 1) + "," + cls + "," + nIn + "," + area + "," + eqD + "," +
                        feretV + "," + minF + "," + mean + "," + intden + "," + d2s(satFrac, 3), grpF);
        }
        File.append((g + 1) + "," + cls + "," + nIn + "," + area + "," + mean + "," + intden, infoF);
        if (cls == "spheroid") { selectImage(sphID); roiManager("select", g); setColor(255); fill(); run("Select None"); }
        if (cls == "single cell") { selectImage(isoID); roiManager("select", g); setColor(255); fill(); run("Select None"); }
        selectImage(gchkID);
        roiManager("select", g);
        Overlay.addSelection(colour, 2);
        run("Select None");
    }

    if (nG > 0) roiManager("Save", sDir + base + "_groups_rois.zip");
    selectImage(sphID);
    saveAs("Tiff", sDir + base + "_spheroid_mask.tif");
    close();
    selectImage(isoID);
    saveAs("Tiff", sDir + base + "_singlecells_mask.tif");
    close();
    selectImage(gchkID);
    roiManager("Deselect");
    run("Select None");
    run("Flatten");
    saveAs("Jpeg", sDir + base + "_groups_check.jpg");
    close();
    selectImage(gchkID); close();
    selectImage(cmID); close();
    roiManager("reset");
    return "" + nSph + "|" + nClu + "|" + nIso + "|" + d2s(sphArea, 0);
}

// ================================ object-aware background subtraction ================================
// Estimates the background ONLY from pixels that are not part of an object (objects = brighter than the image median by
// bgMargin, grown by bgMaskGrow px), using a normalised Gaussian blur, so bright spheroids with a dimmer core keep that core.
// (A rolling ball treats a dim core surrounded by a bright rim as background and subtracts it to black.)
// Returns the ID of the new background-subtracted 32-bit image titled "sub".
function subtractMaskedBackground(srcID) {
    selectImage(srcID);
    run("Select None");
    run("Duplicate...", "title=mk");
    run("Gaussian Blur...", "sigma=3");
    med = getValue("Median");
    setThreshold(med + bgMargin, 1000000000);
    run("Convert to Mask");
    if (bgMaskGrow > 0) run("Maximum...", "radius=" + bgMaskGrow);
    run("32-bit");
    run("Divide...", "value=255");
    run("Multiply...", "value=-1");
    run("Add...", "value=1");
    wID = getImageID();
    imageCalculator("Multiply create 32-bit", srcID, wID);
    numID = getImageID();
    run("Gaussian Blur...", "sigma=" + bgSigma);
    run("Add...", "value=" + (0.001 * med));
    selectImage(wID);
    run("Gaussian Blur...", "sigma=" + bgSigma);
    run("Add...", "value=0.001");
    imageCalculator("Divide create 32-bit", numID, wID);
    bgID = getImageID();
    imageCalculator("Subtract create 32-bit", srcID, bgID);
    subOutID = getImageID();
    run("Min...", "value=0");
    rename("sub");
    selectImage(numID); close();
    selectImage(wID); close();
    selectImage(bgID); close();
    return subOutID;
}

// ================================ EthD spot detector ================================
// Finds small bright spots in the background-subtracted red image and leaves them as ROIs in the ROI Manager.
// Returns "nSpots|threshold|noiseSD".
function detectSpots(subID) {
    selectImage(subID);
    run("Select None");
    run("Duplicate...", "title=sa");
    run("Gaussian Blur...", "sigma=" + spotSigma);
    saID = getImageID();
    selectImage(subID);
    run("Duplicate...", "title=sb");
    run("Gaussian Blur...", "sigma=" + (spotSigma * 2.5));
    sbID = getImageID();
    imageCalculator("Subtract create 32-bit", saID, sbID);
    dID = getImageID();
    selectImage(saID); close();
    selectImage(sbID); close();
    // robust noise level of the band-passed image: median and median absolute deviation
    selectImage(dID);
    med = getValue("Median");
    run("Duplicate...", "title=dabs");
    run("Subtract...", "value=" + med);
    run("Abs");
    mad = getValue("Median");
    close();
    sigR = 1.4826 * mad;
    if (sigR < 0.05) sigR = 0.05;
    thr = maxOf(med + spotK * sigR, spotMinContrast);
    selectImage(dID);
    setThreshold(thr, 1000000000);
    roiManager("reset");
    run("Analyze Particles...", "size=" + spotMinArea + "-" + spotMaxArea + " pixel add");
    resetThreshold();
    n = roiManager("count");
    selectImage(dID); close();
    return "" + n + "|" + d2s(thr, 2) + "|" + d2s(sigR, 3);
}

// ================================ EthD inside the calcein groups ================================
// Opens the calcein group outlines, measures the (leak-corrected, background-subtracted) red image inside each one and writes
// groups_ethd.csv. Returns a short note.
function measureGroupsEthD(redID, zipF, infoF, prefix, ethdFile, calFile, sDir) {
    roiManager("reset");
    roiManager("Open", zipF);
    nG = roiManager("count");
    info = split(File.openAsString(infoF), "\n");     // line 0 = header, line g+1 = group g+1
    outF = sDir + "groups_ethd.csv";
    for (g = 0; g < nG; g++) {
        f = split(replace(info[g + 1], "\r", ""), ",");   // GroupID,Class,NumCells,Area,Mean,IntDen
        if (f[1] == "no cell") continue;
        selectImage(redID);
        roiManager("select", g);
        List.setMeasurements;
        inArea = List.getValue("Area");
        inMean = List.getValue("Mean");
        run("Enlarge...", "enlarge=" + ringPx + " pixel");
        List.setMeasurements;
        bigArea = List.getValue("Area");
        bigMean = List.getValue("Mean");
        run("Select None");
        ringMean = inMean;
        if (bigArea > inArea) ringMean = (bigMean * bigArea - inMean * inArea) / (bigArea - inArea);
        netMean = inMean - ringMean;
        cMean = parseFloat(f[4]);
        ratio = "NA";
        if (cMean > 0) ratio = netMean / cMean;
        File.append(prefix + "," + calFile + "," + ethdFile + "," + (g + 1) + "," + f[1] + "," + f[2] + "," + inArea + "," +
                    f[4] + "," + f[5] + "," + inMean + "," + ringMean + "," + netMean + "," + (netMean * inArea) + "," + ratio, outF);
    }
    roiManager("reset");
    return "EthD measured in " + nG + " calcein groups";
}

// ================================ shrink the cell masks ================================
// Shrinks every ROI in the ROI Manager by px pixels (removes the blurry rim Cellpose tends to include). ROIs that vanish are deleted.
function shrinkRois(px) {
    n = roiManager("count");
    for (i = n - 1; i >= 0; i--) {
        roiManager("select", i);
        run("Enlarge...", "enlarge=-" + px + " pixel");
        if (selectionType() == -1) roiManager("delete");
        else roiManager("update");
    }
    run("Select None");
}

// ================================ label image -> ROI Manager ================================
function labelsToRois(labID) {
    selectImage(labID);
    getStatistics(la, lmean, lmin, lmax);
    print("  Label image: max label = " + lmax);
    roiManager("reset");
    for (v = 1; v <= lmax; v++) {
        selectImage(labID);
        if (bitDepth() == 32) setThreshold(v - 0.5, v + 0.5); else setThreshold(v, v);
        run("Create Selection");
        if (selectionType() != -1) roiManager("Add");
        run("Select None");
    }
    selectImage(labID);
    resetThreshold();
}

function writeImageRow(imgF, prefix, channel, fname, nCells, nKept, nSingle, nCluster, keptArea, keptInt, imgArea, bK, bA, bR, note) {
    note = replace(note, ",", ";");      // a comma in the note would shift the CSV columns
    line = prefix + "," + channel + "," + fname + "," +
           nCells + "," + nKept + "," + nSingle + "," + nCluster + "," + keptArea + "," + keptInt + "," + imgArea + "," +
           bK + "," + bA + "," + bR + "," + note;
    File.append(line, imgF);
}

// ================================ analyse one image ================================
function analyseImage(path, channel, greenPath, note0, prefix, sDir, sDay) {
    fname = File.getName(path);
    snap = replace(replace(fname, ".tiff", ""), ".tif", "");
    cellsF = sDir + "cells_measurements.csv";
    imgF   = sDir + "per_image_summary.csv";
    note = note0;
    bK = ""; bA = ""; bR = "";
    print("Processing: " + prefix + " | " + fname + " [" + channel + "]");

    // a completely black red image carries no information (not even background): record NA, not "0 dead cells"
    if (channel == "EthD" && indexOf(note0, "completely dark") >= 0) {
        print("  WARNING: " + fname + " is completely dark - EthD not analysed (recorded as NA)");
        if (!(cellposeMode == "external" && stage == 1)) writeImageRow(imgF, prefix, channel, fname, "NA", "NA", "NA", "NA", "NA", "NA", "NA", "", "", "", note0 + " | EthD image is completely dark: NOT analysed");
        return;
    }

    // per-day settings
    bgUse = bgSize; diamUse = cellDiameter; grpUse = 0; bgDay = bgSize;
    found = false;
    for (d = 0; d < dayNames.length; d++) {
        if (toUpperCase(sDay) == toUpperCase(dayNames[d])) { bgUse = dayBgSize[d]; bgDay = dayBgSize[d]; diamUse = dayDiameter[d]; grpUse = dayGroup[d]; found = true; }
    }
    if (!found) print("  WARNING: no per-day settings for '" + sDay + "', using defaults");
    if (channel == "EthD" && ethdUseSingleCellSettings) { bgUse = bgSizeEthD; diamUse = cellDiameter; }

    lc = split(loadChannel(path, channel, "src"), "|");
    srcID = parseFloat(lc[0]);
    note = note + " | " + lc[1];
    selectImage(srcID);
    getPixelSize(u0, pw0, ph0);

    // ---- EthD: subtract calcein bleed-through using this set's own green image ----
    if (channel == "EthD" && doBleedCorrection) {
        if (greenPath == "") {
            note = note + " | no unique green image in set: NO bleed-through correction";
        } else {
            lg = split(loadChannel(greenPath, "calcein", "grn"), "|");
            gID = parseFloat(lg[0]);
            selectImage(srcID); w1 = getWidth(); h1 = getHeight();
            selectImage(gID);   w2 = getWidth(); h2 = getHeight();
            if (w1 != w2 || h1 != h2) {
                note = note + " | red and green differ in size: NO bleed-through correction";
                selectImage(gID); close();
            } else {
                fb = split(fitBleed(gID, srcID), "|");
                bleedK = parseFloat(fb[0]); bleedA = parseFloat(fb[1]); bleedR = parseFloat(fb[2]);
                bK = bleedK; bA = bleedA; bR = bleedR;
                if (bleedK > maxBleedK) {
                    note = note + " | LEAKAGE TOO STRONG (k=" + d2s(bleedK, 3) + " r=" + d2s(bleedR, 2) + "): EthD NOT analysed";
                    print("  WARNING: red image is mostly calcein leaking through (k=" + d2s(bleedK, 3) + "), EthD not analysed");
                    selectImage(gID); close();
                    writeImageRow(imgF, prefix, channel, fname, "NA", "NA", "NA", "NA", "NA", "NA", w1 * pw0 * h1 * ph0, bK, bA, bR, note);
                    run("Close All");
                    return;
                }
                if (bleedR >= minBleedCorr && bleedK > 0) {
                    selectImage(gID);
                    run("Multiply...", "value=" + bleedK);
                    run("Add...", "value=" + bleedA);
                    imageCalculator("Subtract create 32-bit", srcID, gID);
                    resID = getImageID();
                    run("Min...", "value=0");
                    selectImage(srcID); close();
                    selectImage(gID); close();
                    selectImage(resID);
                    rename("src");
                    srcID = resID;
                    setPixelSize(pw0);
                    note = note + " | bleed-through corrected k=" + d2s(bleedK, 4) + " r=" + d2s(bleedR, 2);
                } else {
                    selectImage(gID); close();
                    note = note + " | bleed fit poor (r=" + d2s(bleedR, 2) + " k=" + d2s(bleedK, 4) + "): NO correction";
                }
            }
        }
    }

    // ---- spheroid-level EthD: measure the red signal inside each calcein spheroid / cluster / single-cell mask ----
    // (measured on the leak-corrected red image with a local background taken from a ring around each group)
    if (channel == "EthD" && greenPath != "") {
        calName = File.getName(greenPath);
        calBase = "calcein_" + replace(replace(calName, ".tiff", ""), ".tif", "");
        zipF = sDir + calBase + "_groups_rois.zip";
        infoF = sDir + calBase + "_groups_info.csv";
        if (File.exists(zipF) && File.exists(infoF)) {
            selectImage(srcID);
            run("Select None");
            note = note + " | " + measureGroupsEthD(srcID, zipF, infoF, prefix, fname, calName, sDir);
        } else {
            note = note + " | no calcein groups for this set (nothing to measure EthD in)";
        }
    }

    // ---- background subtraction (rolling ball / sliding paraboloid, smoothing disabled) ----
    selectImage(srcID);
    if (channel == "calcein" && bgMethodCalcein == "masked") {
        subID = subtractMaskedBackground(srcID);
    } else {
        run("Duplicate...", "title=sub");
        run("Subtract Background...", "rolling=" + bgUse + " sliding disable");
        subID = getImageID();
    }
    setPixelSize(pw0);
    getPixelSize(pu, pw, ph);
    imgArea = getWidth() * pw * getHeight() * ph;

    useSpots = (channel == "EthD" && ethdMethod == "spots");

    // ---- "no cells" check: enough pixels above the peak threshold? ----
    peakThr = minPeakCalcein;
    if (channel == "EthD") peakThr = minPeakEthD;
    selectImage(subID);
    getHistogram(hv, hc, 256, 0, 256);
    nAbove = 0;
    for (b = floor(peakThr); b < 256; b++) nAbove = nAbove + hc[b];
    note = note + " | px above " + peakThr + ": " + nAbove;
    if (nAbove < minPixAbove && !useSpots) {
        if (!(cellposeMode == "external" && stage == 1)) writeImageRow(imgF, prefix, channel, fname, 0, 0, 0, 0, 0, 0, imgArea, bK, bA, bR, note + " | NO SIGNAL");
        print("  No signal -> recorded 0 cells");
        run("Close All");
        return;
    }

    labID = 0; labTitle = "";
    if (useSpots) {
        // ---- EthD: small bright spots (stained nuclei) ----
        roiManager("reset");
        sp = split(detectSpots(subID), "|");
        note = note + " | spots: " + sp[0] + " found (threshold " + sp[1] + ", noise SD " + sp[2] + ")";
        print("  EthD spots: " + sp[0] + " (threshold " + sp[1] + ", noise SD " + sp[2] + ")");
    } else {
        // ---- segmentation copy (noise floor removed, x10 so small values survive the 16-bit conversion) ----
        selectImage(subID);
        run("Duplicate...", "title=seg");
        floorV = segFloorCalcein;
        if (channel == "EthD") floorV = segFloorEthD;
        if (floorV > 0) run("Subtract...", "value=" + floorV);
        run("Min...", "value=0");
        run("Multiply...", "value=10");
        run("16-bit");
        segID = getImageID();

        // ---- optional focus filter: black out out-of-focus objects in the copy Cellpose sees ----
        if (focusFilter) {
            ff = split(removeOutOfFocus(subID, segID, sDir + channel + "_" + snap + "_focus_check.jpg"), "|");
            note = note + " | focus filter: " + ff[0] + " objects, " + ff[1] + " removed as out of focus";
            print("  focus filter: " + ff[0] + " objects, " + ff[1] + " removed");
        }

        if (cellposeMode == "external") {
            // ---- Cellpose runs outside Fiji: save the image (step 1) or read the masks it made (step 2) ----
            segRoot = rootDir + "ToSegment" + File.separator;
            segFolder = segRoot + "d" + diamUse + File.separator;
            pfx = substring(prefix, 0, lastIndexOf(prefix, ","));          // Group,Day,Scaffold,Set (without the complete flag)
            segBase = replace(replace(replace(pfx, ",", "_"), " ", ""), "__", "_") + "_" + snap;
            segPath = segFolder + segBase + ".tif";
            if (stage == 1) {
                File.makeDirectory(segRoot);
                File.makeDirectory(segFolder);
                selectImage(segID);
                saveAs("Tiff", segPath);
                print("  saved for Cellpose: " + segPath);
                run("Close All");
                return;
            }
            maskPath = segFolder + segBase + "_cp_masks.tif";
            if (!File.exists(maskPath) && File.exists(segFolder + segBase + "_masks.tif")) maskPath = segFolder + segBase + "_masks.tif";
            if (!File.exists(maskPath)) {
                print("  ERROR: no Cellpose mask found for " + fname + " (looked for " + maskPath + ")");
                writeImageRow(imgF, prefix, channel, fname, "NA", "NA", "NA", "NA", "NA", "NA", imgArea, bK, bA, bR, note + " | ERROR no Cellpose mask file");
                run("Close All");
                return;
            }
            open(maskPath);
            labID = getImageID();
            labTitle = getTitle();
            wl = getWidth(); hl = getHeight();
            selectImage(subID);
            if (wl != getWidth() || hl != getHeight()) {
                print("  ERROR: mask size does not match the image for " + fname);
                writeImageRow(imgF, prefix, channel, fname, "NA", "NA", "NA", "NA", "NA", "NA", imgArea, bK, bA, bR, note + " | ERROR mask size differs from image");
                run("Close All");
                return;
            }
            print("  read Cellpose masks: " + segBase + "_cp_masks.tif");
            labelsToRois(labID);
        } else {
            // ---- Cellpose ----
            roiManager("reset");
            before = newArray(nImages);
            for (q = 1; q <= nImages; q++) { selectImage(q); before[q - 1] = getImageID(); }
            selectImage(segID);
            cpArgs = "cp_model=cyto3 custom_model= cell_diameter=" + diamUse + " cyto_channel=1 nuclei_channel=None min_size=" + minSize +
                     " normalize=true resample=true return_rois=" + cellposeReturnRois + " cellprob_threshold=" + cellprobThr + " flow_threshold=" + flowThr +
                     " tile_overlap=0.1 niter=0 compute_flows=false shuffle=true labeltype=16-bit mode_3d=None stitch_threshold=0.0" +
                     " flow3d_smooth=0 torchversion=cpu usegpu=false";
            run(cellposeCmd, cpArgs);

            labID = 0; labTitle = "";
            for (q = 1; q <= nImages; q++) {
                selectImage(q);
                isNew = true;
                for (w = 0; w < before.length; w++) if (before[w] == getImageID()) isNew = false;
                if (isNew && bitDepth() != 24) { labID = getImageID(); labTitle = getTitle(); }
            }
            if (labID == 0 && roiManager("count") == 0) {
                print("  ERROR: Cellpose returned no new image and no ROIs for " + fname + " - skipped");
                writeImageRow(imgF, prefix, channel, fname, 0, 0, 0, 0, 0, 0, imgArea, bK, bA, bR, note + " | ERROR no Cellpose output");
                run("Close All");
                return;
            }
            if (labID != 0) print("  Cellpose returned image: " + labTitle);
            if (cellposeReturnRois != "true") labelsToRois(labID);
        }
    }
    if (!useSpots && maskShrinkPx > 0 && roiManager("count") > 0) shrinkRois(maskShrinkPx);
    nCells = roiManager("count");

    base = channel + "_" + snap;
    nKept = 0; nSingle = 0; nCluster = 0; keptArea = 0; keptInt = 0;
    cutMean = minCellMeanCalcein;
    if (channel == "EthD") cutMean = minCellMeanEthD;
    if (nCells > 0) {
        // measure each ROI on its own (does not depend on the Results table keeping all rows)
        for (r = 0; r < nCells; r++) {
            selectImage(subID);
            roiManager("select", r);
            List.setMeasurements;
            area = List.getValue("Area");
            mean = List.getValue("Mean");
            intden = List.getValue("IntDen");
            feretV = List.getValue("Feret");
            minF = List.getValue("MinFeret");
            run("Select None");
            eqD = 2 * sqrt(area / PI);
            keep = 1;
            if (mean < cutMean || area < minArea) keep = 0;
            if (channel == "EthD" && maxAspectEthD > 0) {
                if (minF > 0 && feretV / minF > maxAspectEthD) keep = 0;
            }
            nKept = nKept + keep;
            if (keep == 1) {
                keptArea = keptArea + area;
                keptInt = keptInt + intden;
                if (eqD < clusterEqDiam) nSingle = nSingle + 1; else nCluster = nCluster + 1;
            }
            line = prefix + "," + channel + "," + fname + "," + (r + 1) + "," +
                   area + "," + mean + "," + intden + "," +
                   feretV + "," + minF + "," + eqD + "," + keep;
            File.append(line, cellsF);
        }
        roiManager("Save", sDir + base + "_rois.zip");
        // QC picture: the blank-subtracted image with the cell outlines drawn on it
        selectImage(subID);
        run("Duplicate...", "title=chk");
        run("Enhance Contrast", "saturated=0.35");
        run("From ROI Manager");
        run("Flatten");
        saveAs("Jpeg", sDir + base + "_check.jpg");
        // group touching cells into spheroids / clusters (calcein only; cells are still in the ROI Manager)
        if (channel == "calcein" && grpUse == 1) {
            gs = split(groupSpheroids(subID, nCells, prefix, channel, fname, sDir, base, srcID), "|");
            note = note + " | grouping: " + gs[0] + " spheroids, " + gs[1] + " small clusters, " + gs[2] + " isolated cells";
            print("  grouping: " + gs[0] + " spheroids (total area " + gs[3] + "), " + gs[1] + " small clusters, " + gs[2] + " isolated cells");
        }
    } else {
        print("  No cells found.");
    }
    writeImageRow(imgF, prefix, channel, fname, nCells, nKept, nSingle, nCluster, keptArea, keptInt, imgArea, bK, bA, bR, note);

    if (labID != 0) {
        selectImage(labID);
        saveAs("Tiff", sDir + base + "_labels.tif");
    }
    run("Close All");
    roiManager("reset");
}

// ================================ combined files in ROOT ================================
function buildMasters() {
    if (!File.exists(registry)) { print("No sets were processed, no combined files written."); return; }
    dirArr = split(File.openAsString(registry), "\n");
    kinds = newArray("cells_measurements.csv", "per_image_summary.csv", "set_report.csv", "groups_measurements.csv", "groups_ethd.csv");
    for (kk = 0; kk < kinds.length; kk++) {
        acc = "";
        first = true;
        for (s = 0; s < dirArr.length; s++) {
            f = replace(dirArr[s], "\r", "") + kinds[kk];
            if (!File.exists(f)) continue;
            lines = split(File.openAsString(f), "\n");
            for (j = 0; j < lines.length; j++) {
                ln = replace(lines[j], "\r", "");
                if (lengthOf(ln) == 0) continue;
                if (j == 0 && !first) continue;      // header only once
                acc = acc + ln + "\n";
            }
            first = false;
        }
        outP = rootDir + "ALL_" + kinds[kk];
        if (File.exists(outP)) ok = File.delete(outP);
        if (File.exists(outP)) {
            outP = rootDir + "ALL_NEW_" + kinds[kk];
            print("WARNING: could not replace ALL_" + kinds[kk] + " (open in Excel?). Saved as " + outP + " instead.");
        }
        File.saveString(acc, outP);
    }
}
