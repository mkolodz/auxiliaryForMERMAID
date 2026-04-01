Contents:

Mapping:
mappingSimVsExpAdjusted_correct.m - reads in a set of experimental and simulated data (cylinder phantom? but not confirmed), remaps the experimental data to match the simulation mapping, plots 2D histograms of counts per LOR and counts per pixel for both simulation and experiment. ALSO, at the end uses the remapped data as input to normalization components estimation (still work in progress). Needs the two functions: createSetOfEquivalentLORs_v3() and estimate_normalization_components_v2().
mappingSimToExpOnly.m - producing a mapping table between simulation and experiment.

Normalization:
createSetOfEquivalentLORs_v3.m - more LORs taken into account than in  createSetOfEquivalentLORs_v2.m, v3 is probablu correct according to Kinouchi et al.

Time alignment:
timeAlignmentAnalysis.m - ongoing work, so far only read in preprocessed data from time alignment measurement and plot 2D histogram of counts per LOR with raw indices.

PETsys parameter scan for new MERMAID modules:
newPETModules_parameterScan.m - plot comparison spectra for various parameters of PETsys. requires mermaid code to work.

