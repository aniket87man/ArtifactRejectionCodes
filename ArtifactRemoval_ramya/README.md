

This folder contains MATLAB code developed for artifact removal in ICMS-evoked local field potential (LFP) recordings.

The project focuses on removing stimulation-locked artifacts from LFP signals while preserving the underlying neural signal structure. Several artifact-removal approaches are implemented and compared, including ERP/template-based methods, PCA-based subtraction, SMARTALite variants, and a SMARTA-inspired adaptive artifact-removal method.

## Folder contents

```text
ArtifactRemoval_ramya/
├── methods/
│   └── MATLAB function files for artifact-removal methods
├── demos/
│   └── Demo, validation, and plotting scripts
├── Thesis.pdf
├── poster_final.pptx
└── README.md
```

## Method files

The `methods/` folder contains the main artifact-removal functions.

### ERP and template-based methods

```text
ERPSubtraction.m
ERPShifted.m
ERPScaled.m
ERPShiftedScaled.m
ERPAligned.m
ERPAlignedPulsewise.m
PulsewiseTemplate.m
```

These methods use stimulation-locked template subtraction. The simpler versions subtract an average ERP artifact template, while shifted, scaled, and aligned versions account for timing shifts and amplitude changes. Pulse-wise methods operate around individual stimulation pulses.

### PCA-based method

```text
PCATemplate.m
```

This method removes artifact-related low-dimensional structure using principal components estimated from the artifact window. It is useful as a stronger baseline, but it can be aggressive if stimulation-locked neural activity overlaps with artifact-related components.

### SMARTA-inspired methods

```text
SMARTALite.m
SMARTAShrinkLite.m
SMARTAFull.m
```

These methods use adaptive pulse-template estimation. Instead of using one fixed artifact template, they estimate local templates from similar pulse segments.

`SMARTAFull.m` is the main SMARTA-inspired method used in the final validation. It uses high-pass similarity features and optimal-shrinkage based denoising for improved artifact-template estimation.

## Demo and validation scripts

The `demos/` folder contains the numbered scripts used during method development, validation, and plotting.

The demo scripts document the development sequence from simple ERP subtraction to PCA-based methods, SMARTALite variants, over-cleaning diagnostics, semi-synthetic validation, and final across-electrode validation.

### Main final validation script

```text
demos/demo_29_FinalValidation.m
```

This is the main final validation script. It compares the final set of methods across selected V1 electrodes.

### Plotting and table-generation script

```text
demos/demo_23_FinalPlotsAndTables.m
```

This script generates plots and summary tables from saved validation metrics.

### Other important scripts

```text
demo_21_OverCleaningDiagnostics.m
demo_22_FinalValidationMetrics.m
demo_24_SemiSyntheticValidation.m
demo_25_TimeSeriesBeforeAfterCleaning.m
demo_28_OverCleaningAndPreservationValidation.m
demo_31_CleanFinalValidationPresentationPlots.m
```

These scripts contain additional diagnostics, metric calculations, semi-synthetic validation, before/after time-series visualization, and final summary plotting.

## Final comparison methods

The final validation compares the following methods:

```text
RawHighStim
ERPSubtraction
PCATemplate_K10
SMARTALite_Ensemble_K3K5
SMARTAFull_K3_hp100
```

The main final method is:

```text
SMARTAFull_K3_hp100
```

where `K = 3` is the number of nearest-neighbor pulse segments used for local template estimation, and `hp100` indicates a 100 Hz high-pass representation used for similarity matching.

## External SMARTA dependency

`SMARTAFull.m` uses optimal-shrinkage based denoising. For this functionality, the following external SMARTA-related files are required on the MATLAB path:

```text
optimal_shrinkage_color_fast.m
createPseudoNoise.m
```

These files are not included in this folder because they are from the original SMARTA-related codebase.

Please download the required external files from the original SMARTA repository:

```text
https://github.com/z123x698c547/Artifact-removal.git

```

After downloading the required external files, add their folder to the MATLAB path before running `SMARTAFull.m`.



Before running the demo scripts, update the local data paths in the scripts according to the machine and folder where the data is stored.

## Notes for users

- `demo_29_FinalValidation.m` should be used as the main final validation script.
- The numbered demo scripts are retained to show the development and validation sequence.
- Local paths in the scripts may need to be edited before running.
- `SMARTAFull.m` is a SMARTA-inspired implementation adapted for this ICMS-LFP artifact-removal project.
