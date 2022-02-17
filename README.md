Sample code to reproduce results of the manuscript *"Exploring MEG brain fingerprints: Evaluation, pitfalls, and interpretations"* by Sareen et al., NeuroImage 2021  , https://doi.org/10.1016/j.neuroimage.2021.118331. 

# MEG Fingerprints

Run the main script: meg_id_icc.m. The code generates Fig.2 and Fig. 3  of the paper (only the Intra-class correlation matrices) for the alpha and beta bands, for 3 MEG Functional Connectivity measures (i.e., AECc, PLM, wPLI). More comment and details are provided within the main script. Also, some sample connectomes are included in the repo to make the code standalone. 
The main script "perm_main" runs the Permutation Testing (PT) framework to assess the statistical significance of the differential identifiability and success rate scores, as explained in (Sareen et. al 2021, NeuroImage).


# MEG Source reconstruction

The code to perform MEG source reconstruction to the 148 cortical centroids of the Destrieux atlas from preprocessed HCP data is availabe in the folder "MEG_source_recon". The scripts assume the HCP data structure and use the following softwares:
- connectome workbench: https://www.humanconnectome.org/software/get-connectome-workbench
- megconnectome, fieldtrip: https://www.humanconnectome.org/software/hcp-meg-pipelines
- cifti-matlab: https://github.com/Washington-University/cifti-matlab 

We used FieldTrip functions and the LCMV beamformer formulation. The main script is batch_ResampleParcellation_SourceReconstruction.m: you can use it by adapting paths in the first section.


# MEG Functional connectivity measures

The code to compute the functional connectomes is available in the folder "FCMethods". Specifically, the following FC metrics are included:

Amplitude coupling methods:

- AEC and AECc [amplitudeenvelopecorrelation.m]

Phase-coulping methods 
- wPLI [wPLI_adjmat.m]
- PLI [PLI.m]
- PLM [connectivity_plm.m] 
- PLV [PLV.m]

*Code Authors*: Ekansh Sareen, Alessandra Griffa, Enrico Amico.

Any comments/queries can be sent to: enrico.amico@epfl.ch

version 1.0 (15 July, 2021)

**Please cite us**! 
Ekansh Sareen, Sélima Zahar, Dimitri Van De Ville, Anubha Gupta, Alessandra Griffa, Enrico Amico, "Exploring MEG brain fingerprints: Evaluation, pitfalls, and interpretations", NeuroImage, Volume 240, 2021, 118331, ISSN 1053 8119, https://doi.org/10.1016/j.neuroimage.2021.118331.
