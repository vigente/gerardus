# Table of Contents

   * [Overview](#overview)
      * [Source code](#source-code)
      * [Data](#data)
      * [Manual landmarks and masks](#manual-landmarks-and-masks)
   * [Summary of scripts](#summary-of-scripts)
      * [paper_figures.m](#paper_figuresm)
      * [blockface_mouse.m](#blockface_mousem)
      * [manual_landmarks.m](#manual_landmarksm)
      * [histology_blockface_mouse.m](#histology_blockface_mousem)
      * [study_stopping_criterion_rigid_refinement.m](#study_stopping_criterion_rigid_refinementm)
      * [study_stopping_criterion_bsp_refinement.m](#study_stopping_criterion_bsp_refinementm)
      * [histology_noblockface_mouse.m](#histology_noblockface_mousem)
      * [histology_refinement_without_diffusion.m](#histology_refinement_without_diffusionm)

Created by [gh-md-toc](https://github.com/ekalinin/github-markdown-toc)

# Overview

## Source code

In this directory we provide the Matlab scripts that we used for the experiments, figures and results in
the paper

Casero et al. "Transformation Diffusion Reconstruction of
Three-Dimensional Histology Volumes from Two-Dimensional Image
Stacks", Medical Image Analysis, 2017, DOI: [10.1016/j.media.2017.03.004](https://doi.org/10.1016/j.media.2017.03.004).

http://www.sciencedirect.com/science/article/pii/S1361841517300397

These scripts depend on [Matlab code from the gerardus project](https://github.com/vigente/gerardus).
You can clone the code and add the toolboxes to your Matlab path following the "Quickstart" section in the [instructions to build and install Gerardus](https://github.com/rcasero/gerardus/wiki/Build-instructions). (Basically, you'll need to clone the gerardus project, run `gerardus/matlab/add_gerardus_paths.m`, and install the [elastix binaries](http://elastix.isi.uu.nl/).

Note: Strictly speaking, the paper above was run with [Gerardus tag "casero2017_3d_histology_transformation_diffusion_reconstruction"](https://github.com/vigente/gerardus/releases/tag/casero2017_3d_histology_transformation_diffusion_reconstruction). 
However, it's quite likely that later versions of Gerardus will work with these scripts too. It's worth giving a try to the current master version.

## Data

The original histology stack can be downloaded from ["Serial histology and blockface images. Full mouse heart, 600 sections"](https://ora.ox.ac.uk/objects/uuid:75f09b9e-e7d5-4b48-9519-619177cea1ef), DOI: [10.5287/bodleian:bpM4PmPvo](https://doi.org/10.5287/bodleian:bpM4PmPvo).

The reconstructed histology with our method can be downloaded from ["3D histology reconstruction with external blockface reference. Full mouse heart, 600 sections"](https://ora.ox.ac.uk/objects/uuid:716fe2ef-f965-40ff-8da4-8a5ad48f4aea), DOI: [10.5287/bodleian:o8eNyrzbX](https://doi.org/10.5287/bodleian:o8eNyrzbX).

## Manual landmarks and masks

In the `hand_tracing/` directory we provide landmarks and masks
that we traced by hand, so that readers can reproduce our results. The
data (blockface and histology images) are not provided in this
directory, but in the location stated in the paper.

# Summary of scripts

A list of all main scripts and their purpose can be found below. We
ran the scripts on Matlab 8.3.0.532 (R2014a).

The scripts are not meant to be invoked by name and run fully, as they
contain debug code, code to hand trace landmarks and masks,
etc. Instead, the user should run individual sections understanding
what each part is doing.

Note also that variables that point to the location of the data will
need to be edited to suit the local configuration of each user.

## paper_figures.m

Script to generate misc figures:

* **Figure 10.** TDR/ATDR smoothing kernel Ψ for several number of repeated 
convolutions or stack diffusion sweeps m.
 * `%% repeat convolution of solution kernels`

* **Figure 4.** Transformation Diffusion Reconstruction (TDR) applied to a 1D 
translation synthetic example. 
 * `%% visualization of transformation diffusion with vertical translation (solution is sine wave)`
 * `%% synthetic problem, comparison with Gaffling 2015`
 * `%% synthetic problem adding drift, comparison with Gaffling 2015`

## blockface_mouse.m

Script to align, preprocess, etc the two blockface volumes (0 deg and 55 deg) of
the whole mouse heart.

## manual_landmarks.m

Script to hand trace landmarks in pairs of

* blockface / histology slices
* histology / histology slices

to validate the reconstruction algorithms.

## histology_blockface_mouse.m

Script to reconstruct Mouse Q53 heart's low-res and hi-res histology,
using blockface as external reference.

It uses the auxiliary `aux_histology2blockface.m` script to perform the 
histology-blocface registration.

Note that when using Transformation Diffusion, this script only computes the 
first stack sweep (M = 1, the one that involves registrations). Further 
refinements with more stack sweeps (M>1) are done in the scripts below. 

The reason to do it this way is for clarity of the code. In a real
reconstruction, you would just choose one value of M. But for the paper, we had
to test and generate figure for the results with different values of M.

## study_stopping_criterion_rigid_refinement.m

Script to check how landmark error changes with the number of stack
sweeps for low-res rigid refinement.

This script works with the output of histology_blockface_mouse.m, section
"Intra-histology refinement: rigid diffusion".

## study_stopping_criterion_bsp_refinement.m

Script to check how landmark error changes with the number of stack
sweeps for low-res and hi-res B-spline refinement.

This script works with the output of histology_blockface_mouse.m, section
"Intra-histology registration refinement: B-spline diffusion" (low-res)
and section "Refine reconstruction in high resolution" (hi-res).

## histology_noblockface_mouse.m

Script to reconstruct Mouse Q53 heart's low-res histology, without an
external reference, to compare to "histology_blockface_mouse.m". This
corresponds to experiment "2.2.1.	Low resolution reconstruction without
external blockface reference" in the paper.

## histology_refinement_without_diffusion.m

Script to replace the transformation diffusion low-res B-spline
refinement in "histology_blockface_mouse.m" by refinement using
registration sweeps of the stack, without using diffusion. This
result corresponds to the "Reg. only" curve in Fig. 8 in the paper.
