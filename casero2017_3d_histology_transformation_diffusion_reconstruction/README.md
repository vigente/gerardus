# Overview

In this "src/" directory we provide the Matlab scripts and
dependencies that we used for the experiments, figures and results in
the paper

Casero et al. "Transformation Diffusion Reconstruction of
Three-Dimensional Histology Volumes from Two-Dimensional Image
Stacks", Medical Image Analysis, 2017.

http://www.sciencedirect.com/science/article/pii/S1361841517300397

In the "src/hand_tracing/" directory we provide landmarks and masks
that we traced by hand, so that readers can reproduce our results. The
data (blockface and histology images) are not provided in this
directory, but in the location stated in the paper.

A list of all main scripts and their purpose can be found below. We
ran the scripts on Matlab 8.3.0.532 (R2014a).

The scripts are not meant to be invoked by name and run fully, as they
contain debug code, code to hand trace landmarks and masks,
etc. Instead, the user should run individual sections understanding
what each part is doing.

Note also that variables that point to the location of the data will
need to be edited to suit the local configuration of each user.

# Summary of scripts

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
% corresponds to the "Reg. only" curve in Fig. 8 in the paper.



grab_src_code_dependencies.m

histology_refinement_only_diffusion.m

validate_reconstruction.m
validate_reconstruction_noblockface.m
validate_reconstruction_onlydiff.m




* 'histology_blockface_mouse.m'

% Script to reconstruct one half of Mouse Q53 heart's low-res histology,
% using blockface as external reference, and detail of Left Ventricle (LV)
% free wall.
%
%   - Correct alignment of blockface stack.
%   - Pre-process blockface images.
%   - Rigid alignment of low-res histology to blockface.
%   - Rigid diffusion refinement of low-res intra-histology.
%   - B-spline diffusion refinement of low-res intra-histology.
%     (4 registration sweeps, 5 neighbor updates each).
%   - B-spline diffusion refinement of hi-res intra-histology, 
%     detail of LV free wall (4 registration sweeps, 
%     100 neighbor updates each).

* 'histology_noblockface_mouse.m'

% Script to reconstruct one half of Mouse Q53 heart's low-res histology,
% without an external reference, to compare to
% "histology_blockface_mouse.m".
%
%   - Sequential pre-alignment of low-res histology.
%   - Rigid diffusion refinement of low-res intra-histology.
%   - B-spline diffusion refinement of low-res intra-histology 
%     (4 registration sweeps, 5 neighbor updates each).

* 'histology_refinement_only_diffusion.m'

% Script to replace the low-res B-spline refinement in
% "histology_blockface_mouse.m":
%
%   - B-spline diffusion refinement of low-res intra-histology.
%     (1 registration sweeps, 1, 20, 40, 60, 80 and 100 neighbor updates).

* 'histology_refinement_without_diffusion.m'

% Script to replace the low-res B-spline refinement in
% "histology_blockface_mouse.m":
%
%   - B-spline diffusion refinement of low-res intra-histology.
%     (20 registration sweeps, no neighbor updates).

* 'manual_landmarks.m'

% Script to hand trace landmarks in pairs of
%   - blockface / histology slices
%   - histology / histology slices
% to validate the reconstruction algorithms.

* 'toy_examples.m'

% Obsolete, not used in final paper.

* 'validate_reconstruction.m'

% This is a validation script. We assume:
%
%   * a low-res histology reconstruction has been computed with
%     histology_blockface_mouse.m, histology_refinement_without_diffusion.m
%     or a similar script.
%
%   * manual landmarks between blockface/histology and histology/histology
%     slices have been traced using script manual_landmarks.m.
%
% Then, this script will apply the reconstruction transforms to the
% landmarks, and compute how far they are from each other at each
% refinement step.
%
% All the user has to do is choose below which transforms are going to be
% validated.
%
% This script can be reused to validate, e.g.
%   * a naive algorithm without diffusion
%   * our diffusion reconstruction
%   * an algorithm with only 1 registration sweep and a lot of diffusion
%   * ...
%
% Validation of the "diffusion only" experiment is in another script
% (validate_reconstruction_onlydiff.m).

* 'validate_reconstruction_noblockface.m'

% This is a validation script. We assume that
%
%   * a histology reconstruction has been computed with
%     histology_noblockface_mouse.m
%
%   * manual landmarks between histology/histology
%     slices have been traced using script manual_landmarks.m.
%
% Then, this script will apply the reconstruction transforms to the
% landmarks, and compute how far they are from each other at each
% refinement step.
%
% Note that we are validating the reconstruction without blockface, but we
% still need the blockface landmarks so that we can measure how much the
% reconstruction differs from the true anatomy.

* 'validate_reconstruction_onlydiff.m'

% Similar script to validate_reconstruction.m, but this one is for the
% "diffusion only" experiment, because instead of having a single vector of
% transforms, we have a cell vector with the transforms for each level of
% diffusion.

* 'paper_figures.m'

% Script to generate the plots/figures used in the paper using the
% images and validation data generated with the scripts above.
