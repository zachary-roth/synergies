# synergies

## Project abstract

De Groote et al. (2014) showed that solving the muscle redundancy problem by optimizing performance assuming independent recruitment of muscles in healthy adults resulted in muscle activations that on average could be explained by four muscle synergies1 . Those synergies are similar to the ones derived from measured experimentally using EMG.

Children with Cerebral Palsy (CP) however are known to have reduced motor control and use less synergies during walking2,3 . It is unclear if optimizing performance assuming independent recruitment of muscles to solve the redundancy problem in these children will result in similar synergies measured experimentally.

## The goals of this project are:

1. Determine the differences between muscle synergies based on estimated and measured muscle activity during walking in a child with CP.

2. Determine muscle synergies from measured muscle activity for a group of children with CP (8 subjects). Compare them and select 3 relevant subjects for further study (repeat analysis in 1).

3. Determine model sensitivity (marker-based model vs MRI model) towards estimated muscle synergies (subjects depending on available MRI models).

4. If time allows it, determine the differences between estimated and measured muscle synergies during functional movements (squat, sit-to-stand-to-sit and counter movement jump) in a child with CP.

## synergies1:  estimated vs measured muscle synergies in 1 participant, **Generic** Musculoskeletal Model
___
### Workflow

1. run `MoveData1.m`
  - process EMG data
      - trim trials based on ICs
      - demean
      - notch Filter (50Hz)
      - bandpass Filter (20-400 Hz)
      - rectify
      - lowpass Filter (15Hz)
  - find the max EMG values across all movements and trials within one participant
      - normalize based on max values
      - resample to 0-100% of movement trial
  - find the calculated muscle activations using the [Muscle Redundancy Solver](https://github.com/KULeuvenNeuromechanics/MuscleRedundancySolver)
  - resample to 0-100% of movement trial
  - write the resampled data for all subjects to the `MoveData.m` structure, which contains the meta data and EMG data for all subjects
2. run `nmf1.m`
  - concatenate the resampled movement data into one trial per movement
  - create a `calcReduced` dataset which contains only the calculated activations corresponding to the observed activations
  - use **non-negative matrix factorization** to calculate 1-6 synergies per movement
  - calculate the **variance-accounted-for** for each number of synergies
  - select the lowest number of synergies which account for >90% of the variance
  - write the synergies and VAF data to the `nmf.mat` structure
3. run `Visualizations1.m`
  - visualize the lowpass filtered EMG data (look for large spikes)
  - visualize the observed activations vs the reconstructed activations (look for correspondance)
  - visualize the synergy activations and muscle weightings for each number of synergies (1-6)
  - visualize the VAF overall and the VAF per muscle for each number of synergies (1-6)

`Figures1.m` is a hard-coded script to generate the plots used in presentations

## synergies2:  measured muscle synergies in 7 participants
___
### Workflow

1. extract trials meta data from protocol pdf *by hand* (CPx_trials.xlsx)
2. run `extract_EMG.m`: extracts the raw EMG signals from the .c3d files as .xlsx files
3. run `get_ICs.m`: 
  - automatically extract the start and stop times of non-gait movement trials
  - automatically extract the first intital contact (IC1) of gait trials
  - extract IC2 using a GUI based on heel marker data
  - extract marker trajectories (.trc file) to double check heel strike estimate
4. run `MoveData2.m`
  - identical to `MoveData1.m` but with MRS sections removed
5. run `nmf2.m`
 - identical to `nmf1.m` but with `calcReduced` sections removed
6. run `Visualizations1.m`

`Figures2.m` is a hard-coded script to generate the plots used in presentations

## synergies3:  estimated vs measured muscle synergies in 1 participant, **MRI-Based** Musculoskeletal Model
___
### Workflow

1. run `MoveData1.m`
  - Select the MRI-based model
2. run `nmf1.m`
3. run `Visualizations1.m`
