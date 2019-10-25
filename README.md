# skullanalyzer
a concrete way of extracting cranial geometric features

The skullanalyzer loads a 3D triangular mesh cranial model along with a number of user defined landmarks and extracts a set of cranial geometric features, which can be used for further analyses. It requires two input files, one Alias Wavefront file (model_name.obj) containing a triangular mesh of the cranial model to be analyzed and a Meshlab PickedPoints (model_name.pp) sidecar file, which contains the 3D coordinates of a few user-defined landmarks that correspond to the 3D model described in the OBJ file.
An open access testing dataset can be downloaded from https://zenodo.org/record/3519103 and a User Manual and Algorithmic Description is also available at https://doi.org/10.5281/zenodo.3519248.

To compile from source enter:

    $ g++ skullanalyzer.cpp -o skullanalyzer
  
To see the help file enter:

    $ skullanalyzer --help

The supplementary GNU Octave function 'plot_features.m' can be used for visualizing the extracted geometric features. See the User Manual for more details.
