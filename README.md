ChildBrainFlow Pipeline
=======================

ChildBrainFlow is an end-to-end pipeline that performs tractography, t1 reconstruction and connectomics.
It is essentially a merged version of multiple individual pipeline to avoid the handling of inputs/outputs
between flows with some parameters tuned for pediatric brain scans. Here is a list of flows from which
process have been taken:

    1. TractoFlow (https://github.com/scilus/tractoflow.git) [1]
    2. FreeSurfer-Flow (https://github.com/scilus/freesurfer_flow)
    3. Connectoflow (https://github.com/scilus/connectoflow)

> [!NOTE]
> Please note that some steps have been removed from the original pipelines if they were not relevant for pediatric data. If you need some of these steps, please use the original pipelines.

Nextflow
--------
To install nextflow, please see : https://www.nextflow.io/docs/latest/getstarted.html#requirements 

The pipeline export by default a `` parameters.json `` within the output directory to provide a documentation of the parameters used during the execution. For a more detailed report (excluding execution's parameters), the default feature of nextflow `` -with-report <report_name.html> `` can be used to export a html report. Simply had this your command line when launching the pipeline: 

```
nextflow run main.nf --input <input> <other_arguments> -with-report <report_name.html>
```

Apptainer
---------
If you are running this pipeline on Linux, it is recommended to run the pipeline using an apptainer image. 
The pipeline comes with a recipe file (`` /containers/apptainer_recipe.def ``) containing all the required 
dependencies to successfully run every profiles. To build the apptainer image, run this command: 

```
sudo apptainer build <image_name> </path/to/apptainer_recipe.def
```

Docker
------
If you are on MacOS or Windows, you can use Docker to run ChildBrainFlow. The pipeline comes with
a Dockerfile containing all the dependencies required to successfully run every profiles of the pipeline. 
Simply run this command from inside the  `` /containers/ `` folder:

```
docker build -t <container_name> .
```
> [!WARNING]
> Due to the high number of dependencies (ANTs, FSL, MRtrix3, Scilpy, Freesurfer, FastSurfer, etc.), the resulting docker image can be pretty large (~ 40Gb).

Usage
-----
See _USAGE_ or run `` nextflow run main.nf --help `` for more details.

References
----------
If you used this pipeline, please cite :

[1] Theaud, G., Houde, J.-C., Boré, A., Rheault, F., Morency, F., Descoteaux, M.,
        TractoFlow: A robust, efficient and reproducible diffusion MRI pipeline
        leveraging Nextflow & Singularity, NeuroImage,
        https://doi.org/10.1016/j.neuroimage.2020.116889.

[2] Kurtzer GM, Sochat V, Bauer MW Singularity: Scientific containers for mobility of compute. PLoS ONE 12(5)
        (2017): e0177459. https://doi.org/10.1371/journal.pone.0177459

[3] P. Di Tommaso, et al. Nextflow enables reproducible computational workflows. Nature Biotechnology 35,
        316–319 (2017) https://doi.org/10.1038/nbt.3820

[4] Garyfallidis, E., Brett, M., Amirbekian, B., Rokem, A., Van Der Walt, S., Descoteaux, M., Nimmo-Smith, I.,
        2014. Dipy, a library for the analysis of diffusion mri data. Frontiers in neuroinformatics 8, 8.
        https://doi.org/10.3389/fninf.2014.00008

[5] Tournier, J. D., Smith, R. E., Raffelt, D. A., Tabbara, R., Dhollander, T., Pietsch, M., … & Connelly, A.
        (2019). MRtrix3: A fast, flexible and open software framework for medical image processing and visualisation.
        NeuroImage 202, https://doi.org/10.1016/j.neuroimage.2019.116137

[6] Avants, B. B., Tustison, N., & Song, G. (2009). Advanced normalization tools (ANTS). Insight j, 2(365), 1-35.

[7] Jenkinson, M., Beckmann, C.F., Behrens, T.E., Woolrich, M.W., Smith, S.M., 2012. Fsl. Neuroimage 62,
        782–790. https://doi.org/10.1016/j.neuroimage.2011.09.015

[8] Fischl, B., Salat, D.H., Busa, E., Albert, M., Dieterich, M., Haselgrove, C., van der Kouwe, A., Killiany, 
        R., Kennedy, D., Klaveness, S., Montillo, A., Makris, N., Rosen, B., Dale, A.M., 2002. Whole brain 
        segmentation: automated labeling of neuroanatomical structures in the human brain. Neuron 33, 341-355.
        https://doi.org/10.1016/s0896-6273(02)00569-x
