# SharpViSu
Software for correction, evaluation and visualization of single-molecule localization microscopy (SMLM) data. 
<a href="https://github.com/andronovl/SharpViSu/blob/master/SharpViSu%20user%20manual.pdf">  SharpViSu user manual  </a>

<b>   SharpViSu 2 </b> also includes:

<b> • SplitViSu </b> for spectral demixing of SMLM data, acquired with a dichroic image splitter: <a href="https://github.com/andronovl/SharpViSu/blob/master/SplitViSu%20manual.pdf">  SplitViSu user manual  </a>

<b> • ClusterViSu </b> for cluster analysis of SMLM data using Voronoi diagrams and Ripley functions: <a href="https://github.com/andronovl/SharpViSu/blob/master/ClusterViSu%20user%20manual.pdf">  ClusterViSu user manual  </a>

<br>

## Installation 

The software can be used either as a standalone application for Windows or as an application for MATLAB 

<br>

<b><i> Standalone SharpViSu with SplitViSu & ClusterViSu for Windows (recommended):  </i></b>

Download and run the latest version of <a href=https://github.com/andronovl/SharpViSu/tree/master/Installer> “Installer/SharpViSu_web_installer.exe” </a> and follow instructions. If not already installed, the MATLAB Compiler Runtime (mcr) will be downloaded from the web and installed automatically.

<br>

<b><i> Standalone SplitViSu or ClusterViSu for Windows:  </i></b>

Download and run the latest version of  <a href=https://github.com/andronovl/SharpViSu/tree/master/Installer> “Installer/SplitViSu_web_installer.exe or ”Installer/ClusterViSu_web_installer.exe" </a> and follow instructions. If not already installed, the MATLAB Compiler Runtime (mcr) will be downloaded from the web and installed automatically.

<br>

<b><i> MATLAB application & source code (works under your MATLAB environment):  </i></b>

Add folder <a href=https://github.com/andronovl/SharpViSu/tree/master/SharpViSu> “SharpViSu” </a> to your MATLAB search path. SharpViSu v2 was developed in MATLAB R2021b, but should work also in older versions. Toolboxes required: Image Processing Toolbox, Signal Processing Toolbox, Statistics and Machine Learning Toolbox

<br>

## User manuals 

<a href="https://github.com/andronovl/SharpViSu/blob/master/SharpViSu%20user%20manual.pdf"> <b> SharpViSu User manual </b> </a>

<a href="https://github.com/andronovl/SharpViSu/blob/master/SplitViSu%20manual.pdf"> <b> User manual for spectral demixing of SMLM data in SplitViSu </b> </a>

<a href="https://github.com/andronovl/SharpViSu/tree/master/Data/SplitViSu"> <b> Test data for spectral demixing in SplitViSu </b> </a>

<br>

## References 
Andronov, L., Lutz, Y., Vonesch, J.-L. & Klaholz, B. P. SharpViSu: integrated analysis and segmentation of super-resolution microscopy data. <i>Bioinformatics</i> (2016) <a href="http://bioinformatics.oxfordjournals.org/content/early/2016/03/17/bioinformatics.btw123">doi:10.1093/bioinformatics/btw123</a>

<br>

Andronov, L., Orlov, I., Lutz, Y., Vonesch, J.-L. & Klaholz, B. P. ClusterViSu, a method for clustering of protein complexes by Voronoi tessellation in super-resolution microscopy. <i>Scientific Reports</i> <b>6</b>, 24084 (2016) <a href="http://www.nature.com/articles/srep24084">http://www.nature.com/articles/srep24084</a>

<br>

Andronov L., Michalon J., Ouararhni K., Orlov I., Hamiche A., Vonesch J-L, Klaholz B.P. 3DClusterViSu: 3D clustering analysis of super-resolution microscopy data by 3D Voronoi tessellations. <i>Bioinformatics</i> <b>34</b> (2018) <a href=https://academic.oup.com/bioinformatics/article/34/17/3004/4960045> https://doi.org/10.1093/bioinformatics/bty200 </a>

<br>

Andronov, L., Genthial, R., Hentsch, D., & Klaholz, B. P. A spectral demixing method for high-precision multi-color localization microscopy. <i>bioRxiv</i> (2021) <a href="https://www.biorxiv.org/content/10.1101/2021.12.23.473862v1.full">https://doi.org/10.1101/2021.12.23.473862</a>

<br>

#### Examples of use
Andronov, L., Ouararhni, K., Stoll, I., Klaholz, B.P., Hamiche, A. CENP-A nucleosome clusters form rosette-like structures around HJURP during G1. <i> Nat. Commun. </i> <b>10</b>, 4436 (2019) <a href="https://www.nature.com/articles/s41467-019-12383-3">https://doi.org/10.1038/s41467-019-12383-3</a>

<br>

Lemaître, C., Grabarz, A., Tsouroula, K., Andronov, L., Furst, A., Pankotai, T., Heyer, V., Rogier, M., Attwood, K M., Kessler, P., Dellaire, G., Klaholz, B., Reina-San-Martin, B., Soutoglou, E. Nuclear position dictates DNA repair pathway choice. <i>Genes. Dev.</i> <b>28</b>, 2450-2463 (2014) <a href="http://genesdev.cshlp.org/content/28/22/2450">https://doi.org/10.1101/gad.248369.114</a>

<br>

Andronov L., Vonesch J-L, Klaholz B.P. Practical Aspects of Super-Resolution Imaging and Segmentation of Macromolecular Complexes by dSTORM. <i> In: Poterszman A. (eds) Multiprotein Complexes. Methods in Molecular Biology, </i> <b>2247</b> 271-286 (2021) Humana, New York, NY. <a href="https://link.springer.com/protocol/10.1007/978-1-0716-1126-5_15">https://doi.org/10.1007/978-1-0716-1126-5_15</a>

<br>

Andronov, L., Genthial, R., Hentsch, D., & Klaholz, B. P. A spectral demixing method for high-precision multi-color localization microscopy. <i>bioRxiv</i> (2021) <a href="https://www.biorxiv.org/content/10.1101/2021.12.23.473862v1.full">https://doi.org/10.1101/2021.12.23.473862</a>
