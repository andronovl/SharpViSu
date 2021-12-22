# SharpViSu
Software for correction, evaluation and visualization of single-molecule localization microscopy (SMLM) data. 
<a href="https://github.com/andronovl/SharpViSu/blob/master/SharpViSu%20user%20manual.pdf">  SharpViSu user manual  </a>

<b>   SharpViSu 2 </b> includes:

<b> • SplitViSu </b> for spectral demixing of SMLM data, acquired with a dichroic image splitter: <a href="https://github.com/andronovl/SharpViSu/blob/master/SplitViSu%20manual.pdf">  SplitViSu user manual  </a>

<b> • ClusterViSu </b> for cluster analysis of SMLM data using Voronoi diagrams and Ripley functions: <a href="https://github.com/andronovl/SharpViSu/blob/master/ClusterViSu%20user%20manual.pdf">  ClusterViSu user manual  </a>

<br>

## Installation 

The software can be used either as a standalone application for Windows or as an application for MATLAB 

<br>

<b><i> Standalone SharpViSu with SplitViSu & ClusterViSu for Windows:  </i></b>

Download and run <a href=https://github.com/andronovl/SharpViSu/tree/master/Installer> “Installer/SharpViSu_web_installer.exe” </a> and follow instructions. If not already installed, the MATLAB Compiler Runtime (mcr) will be downloaded from the web and installed automatically.

<br>

<b><i> Standalone SplitViSu or ClusterViSu for Windows:  </i></b>

Download and run <a href=https://github.com/andronovl/SharpViSu/tree/master/Installer> “Installer/SplitViSu_web_installer.exe or ”Installer/SharpViSu_web_installer.exe" </a> and follow instructions. If not already installed, the MATLAB Compiler Runtime (mcr) will be downloaded from the web and installed automatically.

<br>

<b><i> MATLAB application & source code (works from beforehand installed MATLAB environment):  </i></b>

Add folder <a href=https://github.com/andronovl/SharpViSu/tree/master/SharpViSu> “SharpViSu” </a> to the MATLAB search path. SharpViSu 2 was developed in MATLAB R2021b, but should work also in older versions. Toolboxes required: Image Processing Toolbox, Signal Processing Toolbox, Statistics and Machine Learning Toolbox

<br>

## Usage 

<a href="https://github.com/andronovl/SharpViSu/blob/master/SharpViSu%20user%20manual.pdf"> <b> SharpViSu User manual </b> </a>

<a href="https://github.com/andronovl/SharpViSu/blob/master/SplitViSu%20manual.pdf"> <b> User manual for spectral demixing of SMLM data in SplitViSu </b> </a>

<a href="https://github.com/andronovl/SharpViSu/tree/master/Data/SplitViSu"> <b> Test data for spectral demixing in SplitViSu </b> </a>

<br>

## References 
Andronov, L., Lutz, Y., Vonesch, J.-L. & Klaholz, B. P. SharpViSu: integrated analysis and segmentation of super-resolution microscopy data. <i>Bioinformatics</i> (2016) <a href="http://bioinformatics.oxfordjournals.org/content/early/2016/03/17/bioinformatics.btw123">doi:10.1093/bioinformatics/btw123</a>

<br>
Andronov, L., Orlov, I., Lutz, Y., Vonesch, J.-L. & Klaholz, B. P. ClusterViSu, a method for clustering of protein complexes by Voronoi tessellation in super-resolution microscopy. <i>Scientific Reports</i> <b>6</b>, 24084 (2016) <a href="http://www.nature.com/articles/srep24084">http://www.nature.com/articles/srep24084</a>
