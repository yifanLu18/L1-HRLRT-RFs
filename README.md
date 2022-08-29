# L1-HRLRT-RFs

removing surface waves in RFs by using high-resolution Radon transform

more details can be found in https://doi.org/10.1093/gji/ggac260

=======================================================

Thanks to the suggestion from ChenXin, the sub-function tauptime has been added. The program is from MatTaup: A TauP toolkit for MATLAB, which is written and maintained by Qin Li at the University of Washington. 

In the windows operating system, the general installation process is divided into four steps to ensure that the matlab script can be called smoothly:

Step 1: Unzip the compressed package matTaup.rar and put it in the path C:\Program Files\MATLAB\R2018a\java\jar (adjusted according to your matlab installation path);

Step 2: Execute the operation "edit classpath.txt" in matlab (if there is no modification permission, you need to change the permission of the file), add the path
"$matlabroot/java/jar/matTaup/matTaup/lib/matTaup.jar" under the last line";

Step 3: Add the folder to the matlab path to ensure that matlab can successfully find the *.m script;

Step 4: Restart matlab, execute the example "taupTime([],50,'P,S','deg',45.6)", if no error is reported, the installation is complete
