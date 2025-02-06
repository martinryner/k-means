# k-means
Global k-means algorithm using cutting planes

This code was run on MATLAB R2022b on Windows 11.

Usage
1. Add the src and data folder and its subfolders to the Matlab path.
2. Make sure the "parallel computing toolbox" package is installed in MATLAB (or change "parfor" to "for" in the code)
3. Install Gurobi and configure for Matlab
3. Compile the mex files to your platform in src/mexemd
4. Compile the mex files to your patform in src/OpenCL according to instructions in these packages. 
   Alternatively, you can set a global variable forceCPU = true to run on CPU, but this is slower.
4. Run the tests in kmeantest.m
