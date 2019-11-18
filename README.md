# TRA
This code of Time-series-based Reflectance Adjustment (TRA) approach (Windows C/C++ version) is based 
on Window 10 and Microsoft Visual Studio 2017 (or higher). It requires the installation of HDF library
(version 4.2.11) for reading the HLS product (HDF format).<br> <br>  

When you open ***TRA.sln***, you need to check configuration. Choose ***Release*** and ***x64***, and then make sure:<br>

  (1) Use ***MFC in a Shared DLL***;<br> 

  (2) Use ***Multi-Bytes Character Set***;<br> 

  (3) Add the ***include and library directories*** of HDF;<br> 

  (4) Add ***_CTR_SECURE_NO_WARMINGs*** into preprocessor.<br> 

More details about how to use this code are given in the [How to use.PDF](https://github.com/GERSL/TRA/blob/master/How%20to%20use.pdf).<br><br>  

There are ***two functions for running TRA***: one is for a ***single pixel*** and the other is for the ***MGRS tile***.<br>

(1) When running a single pixel, you need to define the ***path of L30 and S30 point*** data (CSV format) and the ***output path*** of fitted parameters (TXT format). The CSV file should be organized as nine columns (***Year, DOY, Blue, Green, Red, NIR, SWIR1, SWIR2, and ClearFlag***). The ClearFlag labeled as 1 means the clear-sky observation.<br>

(2) When running a MGRS tile, you need to create the TXT file for ***all paths of the HLS L30 or S30*** files (use the command: ***dir *S30*.hdf /s/b >pathS30.txt***), and then define the paths of two TXT files and the output path of the HDF file with the fitted parameters.<br> 

-------------
If using this code, please cite the following paper:

Shang, R., Zhu, Z., 2019. [Harmonizing Landsat 8 and Sentinel-2: A time-series-based reflectance adjustment approach](https://www.sciencedirect.com/science/article/pii/S0034425719304584). *Remote Sensing of Environment*. 235, 111439.
