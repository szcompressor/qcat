# Quick Compression Analysis Toolkit (QCAT)

QCAT is a lightweight tool to analyze the data in the context of lossy compression. You can use it to change the size of data file in binary, convert binary data files to texture files and vice versa. You can also plot the raw data file and decompressed data file and visualize the difference by using the executable 'PlotSliceImage' for 2D and 3D datasets in different directions. 

## Installation
./configure --prefix=[install path]; (e.g., ./configure --prefix=/home/sdi/qcat-1.3-install)

make;

make install 

(After compiling the package, you can use the executables in bin/ to do the compression analysis.)

## Quick Start

```bash
[sdi@sdihost qcat-0.1]$ cd qcat-1.3-install/bin
[sdi@sdihost bin]$ ls
changeDataSize  convertBytesToTxtDouble  convertDataToLogDouble  convertFloatToDouble     convertTxtToBytesFloat  generateRandomData   PlotSliceImage	predCR
compareData     convertBytesToTxtFloat   convertDoubleToFloat    convertTxtToBytesDouble  generateIndexData       generateRandomData2  splitComplexData printProperty
calculateSSIM

[sdi@sdihost bin]$ printProperty 
Usage: printDataProperty [dataType] tgtFilePath]
Example: printDataProperty -f testfloat_8_8_128.dat
 
[sdi@sdihost bin]$ ./PlotSliceImage 
Usage: PlotSliceImage <options>
Options:
*DataType:
	-f : single precision
	-d : double precision
* input data file:
	-i <original data file> : specify original data file
	-x <decompressed data file> : specify decompressed data file
* dimensions: 
	-2 <nx> <ny> : dimensions for 2D data such as data[ny][nx]
	-3 <nx> <ny> <nz> : dimensions for 3D data such as data[nz][ny][nx] 
* Operation: 
	-m <mode> : several operations as follows.
		DIFF : plot the difference between original data (specified by -i) and decompressed data (specified by -x)
		INDV : plot the individual dataset (specified using -i)
	-n <domain>: domain space to be plotted.
		ORI : original data domain
		LOG : log_10 domain of the data
	-p <dimension> : along which dimension for plotting the slice. (options: 1, 2 or 3; default setting: 3
	-s <slice number>: slice number if input is a 3D dataset
* Output: 
	-o <output image file> : Specify the output image file (.png format)
Examples:
	PlotSliceImage -f -i /root/usecases/CESM/CLDHGH_1_1800_3600.dat -2 3600 1800 -m INDV -n ORI -o /root/usecases/CESM/CLDHGH_1_1800_3600.png
	PlotSliceImage -f -i /root/usecases/Hurricane/QSNOWf48.dat -3 500 500 100 -m INDV -n ORI -p 3 -s 50 -o /root/usecases/Hurricane/QSNOWf48.png
	PlotSliceImage -f -i /root/usecases/Hurricane/Uf48.dat -3 500 500 100 -m INDV -n ORI -s 50 -o /root/usecases/Hurricane/Uf48.png
	PlotSliceImage -f -i /root/usecases/CESM/CLDHGH_1_1800_3600.dat -x /root/usecases/CESM/CLDHGH_1_1800_3600.dat.sz.out -2 3600 1800 -m DIFF -n ORI -o /root/usecases/CESM/CLDHGH_1_1800_3600-diff.png

[sdi@sdihost bin]$ ./compareData
Usage: compareData [datatype (-f or -d)] [original data file] [decompressed data file]
			-f means single precision; -d means double precision
Example: compareData -f testfloat_8_8_128.dat testfloat_8_8_128.dat.sz.out

[sdi@sdihost bin]$ ./convertTxtToBytesDouble 
Usage: convertTxtToBytesDouble [srcFilePath] [tgtFilePath]
Example: convertTxtToBytesDouble testfloat_8_8_128.txt testfloat_8_8_128.dat

[sdi@sdihost bin]$ ./convertTxtToBytesFloat
Usage: convertTxtToBytesFloat [srcFilePath] [tgtFilePath]
Example: convertTxtToBytesFloat testfloat_8_8_128.txt testfloat_8_8_128.dat

[sdi@sdihost bin]$ ./convertBytesToTxtDouble 
Usage: convertBytesToTxtDouble [srcFilePath] newsize [tgtFilePath]
Example: convertBytesToTxtDouble testfloat_8_8_128.dat 200 testfloat_8_8_128.xls

[sdi@sdihost bin]$ ./convertBytesToTxtFloat 
Usage: convertBytesToTxtFloat [srcFilePath] newsize [tgtFilePath]
Example: convertBytesToTxtFloat testfloat_8_8_128.dat 200 testfloat_8_8_128.xls
Example: convertBytesToTxtFloat testfloat_8_8_128.dat -1 testfloat_8_8_128.xls

[sdi@sdihost bin]$ ./changeDataSize
Usage: changeDataSize [srcFilePath] newsize
Example: changeDataSize testfloat_8_8_128.dat 200

[sdi@sdihost bin]$ ./convertFloatToDouble 
Usage: convertFloatToDouble [srcFilePath] [tgtFilePath]
Example: convertFloatToDouble testfloat_8_8_128.dat testdouble_8_8_128.f64 

[sdi@sdihost bin]$ ./convertDoubleToFloat 
Usage: convertDoubleToFloat [srcFilePath] [tgtFilePath]
Example: convertDoubleToFloat testdouble_8_8_128.dat testfloat_8_8_128.dat

[sdi@sdihost bin]$ ./calculateSSIM 
Usage: calculateSSIM [datatype (-f or -d)] [original data file] [decompressed data file] [dimesions... (from fast to slow)]
			-f means single precision; -d means double precision
Example: calculateSSIM -f CLOUD_100x500x500.dat CLOUD_100x500x500.dat.sz.out 500 500 100

[sdi@sdihost bin]$ ./computePDF
Usage: computePDF [datatype (-f or -d)] [original data file] [decompressed data file]
			-f means single precision; -d means double precision
Example: computePDF -f CLOUD_100x500x500.dat CLOUD_100x500x500.dat.sz.out

[sdi@sdihost bin]$ ./predCR
Usage: predCR [datatype (-f or -d)] [quantBinCapacity] [errorBound] [original data file] [predcted data file]
			-f means single precision; -d means double precision
Example: predCR -f 1024 1E-1 original.dat predicted.dat

```

Summary of particularly useful executables: 
1. printProperty : print the property of the dataset (such as value range and entropy)
2. PlotSliceImage : It helps you to plot the images based on a 2D or 3D scientific data file (stored in binary format such as little endian type).
3. compareData: compare two data files (they are generally expected to be original data file and decompressed data file, respectively) and show the compression related metrics such as psnr.
4. calculateSSIM: calculate the visualization SSIM for the multi-dimensional dataset (original data vs. decompressed data)
5. computePDF: compute the probaility density function (PDF) of compression errors between original data file and decompressed data file. 
6. predCR: If you have the original data file and a prediction data file, then predCR can help you generate the compression ratio in terms of the SZ error-bounded prediction framework. 



Contact: sdi@anl.gov
