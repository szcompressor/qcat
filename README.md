# Quick Compression Analysis (QCA) 0.1

QCA is a lightweight tool to analyze the data in the context of lossy compression. You can use it to change the size of data file in binary, convert binary data files to texture files and vice versa. You can also plot the raw data file and decompressed data file and visualize the difference by using the executable 'PlotSliceImage' for 2D and 3D datasets in different directions. 

## Installation
./configure --prefix=`pwd`/qca-0.1-install
make
make install 

## Test

[sdi@sdihost sda-0.1]$ cd sda-0.1-install/bin
[sdi@sdihost bin]$ ls
changeDataSize           convertBytesToTxtFloat  convertTxtToBytesDouble  generateRandomData   PlotSliceImage
convertBytesToTxtDouble  convertDataToLogDouble  convertTxtToBytesFloat   generateRandomData2  splitComplexData
[sdi@sdihost bin]$ PlotSliceImage 
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

[sdi@sdihost examples]$ compareData
Test case: compareData [datatype (-f or -d)] [original data file] [decompressed data file]
			-f means single precision; -d means double precision
Example: compareData -f testfloat_8_8_128.dat testfloat_8_8_128.dat.sz.out

[sdi@sdihost bin]$ convertTxtToBytesDouble 
Test case: convertTxtToBytesDouble [srcFilePath] [tgtFilePath]
Example: convertTxtToBytesDouble testfloat_8_8_128.txt testfloat_8_8_128.dat
[sdi@sdihost bin]$ convertTxtToBytesFloat
Test case: convertTxtToBytesFloat [srcFilePath] [tgtFilePath]
Example: convertTxtToBytesFloat testfloat_8_8_128.txt testfloat_8_8_128.dat
[sdi@sdihost bin]$ convertBytesToTxtDouble 
Test case: convertBytesToTxtDouble [srcFilePath] newsize [tgtFilePath]
Example: convertBytesToTxtDouble testfloat_8_8_128.dat 200 testfloat_8_8_128.xls
[sdi@sdihost bin]$ convertBytesToTxtFloat 
Test case: convertBytesToTxtFloat [srcFilePath] newsize [tgtFilePath]
Example: convertBytesToTxtFloat testfloat_8_8_128.dat 200 testfloat_8_8_128.xls
Example: convertBytesToTxtFloat testfloat_8_8_128.dat -1 testfloat_8_8_128.xls
[sdi@sdihost bin]$ changeDataSize
Test case: changeDataSize [srcFilePath] newsize
Example: changeDataSize testfloat_8_8_128.dat 200
[sdi@sdihost examples]$ convertFloatToDouble 
Test case: convertFloatToDouble [srcFilePath] [tgtFilePath]
Example: convertFloatToDouble testfloat_8_8_128.dat testdouble_8_8_128.f64 
[sdi@sdihost examples]$ convertDoubleToFloat 
Test case: convertDoubleToFloat [srcFilePath] [tgtFilePath]
Example: convertDoubleToFloat testdouble_8_8_128.dat testfloat_8_8_128.dat

Contact: sdi@anl.gov
