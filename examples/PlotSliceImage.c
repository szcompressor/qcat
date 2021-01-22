#include <stdio.h> 
#include <stdlib.h>
#include <string.h> 
#include "rw.h"
#include "qcat_gnuplot.h"

void usage()
{
	printf("Usage: PlotSliceImage <options>\n");
	printf("Options:\n");
	printf("*DataType:\n");
	printf("	-f : single precision\n");
	printf("	-d : double precision\n");
	printf("* input data file:\n");
	printf("	-i <original data file> : specify original data file\n");
	printf("	-x <decompressed data file> : specify decompressed data file\n");
	printf("* dimensions: \n");
	printf("	-2 <nx> <ny> : dimensions for 2D data such as data[ny][nx]\n");
	printf("	-3 <nx> <ny> <nz> : dimensions for 3D data such as data[nz][ny][nx] \n");
	printf("* Operation: \n");
	printf("	-m <mode> : several operations as follows.\n");
	printf("		DIFF : plot the difference between original data (specified by -i) and decompressed data (specified by -x)\n");
	printf("		INDV : plot the individual dataset (specified using -i)\n");
	printf("	-n <domain>: domain space to be plotted.\n");
	printf("		ORI : original data domain\n");
	printf("		LOG : log_10 domain of the data\n");
	printf("	-p <dimension> : along which dimension for plotting the slice. (options: 1, 2 or 3; default setting: 3\n");
	printf("	-s <slice number>: slice number if input is a 3D dataset\n");
	printf("	-r <min> <max> : specify the value range (i.e., min and max) for the plotting\n");

	printf("* Output: \n");
	printf("	-o <output image file> : Specify the output image file (.png format)\n");
	printf("Examples:\n");
	printf("	PlotSliceImage -f -i /root/usecases/CESM/CLDHGH_1_1800_3600.dat -2 3600 1800 -m INDV -n ORI -o /root/usecases/CESM/CLDHGH_1_1800_3600.png\n");
	printf("	PlotSliceImage -f -i /root/usecases/Hurricane/QSNOWf48.dat -3 500 500 100 -m INDV -n ORI -p 3 -s 50 -o /root/usecases/Hurricane/QSNOWf48.png\n");
	printf("	PlotSliceImage -f -i /root/usecases/Hurricane/Uf48.dat -3 500 500 100 -m INDV -n ORI -s 50 -o /root/usecases/Hurricane/Uf48.png\n");
	printf("	PlotSliceImage -f -i /root/usecases/CESM/CLDHGH_1_1800_3600.dat -x /root/usecases/CESM/CLDHGH_1_1800_3600.dat.sz.out -2 3600 1800 -m DIFF -n ORI -o /root/usecases/CESM/CLDHGH_1_1800_3600-diff.png\n");
	
}

int main(int argc, char * argv[])
{
	int mode = 0; //0: plot original data ; 1: plot difference
	int domain = 0; //0: original domain ; 1: log domain
	size_t sliceNumber = 0;
	int dataType = 0;
	int status = 0;
	char *oriFilePath = NULL, *decFilePath = NULL, *outputFilePath = NULL;
	char gnuDataFilePath[640], gnuScriptFilePath[640];
	size_t r3 = 0, r2 = 0, r1 = 0;
	int plotDim = 3;
	float range_min = 0, range_max = 0;
	
	if(argc==1)
	{
		usage();
		return 0;
	}

	int i = 0;
	for(i=1;i<argc;i++)
	{
		if (argv[i][0] != '-' || argv[i][2])
			usage();
		switch (argv[i][1])
		{
		case '2':
			if (++i == argc || sscanf(argv[i], "%zu", &r1) != 1 ||
				++i == argc || sscanf(argv[i], "%zu", &r2) != 1)
				usage();
			break;
		case '3':
			if (++i == argc || sscanf(argv[i], "%zu", &r1) != 1 ||
				++i == argc || sscanf(argv[i], "%zu", &r2) != 1 ||
				++i == argc || sscanf(argv[i], "%zu", &r3) != 1)
				usage();		
			break;			
		case 'm':
			if (++i == argc)
				usage();
			if(strcmp(argv[i],"DIFF")==0)
				mode = 1;
			else if(strcmp(argv[i], "INDV")==0)
				mode = 0;
			else
			{
				printf("Error: wrong mode - please use either DIFF or INDV.\n");
				usage();
			}
			break;
		case 'p':
			if (++i == argc || sscanf(argv[i], "%d", &plotDim) != 1)
				usage();
			break;
		case 'r':
			if (++i == argc || sscanf(argv[i], "%f", &range_min) != 1 ||
				++i == argc || sscanf(argv[i], "%f", &range_max) != 1)
				usage();
			break;						
		case 'i':
			if (++i == argc)
				usage();
			oriFilePath = argv[i];
			break;
		case 'o':
			if (++i == argc)
				usage();
			outputFilePath = argv[i];	
			break;
		case 'x':
			if (++i == argc)
				usage();
			decFilePath = argv[i];
			break;
		case 'f': 
			dataType = 0;
			break;
		case 'd':
			dataType = 1;	
			break;
		case 's': 
			if (++i == argc)
				usage();
			sliceNumber = atoi(argv[i]);
			break;
		case 'n':
			if (++i == argc)
				usage();
			if(strcmp(argv[i],"ORI")==0)
				domain = 0;
			else if(strcmp(argv[i], "LOG")==0)
				domain = 1;
			else
			{
				printf("Error: wrong domain - please use either ORI or LOG.\n");
				usage();
			}
			break;		
		default: 
			usage();
			break;
		}
	}

	if(checkFileExistance(oriFilePath)==0)
	{
		printf("Error: the input file %s doesn't exist.\n", oriFilePath);
		exit(0);
	}
		
	if(r3!=0)
	{
		if(plotDim==3 && sliceNumber>=r3)
		{
			printf("Error: wrong sliceNumber: %zu. It should be in [0,%zu]\n", sliceNumber, r3-1);
			exit(0);			
		}
		if(plotDim==2 && sliceNumber>=r2)
		{
			printf("Error: wrong sliceNumber: %zu. It should be in [0,%zu]\n", sliceNumber, r2-1);
			exit(0);						
		}
		if(plotDim==1 && sliceNumber>=r1)
		{
			printf("Error: wrong sliceNumber: %zu. It should be in [0,%zu]\n", sliceNumber, r1-1);
			exit(0);									
		}
	}	
		
	size_t nbEle = 0;	
	
	char *dataDir = extractDirFromPath(oriFilePath);
	char tmpDir[1000], cmd[1000];
	sprintf(tmpDir, "%s/.tmp", dataDir);
	sprintf(cmd, "mkdir -p %s", tmpDir);
	system(cmd);
	
	if(dataType==0)
	{
		float* oriData = readFloatData(oriFilePath, &nbEle, &status);
		float* sliceData = NULL;
		if(mode==0) //plot original data
		{	
			sliceData = generateSliceData_float(plotDim, sliceNumber, r3, r2, r1, oriData, domain); 
			sprintf(gnuDataFilePath, "%s.orid", oriFilePath);
			sprintf(gnuScriptFilePath, "%s.orip", oriFilePath);
		}
		else if(mode==1) //plot difference
		{
			if(checkFileExistance(decFilePath)==0 )
			{
				printf("Error: the decompressed data file doesn't exist.\n");
				exit(0);
			}	

			float* decData = readFloatData(decFilePath, &nbEle, &status);
			sliceData = generateSliceDiff_float(plotDim, sliceNumber, r3, r2, r1, oriData, decData, domain);
			free(decData);	
			sprintf(gnuDataFilePath, "%s.difd", oriFilePath);
			sprintf(gnuScriptFilePath, "%s.difp", oriFilePath);			
		}

		switch(plotDim)
		{
			case 3:
				SDA_writeData_genuplotImage(sliceData, dataType, r2, r1, gnuDataFilePath);
				break;
			case 2:
				SDA_writeData_genuplotImage(sliceData, dataType, r3, r1, gnuDataFilePath);
				break;
			case 1:
				SDA_writeData_genuplotImage(sliceData, dataType, r3, r2, gnuDataFilePath);
				break;
		}
			
		free(sliceData);

		char** sliceImageStrs = NULL;
		switch(plotDim)
		{
			case 3:
				sliceImageStrs = genGnuplotScript_sliceImage(gnuDataFilePath, r2, r1, outputFilePath, range_min, range_max);
				break;
			case 2:
				sliceImageStrs = genGnuplotScript_sliceImage(gnuDataFilePath, r3, r1, outputFilePath, range_min, range_max);
				break;
			case 1:
				sliceImageStrs = genGnuplotScript_sliceImage(gnuDataFilePath, r3, r2, outputFilePath, range_min, range_max);
				break;
		}
		
		SDA_writeStrings(10, sliceImageStrs, gnuScriptFilePath);
		for(i=0;i<10;i++)
			free(sliceImageStrs[i]);
		free(sliceImageStrs);
		free(oriData);
	}
	else if(dataType==1)
	{
		double* oriData = readDoubleData(oriFilePath, &nbEle, &status);
		double* sliceData = NULL;
		if(mode==0) //plot original data
		{	
			sliceData = generateSliceData_double(plotDim, sliceNumber, r3, r2, r1, oriData, domain); 
			sprintf(gnuDataFilePath, "%s.orid", oriFilePath);
			sprintf(gnuScriptFilePath, "%s.orip", oriFilePath);			
		}
		else if(mode==1) //plot difference
		{
			if(checkFileExistance(decFilePath)==0 )
			{
				printf("Error: the decompressed data file doesn't exist.\n");
				exit(0);
			}	

			double* decData = readDoubleData(decFilePath, &nbEle, &status);
			sliceData = generateSliceDiff_double(plotDim, sliceNumber, r3, r2, r1, oriData, decData, domain);		
			free(decData);
			sprintf(gnuDataFilePath, "%s.difd", oriFilePath);
			sprintf(gnuScriptFilePath, "%s.difp", oriFilePath);			
		}

		switch(plotDim)
		{
			case 3:
				SDA_writeData_genuplotImage(sliceData, dataType, r2, r1, gnuDataFilePath);
				break;
			case 2:
				SDA_writeData_genuplotImage(sliceData, dataType, r3, r1, gnuDataFilePath);
				break;
			case 1:
				SDA_writeData_genuplotImage(sliceData, dataType, r3, r2, gnuDataFilePath);
				break;
		}
			
		free(sliceData);

		char** sliceImageStrs = NULL;
		switch(plotDim)
		{
			case 3:
				sliceImageStrs = genGnuplotScript_sliceImage(gnuDataFilePath, r2, r1, outputFilePath, range_min, range_max);
				break;
			case 2:
				sliceImageStrs = genGnuplotScript_sliceImage(gnuDataFilePath, r3, r1, outputFilePath, range_min, range_max);
				break;
			case 1:
				sliceImageStrs = genGnuplotScript_sliceImage(gnuDataFilePath, r3, r2, outputFilePath, range_min, range_max);
				break;
		}

		SDA_writeStrings(10, sliceImageStrs, gnuScriptFilePath);
		for(i=0;i<10;i++)
			free(sliceImageStrs[i]);
		free(sliceImageStrs);
		free(oriData);
	}	
	
	sprintf(cmd, "cd %s;gnuplot %s;mv %s %s/.tmp; mv %s %s/.tmp", dataDir, gnuScriptFilePath, gnuScriptFilePath, dataDir, gnuDataFilePath, dataDir);
	printf("Image file is plotted and put here: %s\n", outputFilePath);
	system(cmd);
}
