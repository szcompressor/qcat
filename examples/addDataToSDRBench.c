/**
 *  @file addDataToSDRBench.c
 *  @author Sheng Di
 *  @brief Analyze binary data files and generate SDRBench YAML dataset entry
 *  (C) 2024 by Mathematics and Computer Science (MCS), Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <dirent.h>
#include <sys/stat.h>
#include <float.h>
#include "ByteToolkit.h"
#include "rw.h"
#include "qcat_dataAnalysis.h"

#define MAX_FILES 1024
#define MAX_LINE 1024
#define MAX_PATH_LEN 1024
#define MAX_DIMS 5

typedef struct {
	char application_name[256];
	char source[256];
	char source_url[512];
	char source_note[512];
	char type[256];
	char data_type[32];
	size_t dims[MAX_DIMS];
	int ndims;
	char file_extension[32];
} SDRBenchConfig;

static void trim(char *str)
{
	int len = strlen(str);
	while (len > 0 && (str[len-1] == '\n' || str[len-1] == '\r' ||
	                    str[len-1] == ' '  || str[len-1] == '\t'))
		str[--len] = '\0';
	char *start = str;
	while (*start == ' ' || *start == '\t')
		start++;
	if (start != str)
		memmove(str, start, strlen(start) + 1);
}

static void parseConfig(const char *path, SDRBenchConfig *cfg)
{
	memset(cfg, 0, sizeof(SDRBenchConfig));

	FILE *f = fopen(path, "r");
	if (!f) {
		fprintf(stderr, "Error: cannot open config file %s\n", path);
		exit(1);
	}

	char line[MAX_LINE];
	while (fgets(line, sizeof(line), f)) {
		trim(line);
		if (line[0] == '#' || line[0] == '\0')
			continue;

		char *eq = strchr(line, '=');
		if (!eq)
			continue;

		*eq = '\0';
		char *key = line;
		char *value = eq + 1;
		trim(key);
		trim(value);

		if (strcmp(key, "ApplicationData") == 0)
			strncpy(cfg->application_name, value, sizeof(cfg->application_name) - 1);
		else if (strcmp(key, "source") == 0)
			strncpy(cfg->source, value, sizeof(cfg->source) - 1);
		else if (strcmp(key, "source_url") == 0)
			strncpy(cfg->source_url, value, sizeof(cfg->source_url) - 1);
		else if (strcmp(key, "source_note") == 0)
			strncpy(cfg->source_note, value, sizeof(cfg->source_note) - 1);
		else if (strcmp(key, "type") == 0)
			strncpy(cfg->type, value, sizeof(cfg->type) - 1);
		else if (strcmp(key, "data_type") == 0)
			strncpy(cfg->data_type, value, sizeof(cfg->data_type) - 1);
		else if (strcmp(key, "dimensions") == 0) {
			cfg->ndims = 0;
			char *tmp = strdup(value);
			char *tok = strtok(tmp, "x");
			while (tok && cfg->ndims < MAX_DIMS) {
				cfg->dims[cfg->ndims++] = (size_t)atol(tok);
				tok = strtok(NULL, "x");
			}
			free(tmp);
		}
		else if (strcmp(key, "file_extension") == 0)
			strncpy(cfg->file_extension, value, sizeof(cfg->file_extension) - 1);
	}
	fclose(f);

	if (cfg->file_extension[0] == '\0') {
		if (strcmp(cfg->data_type, "float32") == 0)
			strcpy(cfg->file_extension, "f32");
		else if (strcmp(cfg->data_type, "float64") == 0)
			strcpy(cfg->file_extension, "f64");
	}
}

static int hasExtension(const char *filename, const char *ext)
{
	if (ext[0] == '\0')
		return 1;
	const char *dot = strrchr(filename, '.');
	if (!dot)
		return 0;
	return strcmp(dot + 1, ext) == 0;
}

static int compareStrings(const void *a, const void *b)
{
	return strcmp(*(const char **)a, *(const char **)b);
}

static int isRegularFile(const char *path)
{
	struct stat st;
	if (stat(path, &st) != 0)
		return 0;
	return S_ISREG(st.st_mode);
}

int main(int argc, char *argv[])
{
	if (argc < 3) {
		printf("Usage: addDataToSDRBench <sdrbench.conf> <data_folder>\n");
		printf("\nAnalyzes binary data files and generates:\n");
		printf("  1. data-property.txt  - property output for each data file\n");
		printf("  2. dataset_entry.yml  - YAML entry for SDRBench datasets.yml\n");
		printf("\nExample: addDataToSDRBench sdrbench.conf ./CESM-data\n");
		exit(0);
	}

	const char *configPath = argv[1];
	const char *dataFolder = argv[2];

	SDRBenchConfig cfg;
	parseConfig(configPath, &cfg);

	if (cfg.application_name[0] == '\0') {
		fprintf(stderr, "Error: ApplicationData not specified in config\n");
		exit(1);
	}
	if (cfg.data_type[0] == '\0') {
		fprintf(stderr, "Error: data_type not specified in config\n");
		exit(1);
	}
	if (cfg.ndims == 0) {
		fprintf(stderr, "Error: dimensions not specified in config\n");
		exit(1);
	}

	int dataType = (strcmp(cfg.data_type, "float32") == 0) ? QCAT_FLOAT : QCAT_DOUBLE;

	int x = 1;
	char *y = (char *)&x;
	if (*y == 1)
		sysEndianType = 0;
	else
		sysEndianType = 1;

	DIR *dir = opendir(dataFolder);
	if (!dir) {
		fprintf(stderr, "Error: cannot open directory %s\n", dataFolder);
		exit(1);
	}

	char *fileNames[MAX_FILES];
	int nFiles = 0;
	struct dirent *entry;
	while ((entry = readdir(dir)) != NULL && nFiles < MAX_FILES) {
		if (hasExtension(entry->d_name, cfg.file_extension)) {
			char fullPath[MAX_PATH_LEN];
			snprintf(fullPath, sizeof(fullPath), "%s/%s", dataFolder, entry->d_name);
			if (isRegularFile(fullPath)) {
				fileNames[nFiles] = strdup(entry->d_name);
				nFiles++;
			}
		}
	}
	closedir(dir);

	if (nFiles == 0) {
		fprintf(stderr, "Error: no data files found in %s with extension .%s\n",
		        dataFolder, cfg.file_extension);
		exit(1);
	}

	qsort(fileNames, nFiles, sizeof(char *), compareStrings);
	printf("Found %d data files in %s\n", nFiles, dataFolder);

	double *entropy_8bit = (double *)malloc(nFiles * sizeof(double));
	double *entropy_fp = (double *)malloc(nFiles * sizeof(double));
	size_t totalBytes = 0;

	FILE *propFile = fopen("data-property.txt", "w");
	if (!propFile) {
		fprintf(stderr, "Error: cannot create data-property.txt\n");
		exit(1);
	}

	for (int i = 0; i < nFiles; i++) {
		char filePath[MAX_PATH_LEN];
		snprintf(filePath, sizeof(filePath), "%s/%s", dataFolder, fileNames[i]);

		printf("\nProcessing [%d/%d]: %s\n", i + 1, nFiles, fileNames[i]);

		int status = 0;
		size_t nbEle = 0;
		void *data = NULL;

		if (dataType == QCAT_FLOAT)
			data = (void *)readFloatData(filePath, &nbEle, &status);
		else
			data = (void *)readDoubleData(filePath, &nbEle, &status);

		if (status != RW_SCES) {
			fprintf(stderr, "Error: cannot read file %s\n", filePath);
			entropy_8bit[i] = 0;
			entropy_fp[i] = 0;
			free(data);
			continue;
		}

		QCAT_DataProperty *property = computeProperty(dataType, data, nbEle, 3);

		fprintf(propFile, "========== %s ==========\n", fileNames[i]);
		fprintf(propFile, "numOfElem = %zu\n", property->numOfElem);
		fprintf(propFile, "totalDataSize = %zu bytes (%f MB)\n",
		        property->totalByteSize,
		        property->totalByteSize / (1024.0 * 1024.0));
		fprintf(propFile, "min = %f\n", property->minValue);
		fprintf(propFile, "max = %f\n", property->maxValue);
		fprintf(propFile, "valueRange = %f\n", property->valueRange);
		fprintf(propFile, "avgValue = %f\n", property->avgValue);
		if (property->entropy_8bits >= 0)
			fprintf(propFile, "entropy_8bits = %f\n", property->entropy_8bits);
		if (property->entropy_32bits >= 0)
			fprintf(propFile, "entropy_32bits = %f\n", property->entropy_32bits);
		if (property->entropy_64bits >= 0)
			fprintf(propFile, "entropy_64bits = %f\n", property->entropy_64bits);
		fprintf(propFile, "zeromean_variance = %f\n\n", property->zeromean_variance);

		entropy_8bit[i] = property->entropy_8bits;
		if (dataType == QCAT_FLOAT)
			entropy_fp[i] = property->entropy_32bits;
		else
			entropy_fp[i] = property->entropy_64bits;

		totalBytes += property->totalByteSize;

		printf("  numOfElem = %zu\n", (size_t)property->numOfElem);
		printf("  entropy_8bits = %.1f\n", property->entropy_8bits);
		if (dataType == QCAT_FLOAT)
			printf("  entropy_32bits = %.1f\n", property->entropy_32bits);
		else
			printf("  entropy_64bits = %.1f\n", property->entropy_64bits);

		free(property);
		free(data);
	}
	fclose(propFile);
	printf("\nProperty output written to data-property.txt\n");

	double e8_avg = 0, e8_min = DBL_MAX, e8_max = -DBL_MAX;
	double efp_avg = 0, efp_min = DBL_MAX, efp_max = -DBL_MAX;

	for (int i = 0; i < nFiles; i++) {
		e8_avg += entropy_8bit[i];
		if (entropy_8bit[i] < e8_min) e8_min = entropy_8bit[i];
		if (entropy_8bit[i] > e8_max) e8_max = entropy_8bit[i];

		efp_avg += entropy_fp[i];
		if (entropy_fp[i] < efp_min) efp_min = entropy_fp[i];
		if (entropy_fp[i] > efp_max) efp_max = entropy_fp[i];
	}
	e8_avg /= nFiles;
	efp_avg /= nFiles;

	FILE *ymlFile = fopen("dataset_entry.yml", "w");
	if (!ymlFile) {
		fprintf(stderr, "Error: cannot create dataset_entry.yml\n");
		exit(1);
	}

	fprintf(ymlFile, "- name: \"%s\"\n", cfg.application_name);

	if (cfg.source[0])
		fprintf(ymlFile, "  source: \"%s\"\n", cfg.source);
	else
		fprintf(ymlFile, "  source: \"\"\n");

	fprintf(ymlFile, "  source_url: \"%s\"\n", cfg.source_url);

	if (cfg.source_note[0])
		fprintf(ymlFile, "  source_note: \"%s\"\n", cfg.source_note);

	if (cfg.type[0])
		fprintf(ymlFile, "  type: \"%s\"\n", cfg.type);
	else
		fprintf(ymlFile, "  type: \"\"\n");

	char dimStr[128] = "";
	for (int i = 0; i < cfg.ndims; i++) {
		char tmp[32];
		if (i > 0) strcat(dimStr, "x");
		snprintf(tmp, sizeof(tmp), "%zu", cfg.dims[i]);
		strcat(dimStr, tmp);
	}

	const char *precStr = (dataType == QCAT_FLOAT) ? "Single" : "Double";
	fprintf(ymlFile, "  format: |\n");
	fprintf(ymlFile, "    %d fields, %dD, %s\n", nFiles, cfg.ndims, dimStr);
	fprintf(ymlFile, "    %s precision, binary\n", precStr);

	const char *efp_key = (dataType == QCAT_FLOAT) ? "e32" : "e64";
	fprintf(ymlFile, "  entropy:\n");
	fprintf(ymlFile, "    - label: \"\"\n");
	fprintf(ymlFile, "      rows:\n");
	fprintf(ymlFile, "        - { stat: \"avg\", e8: \"%.1f\", %s: \"%.1f\" }\n",
	        e8_avg, efp_key, efp_avg);
	fprintf(ymlFile, "        - { stat: \"min\", e8: \"%.1f\", %s: \"%.1f\" }\n",
	        e8_min, efp_key, efp_min);
	fprintf(ymlFile, "        - { stat: \"max\", e8: \"%.1f\", %s: \"%.1f\" }\n",
	        e8_max, efp_key, efp_max);

	double sizeMB = totalBytes / (1024.0 * 1024.0);
	double sizeGB = totalBytes / (1024.0 * 1024.0 * 1024.0);
	if (sizeGB >= 1.0)
		fprintf(ymlFile, "  size: \"%.2fGB\"\n", sizeGB);
	else
		fprintf(ymlFile, "  size: \"%.0fMB\"\n", sizeMB);

	const char *typeFlag = (dataType == QCAT_FLOAT) ? "-f" : "-d";
	const char *typeStr = (dataType == QCAT_FLOAT) ? "float" : "double";

	char dimArgsSz[256] = "";
	for (int i = cfg.ndims - 1; i >= 0; i--) {
		char tmp[32];
		snprintf(tmp, sizeof(tmp), " %zu", cfg.dims[i]);
		strcat(dimArgsSz, tmp);
	}

	char dimArgsLp[256] = "";
	for (int i = cfg.ndims - 1; i >= 0; i--) {
		char tmp[64];
		snprintf(tmp, sizeof(tmp), " -d %zu", cfg.dims[i]);
		strcat(dimArgsLp, tmp);
	}

	fprintf(ymlFile, "  commands:\n");
	fprintf(ymlFile, "    sz_compress: \"sz -z %s -i %s -M REL -R 1E-2 -%d%s\"\n",
	        typeFlag, fileNames[0], cfg.ndims, dimArgsSz);
	fprintf(ymlFile, "    sz_decompress: \"sz -x %s -i %s -%d%s -s %s.sz -a\"\n",
	        typeFlag, fileNames[0], cfg.ndims, dimArgsSz, fileNames[0]);
	fprintf(ymlFile, "    zfp: \"zfp %s -i %s -z %s.zfp -o %s.zfp.out -%d%s -a 1E-2 -s\"\n",
	        typeFlag, fileNames[0], fileNames[0], fileNames[0], cfg.ndims, dimArgsSz);
	fprintf(ymlFile, "    libpressio: \"pressio -b compressor=$COMP -i %s%s -t %s"
	        " -o rel=1e-2 -m time -m size -M all\"\n",
	        fileNames[0], dimArgsLp, typeStr);
	fprintf(ymlFile, "    libpressio_note: \"where $COMP can be sz, zfp, sz3, mgard, etc.\"\n");
	fprintf(ymlFile, "    zchecker: \"./runZCCase.sh %s REL %s raw-data-dir %s%s\"\n",
	        typeFlag, cfg.application_name, cfg.file_extension, dimArgsSz);

	const char *baseUrl = "https://g-d0cd3f.fd635.8443.data.globus.org/raw-data";
	fprintf(ymlFile, "  links:\n");
	fprintf(ymlFile, "    - { label: \"Dataset\", url: \"%s/%s/SDRBENCH-%s-%s.tar.gz\" }\n",
	        baseUrl, cfg.application_name, cfg.application_name, dimStr);
	fprintf(ymlFile, "    - { label: \"Metadata\", url: \"%s/%s/template_data.txt\" }\n",
	        baseUrl, cfg.application_name);
	fprintf(ymlFile, "    - { label: \"Property\", url: \"%s/%s/SDRBENCH-%s-%s-property.txt\" }\n",
	        baseUrl, cfg.application_name, cfg.application_name, dimStr);

	fclose(ymlFile);

	printf("\nEntropy Summary:\n");
	printf("  8-bit Entropy:  avg=%.1f  min=%.1f  max=%.1f\n", e8_avg, e8_min, e8_max);
	printf("  %d-bit Entropy: avg=%.1f  min=%.1f  max=%.1f\n",
	       (dataType == QCAT_FLOAT) ? 32 : 64, efp_avg, efp_min, efp_max);
	printf("\nYAML entry written to dataset_entry.yml\n");

	free(entropy_8bit);
	free(entropy_fp);
	for (int i = 0; i < nFiles; i++)
		free(fileNames[i]);

	return 0;
}
