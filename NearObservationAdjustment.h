#pragma once
#include "mfhdf.h"
#include "math.h"
#include"omp.h"
#define PI  3.141592653

class NearObservationAdjustment
{
public:
	NearObservationAdjustment();
	~NearObservationAdjustment();

	int Main();

	/*    Get the Year and Day-of-year from the file names.
	*/
	void GetDataSeriesInfoFromList_HLS(CStringList &sDataFileList, int *piDataYear, int *piDataDay);

	/*    Linear regression.
	*/
	double regression(double *x, double *y, int count, double *slope, double *intercept);

	/*    Read the path of all files in the TXT file.
	*/
	void ReadAllFilePathsInTXT(CString sInputTXT, CStringList &sOutputFileList);

	/*    TRA for a single pixel.
	*
	*  Function input:
	*			inputFileName_L30:    the path of L30 point file (CSV).
	*			inputFileName_S30:    the path of S30 point file (CSV).
	*			lineLength_L30:           the number of lines in the L30 point file (exclude head line).
	*			lineLength_S30:            the number of lines in the S30 point file (exclude head line).
	*
	*  Function output:
	*			outputParaFile:             the path of output file (CSV).
	*/
	void TRA_MGRStile(CStringList &sListFile_L30, CStringList &sListFile_S30, int iTotalParts, CString outputFileName);

	/*    TRA for a single pixel.
	*
	*  Function input: The paths of L30 and S30 point file (CSV).
	*  Function output: The fitting parameters of the six bands (TXT).
	*/
	void TRA_PixelCSV(CString inputFileName_L30, CString inputFileName_S30,
		int lineLength_L30, int lineLength_S30, CString outputParaFile);
};

