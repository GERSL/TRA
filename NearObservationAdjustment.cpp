#include "stdafx.h"
#include "NearObservationAdjustment.h"


NearObservationAdjustment::NearObservationAdjustment(){}

NearObservationAdjustment::~NearObservationAdjustment(){}


int NearObservationAdjustment::Main()
{

	/*    TRA for a single pixel. 
	*/
	//CString inputFileName_L30 = "exampleData_L30.csv";
	//CString inputFileName_S30 = "exampleData_S30.csv"; 
	//int lineLength_L30 = 196;
	//int lineLength_S30 = 216;
	//CString outputParaFile = "exampleOutput.csv";
	//TRA_PixelCSV(inputFileName_L30, inputFileName_S30, lineLength_L30, lineLength_S30, outputParaFile);



	/*    TRA for a MGRS tile.
	*/
	CString path_allL30File = "H:\\SR_Data\\HLS\\V1.4\\TIF\\AllL30.txt";
	CString path_allS30File = "H:\\SR_Data\\HLS\\V1.4\\TIF\\AllS30A.txt";
	CString outputFileName = "exampleTile";
	int iTotalParts = 3;   // Seperate the tile into three parts
	CStringList sfinalList_L30, sfinalList_S30;
	ReadAllFilePathsInTXT(path_allL30File, sfinalList_L30);
	ReadAllFilePathsInTXT(path_allS30File, sfinalList_S30);
	TRA_MGRStile(sfinalList_L30, sfinalList_S30, iTotalParts, outputFileName);
	
	return 0;
}



/*    Get the Year and Day-of-year from the file names.
*/
void NearObservationAdjustment::GetDataSeriesInfoFromList_HLS(CStringList &sDataFileList, int *piDataYear, int *piDataDay)
{
	int iCount = sDataFileList.GetCount();
	for (int iF = 0; iF<iCount; iF++)
	{
		CString sFullURLFile = sDataFileList.GetAt(sDataFileList.FindIndex(iF));
		CString sFileName, sDataDate;
		sFileName = sFullURLFile.Right(sFullURLFile.GetLength() - sFullURLFile.ReverseFind('\\') - 1); //HLS.L30.T17RKQ.2015003.v1.3.hdf
		sDataDate = sFileName;
		sDataDate = sDataDate.Right(sDataDate.GetLength() - sDataDate.Find(".") - 1); //L30.T17RKQ.2015003.v1.3.hdf
		sDataDate = sDataDate.Right(sDataDate.GetLength() - sDataDate.Find(".") - 1); //T17RKQ.2015003.v1.3.hdf
		sDataDate = sDataDate.Right(sDataDate.GetLength() - sDataDate.Find(".") - 1); //2015003.v1.3.hdf
		int iDataYear = atoi(sDataDate.Mid(0, 4));
		int iDataJulianDay = atoi(sDataDate.Mid(4, 3));
		sDataDate = sDataDate.Right(sDataDate.GetLength() - sDataDate.Find(".") - 1);

		int iTileH = atoi(sDataDate.Mid(1, 2));
		int iTileV = atoi(sDataDate.Mid(4, 2));

		piDataYear[iF] = iDataYear;
		piDataDay[iF] = iDataJulianDay;
	}//***iF
}

/*    Read the path of all files in the TXT file.
*/
void NearObservationAdjustment::ReadAllFilePathsInTXT(CString sInputTXT, CStringList &sOutputFileList)
{
	CString sFileName, sDataDate;
	char sCompositeFile[1024];
	//CStringList sFilesList_L30;
	FILE *fp_L30;
	fp_L30 = fopen(sInputTXT, "r");
	if (!fp_L30) return;
	for (;;)
	{
		if (feof(fp_L30))break;
		if (fscanf(fp_L30, "%s", sCompositeFile) < 1)continue;
		sOutputFileList.AddTail(sCompositeFile);
	}
	fclose(fp_L30);
	return;
}

/*    Linear regression.
*/
double NearObservationAdjustment::regression(double *x, double *y, int count, double *slope, double *intercept)
{
	double meanX, meanY, SSXX, SSYY, SSXY;
	int i, j;
	double a, b, r_square;
	meanX = 0;
	meanY = 0;
	SSXX = 0;
	SSYY = 0;
	SSXY = 0;
	for (i = 0; i < count; i++)
	{
		meanX = meanX + x[i];
		meanY = meanY + y[i];
	}
	meanX = meanX / count;
	meanY = meanY / count;

	for (i = 0; i < count; i++)
	{
		SSXX = SSXX + pow(x[i] - meanX, 2);
		SSYY = SSYY + pow(y[i] - meanY, 2);
		SSXY = SSXY + (x[i] - meanX)*(y[i] - meanY);
	}
	*slope = 9999;
	*intercept = 9999;

	if (SSXX == 0 || SSYY == 0)
		r_square = 0;
	else
		r_square = pow(SSXY, 2) / SSXX / SSYY;

	if (SSXX != 0)
	{
		b = SSXY / SSXX;
		a = meanY - b * meanX;
		if (a != 0)
		{
			*slope = b;
			*intercept = a;
		}
		else
		{
			*slope = b;
			*intercept = 0;
		}
	}
	return r_square;
}



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
void NearObservationAdjustment::TRA_PixelCSV(CString inputFileName_L30, CString inputFileName_S30,
	 int lineLength_L30,int lineLength_S30, CString outputParaFile)
{
	// The transformation parameters in the HLS V1.4 product (See HLS.v1.4.UserGuide.pdf).
	double paraBandAdjust_slope_V14A[6] = { 0.9778, 1.0053, 0.9765, 0.9983, 0.9987, 1.003 };
	double paraBandAdjust_intercept_V14A[6] = { -0.004, -0.0009, 0.0009, -0.0001, -0.0011, -0.0012 };
	CString fileHead;

	//---------------------------------------------------------
	//-------------- Read S30 of a point -------------
	//---------------------------------------------------------
	int *piS30Year = (int *)malloc(lineLength_S30 * sizeof(int));
	memset(piS30Year, 0, lineLength_S30 * sizeof(int));
	int *piS30Doy = (int *)malloc(lineLength_S30 * sizeof(int));
	memset(piS30Doy, 0, lineLength_S30 * sizeof(int));
	int *piS30QA = (int *)malloc(lineLength_S30 * sizeof(int));
	memset(piS30QA, 0, lineLength_S30 * sizeof(int));
	int **piS30SR = new int*[lineLength_S30];
	for (int iL = 0; iL < lineLength_S30; iL++)
	{
		piS30SR[iL] = new int[6];
		memset(piS30SR[iL], 0, sizeof(int) * 6);
	}
	FILE * fpw_s = fopen(inputFileName_S30, "r");
	fscanf(fpw_s, "%s\n", &fileHead); // Skip head
	for (int iL = 0; iL < lineLength_S30; iL++)
	{
		fscanf(fpw_s, "%d,%d,%d,%d,%d,%d,%d,%d,%d\n",
			&piS30Year[iL], &piS30Doy[iL], &piS30SR[iL][0],
			&piS30SR[iL][1], &piS30SR[iL][2], &piS30SR[iL][3],
			&piS30SR[iL][4], &piS30SR[iL][5], &piS30QA[iL]);

		//printf("%d\n", piS30SR[iL][4]);
	}
	fclose(fpw_s);

	//---------------------------------------------------------
	//-------------- Read L30 of a point -------------
	//---------------------------------------------------------
	int *piL30Year = (int *)malloc(lineLength_L30 * sizeof(int));
	memset(piL30Year, 0, lineLength_L30 * sizeof(int));
	int *piL30Doy = (int *)malloc(lineLength_L30 * sizeof(int));
	memset(piL30Doy, 0, lineLength_L30 * sizeof(int));
	int *piL30QA = (int *)malloc(lineLength_L30 * sizeof(int));
	memset(piL30QA, 0, lineLength_L30 * sizeof(int));
	int **piL30SR = new int*[lineLength_L30];
	for (int iL = 0; iL < lineLength_L30; iL++)
	{
		piL30SR[iL] = new int[6];
		memset(piL30SR[iL], 0, sizeof(int) * 6);
	}
	FILE *fpw_p= fopen(inputFileName_L30, "r");
	fscanf(fpw_p, "%s\n", &fileHead); // Skip head
	for (int iL = 0; iL < lineLength_L30; iL++)
	{
		fscanf(fpw_p, "%d,%d,%d,%d,%d,%d,%d,%d,%d\n",
			&piL30Year[iL], &piL30Doy[iL], &piL30SR[iL][0],
			&piL30SR[iL][1], &piL30SR[iL][2], &piL30SR[iL][3],
			&piL30SR[iL][4], &piL30SR[iL][5], &piL30QA[iL]);

		//printf("%d\n", piL30SR[iL][4]);
	}
	fclose(fpw_p);




	//	Create output file
	FILE * outFile = fopen(outputParaFile, "w");
	fprintf(outFile, "BandIndex,Slope,Intercept,R2\n");

	//---------------------------------------------------------
	//-------------- Process for each band -------------
	//---------------------------------------------------------
	for (int iB = 0; iB < 6; iB++)
	{
		int iTotalCount = 365 * 4 + 121; // from 20150101 to 20190430
		double *piL30 = (double *)malloc(iTotalCount * sizeof(double));
		memset(piL30, 0, iTotalCount * sizeof(double));
		for (int iL = 0; iL < lineLength_L30; iL++)
		{
			int tmpIndex_L30 = (piL30Year[iL] - 2015) * 365 + piL30Doy[iL] - 1 + int((piL30Year[iL] - 2015 + 3) / 4);
			if (piL30QA[iL] == 1)		piL30[tmpIndex_L30] = piL30SR[iL][iB] * 0.0001;
		}

		int16 OrigRMSD = 0, LinearRMSD = 0;
		double OrigRMSDsum = 0.0, LinearRMSDsum = 0.0;
		int OrigRMSDcount = 0, LinearRMSDcount = 0, outFlag = 0;

		//----------------------------------------------------------------------
		//-------------- One-day match of L30 and S30 --------------
		//----------------------------------------------------------------------
		double *piS30 = (double *)malloc(lineLength_S30 * sizeof(double));
		memset(piS30, 0, lineLength_S30 * sizeof(double));
		double *piL30inter = (double *)malloc(lineLength_S30 * sizeof(double));
		memset(piL30inter, 0, lineLength_S30 * sizeof(double));
		for (int iL = 0; iL < lineLength_S30; iL++)
		{
			if (piS30Year[iL] == 0) continue;

			if (piS30QA[iL] == 1)		piS30[iL] = piS30SR[iL][iB] * 0.0001;

			int tmpIndex_S30 = (piS30Year[iL] - 2015) * 365 + piS30Doy[iL] - 1 + int((piS30Year[iL] - 2015 + 3) / 4);

			// Find valid L30 observation within 1 days of S30
			if (piL30[tmpIndex_S30] > 0) { piL30inter[iL] = piL30[tmpIndex_S30]; continue; }
			if (piL30[tmpIndex_S30 - 1] > 0) { piL30inter[iL] = piL30[tmpIndex_S30 - 1]; continue; }
			if (piL30[tmpIndex_S30 + 1] > 0) { piL30inter[iL] = piL30[tmpIndex_S30 + 1]; continue; }
		}

		//----------------------------------------------------------------------
		//-------- Pixel-level fitting between S30 and L30 -------
		//-----------------------------------------------------------------------
		double *piFittingX = new double[lineLength_S30];
		memset(piFittingX, 0, lineLength_S30 * sizeof(double));
		double *piFittingY = new double[lineLength_S30];		
		memset(piFittingY, 0, lineLength_S30 * sizeof(double));
		int iFittingCount = 0;
		for (int iL = 0; iL < lineLength_S30; iL++)
		{			
			if (piS30Year[iL] == 0) continue;		//Non-date screening			
			if (piS30[iL] <= 0.0 || piL30inter[iL] <= 0.0) continue;	//Invalid screening

			// Original sentinel----> fitting X
			double tempX = (piS30[iL] - paraBandAdjust_intercept_V14A[iB]) / paraBandAdjust_slope_V14A[iB];
			double tempY = piL30inter[iL];
			if (tempX <= 0) continue; //Invalid original screening

			// Use the blue band reflectance to screen remaining clouds (Zhang et al., 2018)
			if (iB == 0) 
			{				
				if (fabs(tempX - tempY) > 0.5*(tempX + tempY)) piS30QA[iL] = 2;
			}
			if (piS30QA[iL] == 2) continue; //QA screening

			//  Save valid fitting X and Y
			piFittingX[iFittingCount] = tempX;
			piFittingY[iFittingCount] = tempY;  
			iFittingCount++;  //fitting count
		}//iL

		 /***** calculate the regression parameter *****/
		double R_square, slope = 9999, intercept = 9999;
		if (iFittingCount >= 4)	{
			R_square = regression(piFittingX, piFittingY, iFittingCount, &slope, &intercept);
		}			
		else 	{
			//----------------------------------------------------------------------------
			//-------------- Interpolating match of L30 and S30 --------------
			//----------------------------------------------------------------------------
			memset(piS30, 0, lineLength_S30 * sizeof(double));
			memset(piL30inter, 0, lineLength_S30 * sizeof(double));
			for (int iL = 0; iL < lineLength_S30; iL++)
			{
				if (piS30Year[iL] == 0) continue;

				if (piS30QA[iL] == 1)		piS30[iL] = piS30SR[iL][iB] * 0.0001;

				int tmpIndex_S30 = (piS30Year[iL] - 2015) * 365 + piS30Doy[iL] - 1 + int((piS30Year[iL] - 2015 + 3) / 4);

				// Find valid L30 observation within 1 days of S30
				if (piL30[tmpIndex_S30] > 0) { piL30inter[iL] = piL30[tmpIndex_S30]; continue; }
				if (piL30[tmpIndex_S30 - 1] > 0) { piL30inter[iL] = piL30[tmpIndex_S30 - 1]; continue; }
				if (piL30[tmpIndex_S30 + 1] > 0) { piL30inter[iL] = piL30[tmpIndex_S30 + 1]; continue; }

				// Find valid L30 observation within 16 days of S30 in backward direction
				int iStep = 2;
				double interLeft = 0.0;
				int interLeft_pos = 0;
				while (iStep < 16)
				{
					if (piL30[tmpIndex_S30 - iStep] > 0)
					{
						interLeft = piL30[tmpIndex_S30 - iStep];
						interLeft_pos = tmpIndex_S30 - iStep;
						break;
					}
					iStep = iStep + 1;
				}

				// Find valid L30 observation within 16 days of S30 in forward direction
				iStep = 2;
				double interRight = 0.0;
				int interRight_pos = 0;
				while (iStep < 16)
				{
					if (piL30[tmpIndex_S30 + iStep] > 0)
					{
						interRight = piL30[tmpIndex_S30 + iStep];
						interRight_pos = tmpIndex_S30 + iStep;
						break;
					}
					iStep = iStep + 1;
				}

				// Linear interpolation
				if (interLeft != 0.0 && interRight != 0.0)
				{
					piL30inter[iL] = interLeft + (interRight - interLeft) * (tmpIndex_S30 - interLeft_pos) / (interRight_pos - interLeft_pos);
				}
			}

			//----------------------------------------------------------------------
			//-------- Pixel-level fitting between S30 and L30 -------
			//-----------------------------------------------------------------------
			memset(piFittingX, 0, lineLength_S30 * sizeof(double));
			memset(piFittingY, 0, lineLength_S30 * sizeof(double));
			iFittingCount = 0;
			for (int iL = 0; iL < lineLength_S30; iL++)
			{
				if (piS30Year[iL] == 0) continue;		//Non-date screening			
				if (piS30[iL] <= 0.0 || piL30inter[iL] <= 0.0) continue;	//Invalid screening

				// Original sentinel----> fitting X
				double tempX = (piS30[iL] - paraBandAdjust_intercept_V14A[iB]) / paraBandAdjust_slope_V14A[iB];
				double tempY = piL30inter[iL];
				if (tempX <= 0) continue; //Invalid original screening

				// Use the blue band reflectance to screen remaining clouds (Zhang et al., 2018)
				if (iB == 0)
				{
					if (fabs(tempX - tempY) > 0.5*(tempX + tempY)) piS30QA[iL] = 2;
				}
				if (piS30QA[iL] == 2) continue; //QA screening

				//  Save valid fitting X and Y
				piFittingX[iFittingCount] = tempX;
				piFittingY[iFittingCount] = tempY;
				iFittingCount++;  //fitting count
			}//iL
			R_square = regression(piFittingX, piFittingY, iFittingCount, &slope, &intercept);
		}

		// Save the output fitting parameters
		fprintf(outFile, "%d,%7.4f,%7.4f,%7.4f\n", iB+1, slope, intercept, R_square);
	
		delete piFittingX;
		delete piFittingY;
		delete piS30;
		delete piL30inter;
		delete piL30;
	}//iB


	delete piS30Year;
	delete piS30Doy;
	delete piS30QA;
	for (int iL = 0; iL < lineLength_S30; iL++) { delete piS30SR[iL]; }
	delete piS30SR;

	delete piL30Year;
	delete piL30Doy;
	delete piL30QA;
	for (int iL = 0; iL < lineLength_L30; iL++) { delete piL30SR[iL]; }
	delete piL30SR;

	fclose(outFile);
	return;
}



/*    TRA for a MGRS tile.
*
*  Function input:
*			sListFile_L30:    the path of all L30  HDF file (TXT).
*			sListFile_S30:    the path of all S30  HDF file (TXT).
*			iTotalParts:        the number of parts deviding a MDRS tile
*
*  Function output:
*			outputFileName:  the name of output file .
*/
void NearObservationAdjustment::TRA_MGRStile(CStringList &sListFile_L30,
	CStringList &sListFile_S30, int iTotalParts,CString outputFileName)
{
	//---------------------- function input -----------------------
	//---- sListFile_L30:   L30 file path list
	//---- sListFile_S30:   S30 file path list
	//------------------------------------------------------------------

	//------------------------------------     blue          green        red          NIR        SWIR-1    SWIR-2 
	CString Landsat8_Bands[] = { "band02" ,"band03","band04","band05" ,"band06","band07" };
	CString Sentinel2_Bands[] = { "B02" ,       "B03",       "B04",      "B8A" ,       "B11",      "B12" };
	CString Landsat8_outputBands[] = { "blue" ,"green","red","NIR" ,"SWIR1","SWIR2" };
	double paraBandAdjust_slope_V14A[6] = { 0.9778, 1.0053, 0.9765, 0.9983, 0.9987, 1.003 };
	double paraBandAdjust_intercept_V14A[6] = { -0.004, -0.0009, 0.0009, -0.0001, -0.0011, -0.0012 };
	int iTotalCount = 365 * 4 + 121; // from 20150101 to 20190430
	int32 start[2] = { 0, 0 };
	int32 end[2] = { 3660, 3660 };  // Total size




	//Get the year and DOY of each L30 file
	int iFilesCount_L30 = sListFile_L30.GetCount();
	int *piDataYear = (int *)malloc(iFilesCount_L30 * sizeof(int));
	int *piDataDate = (int *)malloc(iFilesCount_L30 * sizeof(int));
	memset(piDataYear, 0, iFilesCount_L30 * sizeof(int));
	memset(piDataDate, 0, iFilesCount_L30 * sizeof(int));
	GetDataSeriesInfoFromList_HLS(sListFile_L30, piDataYear, piDataDate);


	//Get the year and DOY of each S30 file
	int iFilesCount_S30 = sListFile_S30.GetCount();
	int *piDataYear_S30 = (int *)malloc(iFilesCount_S30 * sizeof(int));
	int *piDataDate_S30 = (int *)malloc(iFilesCount_S30 * sizeof(int));
	memset(piDataYear_S30, 0, iFilesCount_S30 * sizeof(int));
	memset(piDataDate_S30, 0, iFilesCount_S30 * sizeof(int));
	GetDataSeriesInfoFromList_HLS(sListFile_S30, piDataYear_S30, piDataDate_S30);



	//------------------------------------------------------------------
	// Read one int16 layer (3660*3660) needs 25.56 Mb.
	// For 400 int16 + uint8 layers, it needs ~ 15 Gb.
	// Therefore, a MGRS tile should be separated.
	//------------------------------------------------------------------
	for (int iPart = 1; iPart <= iTotalParts; iPart++)
	{
		printf("Processing part %d ...\n", iPart);

		if (iPart == 1)
		{
			//--------------------------------------------------------------------
			//-------------- Create empty output HDF file -------------
			//--------------------------------------------------------------------
			for (int iB = 0; iB < 6; iB++)
			{
				CString sOutputFileName;
				sOutputFileName.Format("%s.FittingPara%s.hdf", outputFileName, Landsat8_outputBands[iB]);
				printf("Create file: %s\n", sOutputFileName);

				double *pOutDouble = new double[end[0] * end[1]];
				memset(pOutDouble, 0, sizeof(double)*end[0] * end[1]);
				int32 *pOutInt16 = new int32[end[0] * end[1]];
				memset(pOutInt16, 0, sizeof(int32)*end[0] * end[1]);
				uint8 *pUInt8 = new uint8[end[0] * end[1]];
				memset(pUInt8, 0, sizeof(uint8)*end[0] * end[1]);

				int32 sd_idw = SDstart(sOutputFileName, DFACC_CREATE);
				int32 sds_idw = SDcreate(sd_idw, "Slope", DFNT_DOUBLE, 2, end);
				SDwritedata(sds_idw, start, NULL, end, pOutDouble);
				SDendaccess(sds_idw);
				sds_idw = SDcreate(sd_idw, "Intercept", DFNT_DOUBLE, 2, end);
				SDwritedata(sds_idw, start, NULL, end, pOutDouble);
				SDendaccess(sds_idw);
				sds_idw = SDcreate(sd_idw, "FittingCount", DFNT_INT16, 2, end);
				SDwritedata(sds_idw, start, NULL, end, pOutInt16);
				SDendaccess(sds_idw);
				sds_idw = SDcreate(sd_idw, "QA", DFNT_UINT8, 2, end);
				SDwritedata(sds_idw, start, NULL, end, pUInt8);
				SDendaccess(sds_idw);
				SDend(sd_idw);
				delete pOutInt16;
				delete pOutDouble;
				delete pUInt8;
			}
		}

		int32 sizeY = 3660 / iTotalParts;
		int32 start_p[2] = { 0, (iPart - 1)*sizeY };  // start 
		int32 iSrcDim[2] = { 3660, sizeY }; //Size for reading

		
		//---------------------------------------------------------
		//-------------- Read Sentinel-2 QA -------------
		//---------------------------------------------------------
		printf("Reading Sentinel QA (%d files):\n", iFilesCount_S30);
		uint8 **piS30_QA = new uint8*[iFilesCount_S30];
		int32 sd_id, sds_id;
		CString inputFileName_S30;
       #pragma omp parallel for num_threads(16) 
		for (int iL = 0; iL < iFilesCount_S30; iL++)
		{
			piS30_QA[iL] = new uint8[iSrcDim[0] * iSrcDim[1]];
			memset(piS30_QA[iL], 0, sizeof(uint8)*iSrcDim[0] * iSrcDim[1]);

			//Read Sentinel QA
			uint8 *tempQA = new uint8[iSrcDim[0] * iSrcDim[1]];
			memset(tempQA, 0, sizeof(uint8)*iSrcDim[0] * iSrcDim[1]);
			inputFileName_S30 = sListFile_S30.GetAt(sListFile_S30.FindIndex(iL));
			sd_id = SDstart(inputFileName_S30, DFACC_READ);
			if (sd_id < 0) { SDend(sd_id);	delete tempQA; continue; }   // data does not exist, continue	
			printf(" S%03d", iL + 1); 
			sds_id = SDselect(sd_id, SDnametoindex(sd_id, "QA"));
			SDreaddata(sds_id, start_p, NULL, iSrcDim, tempQA);
			SDendaccess(sds_id);
			SDend(sd_id);

			//Simplify  Sentinel-2 QA
			for (int iP = 0; iP < iSrcDim[0] * iSrcDim[1]; iP++)
			{
				uint8 qaCirrusFlag_s, qaCloudFlag_s, qaCShadowFlag_s, qaSnowFlag_s, qaWaterFlag_s, finalFlag_s = 0;
				qaCirrusFlag_s = tempQA[iP] % 2;
				qaCloudFlag_s = (tempQA[iP] >> 1) % 2;
				qaCShadowFlag_s = (tempQA[iP] >> 3) % 2;
				qaSnowFlag_s = (tempQA[iP] >> 4) % 2;
				qaWaterFlag_s = (tempQA[iP] >> 5) % 2;
				if (qaCirrusFlag_s == 0 && qaCloudFlag_s == 0 && qaCShadowFlag_s == 0 && qaSnowFlag_s == 0 && qaWaterFlag_s == 0) piS30_QA[iL][iP] = 1;
			}//iP
			delete tempQA;
		}//iL

		//---------------------------------------------------------
		//-------------- Read Landsat-8 QA -------------
		//---------------------------------------------------------		
		printf("\nReading Landsat QA (%d files):\n", iFilesCount_L30);
		uint8 **piL30_QA = new uint8*[iFilesCount_L30];
		CString  inputFileName_L30;
		int32 sd_id_L30, sds_id_L30;
		#pragma omp parallel for num_threads(16) 
		for (int iL = 0; iL < iFilesCount_L30; iL++)
		{
			piL30_QA[iL] = new uint8[iSrcDim[0] * iSrcDim[1]];
			memset(piL30_QA[iL], 0, sizeof(uint8)*iSrcDim[0] * iSrcDim[1]);

			//Read Landsat QA
			uint8 *tempQA = new uint8[iSrcDim[0] * iSrcDim[1]];
			memset(tempQA, 0, sizeof(uint8)*iSrcDim[0] * iSrcDim[1]);
			inputFileName_L30 = sListFile_L30.GetAt(sListFile_L30.FindIndex(iL));
			sd_id_L30 = SDstart(inputFileName_L30, DFACC_READ);
			if (sd_id_L30 < 0) { SDend(sd_id_L30);	delete tempQA; continue; }   // data does not exist, continue	
			printf(" L%03d", iL + 1); 
			sds_id_L30 = SDselect(sd_id_L30, SDnametoindex(sd_id_L30, "QA"));
			SDreaddata(sds_id_L30, start_p, NULL, iSrcDim, tempQA);
			SDendaccess(sds_id_L30);
			SDend(sd_id_L30);

			//Simplify  Landsat QA
			for (int iP = 0; iP < iSrcDim[0] * iSrcDim[1]; iP++)
			{
				uint8 qaCirrusFlag, qaCloudFlag, qaACloudFlag, qaCShadowFlag, qaSnowFlag, qaWaterFlag, qaAerosolFlag, finalFlag = 0;
				qaCirrusFlag = tempQA[iP] % 2;
				qaCloudFlag = (tempQA[iP] >> 1) % 2;
				qaACloudFlag = (tempQA[iP] >> 2) % 2;
				qaCShadowFlag = (tempQA[iP] >> 3) % 2;
				qaSnowFlag = (tempQA[iP] >> 4) % 2;
				qaWaterFlag = (tempQA[iP] >> 5) % 2;
				qaAerosolFlag = (tempQA[iP] >> 6) % 4; 
				if (qaCirrusFlag == 0 && qaCloudFlag == 0 && qaACloudFlag == 0 && qaCShadowFlag == 0 &&
					qaSnowFlag == 0 && qaWaterFlag == 0 && qaAerosolFlag < 3) piL30_QA[iL][iP] = 1;
			}//iP
			delete tempQA;
		}//iF
		printf("\n");

		//---------------------------------------------------------
		//-------------- Process for each band -----------
		//---------------------------------------------------------
		for (int iB = 0; iB < 6; iB++)
		{
			printf("Processing %s band...\n", Landsat8_outputBands[iB]);
			
			// Define output array
			double *pOutSlope = new double[iSrcDim[0] * iSrcDim[1]];
			double *pOutIntercept = new double[iSrcDim[0] * iSrcDim[1]];
			for (int iP = 0; iP < iSrcDim[0] * iSrcDim[1]; iP++)
			{
				pOutSlope[iP] = 9999;  pOutIntercept[iP] = 9999;
			}
			int16 *pOutFittingCount = new int16[iSrcDim[0] * iSrcDim[1]];
			memset(pOutFittingCount, 0, sizeof(int16)*iSrcDim[0] * iSrcDim[1]);
			uint8 *pOutFittingQA = new uint8[iSrcDim[0] * iSrcDim[1]];
			memset(pOutFittingQA, 0, sizeof(uint8)*iSrcDim[0] * iSrcDim[1]);

		

			//---------------------------------------------------------
			//-------- Read Sentinel-2 reflectance ---------
			//---------------------------------------------------------
			printf("Read Sentinel-2 reflectance (%d files)...\n", iFilesCount_S30);
			int16 **piS30_Band = new int16*[iFilesCount_S30];
			#pragma omp parallel for num_threads(16) 
			for (int iF = 0; iF < iFilesCount_S30; iF++)
			{
				piS30_Band[iF] = new int16[iSrcDim[0] * iSrcDim[1]];
				memset(piS30_Band[iF], 0, sizeof(int16)*iSrcDim[0] * iSrcDim[1]);

				//Read Sentinel
				inputFileName_S30 = sListFile_S30.GetAt(sListFile_S30.FindIndex(iF));
				sd_id = SDstart(inputFileName_S30, DFACC_READ);
				if (sd_id < 0) { SDend(sd_id);	continue; }
				printf(" S%03d", iF + 1); //piDataYear_S30[iF], piDataDate_S30[iF]);
				sds_id = SDselect(sd_id, SDnametoindex(sd_id, Sentinel2_Bands[iB]));
				SDreaddata(sds_id, start_p, NULL, iSrcDim, piS30_Band[iF]);
				SDendaccess(sds_id);
				SDend(sd_id);
			}//iF
			printf("\n");

			//---------------------------------------------------------
			//-------- Read Landsat-8 reflectance ---------
			//---------------------------------------------------------
			printf("Read Landsat-8 reflectance (%d files)...\n", iFilesCount_L30);
			int16 **piL30_Band = new int16*[iFilesCount_L30];
			#pragma omp parallel for num_threads(16) 
			for (int iF = 0; iF < iFilesCount_L30; iF++)
			{
				piL30_Band[iF] = new int16[iSrcDim[0] * iSrcDim[1]];
				memset(piL30_Band[iF], 0, sizeof(int16)*iSrcDim[0] * iSrcDim[1]);

				//Read Landsat
				inputFileName_L30 = sListFile_L30.GetAt(sListFile_L30.FindIndex(iF));
				sd_id_L30 = SDstart(inputFileName_L30, DFACC_READ);
				if (sd_id_L30 < 0) { SDend(sd_id_L30);	continue; }
				printf(" L%03d", iF + 1); //piDataYear[iF], piDataDate[iF]);
				sds_id_L30 = SDselect(sd_id_L30, SDnametoindex(sd_id_L30, Landsat8_Bands[iB]));
				SDreaddata(sds_id_L30, start_p, NULL, iSrcDim, piL30_Band[iF]);
				SDendaccess(sds_id_L30);
				SDend(sd_id_L30);
			}//iF
			printf("\nReading ends, pixel-level processing starts...\n");

			 //----------------------------------------------------------------------
			//-------- Pixel-level fitting between S30 and L30 -------
			//-----------------------------------------------------------------------
			for (int iP = 0; iP < iSrcDim[0] * iSrcDim[1]; iP++)
			{
				if (iP % (iSrcDim[0] * iSrcDim[1] / 100 - 1) == 0) printf("*%d*", iP / (iSrcDim[0] * iSrcDim[1] / 100 - 1));

				// Extract Landsat-8 time series of a pixel
				// ------ Use a daily array to save Landsat-8 time series
				double *piL30 = (double *)malloc(iTotalCount * sizeof(double));
				memset(piL30, 0, iTotalCount * sizeof(double));
				for (int iL = 0; iL < iFilesCount_L30; iL++)
				{
					//Invalid screening of Landsat-8
					if (piL30_QA[iL][iP] != 1 || piL30_Band[iL][iP] <= 0 || piL30_Band[iL][iP] > 10000) continue; 
					int tmpIndex_L30 = (piDataYear[iL] - 2015) * 365 + (piDataDate[iL] - 1)
						+ int((piDataYear[iL] - 2015 + 2) / 4);
					piL30[tmpIndex_L30] = piL30_Band[iL][iP] * 0.0001;
				}

				// Extract Sentinel-2 time series of a pixel
				double *piS30 = (double *)malloc(iFilesCount_S30 * sizeof(double));
				memset(piS30, 0, iFilesCount_S30 * sizeof(double));
				for (int iL = 0; iL < iFilesCount_S30; iL++)
				{
					//Invalid screening of Sentinel-2
					if (piS30_QA[iL][iP] != 1 || piS30_Band[iL][iP] <= 0 || piS30_Band[iL][iP] > 10000) continue;
					piS30[iL] = piS30_Band[iL][iP] * 0.0001;
				}

				//---------------------------------------------------------------
				//-------- Linear matching: first 1D matching -------
				//---------------------------------------------------------------
				// Array used for linear matching (size = iFilesCount_S30)
				double *piL30inter = new double[iFilesCount_S30];
				memset(piL30inter, 0, iFilesCount_S30 * sizeof(double));
				for (int iL = 0; iL < iFilesCount_S30; iL++)
				{
					int tmpIndex_S30 = (piDataYear_S30[iL] - 2015) * 365 + (piDataDate_S30[iL] - 1)
						+ int((piDataYear_S30[iL] - 2015 + 2) / 4);
					//printf("%d\n%d\n", iL, tmpIndex_S30);
										
					// Find valid L30 observation within 1 days of S30
					if (piL30[tmpIndex_S30] > 0) { piL30inter[iL] = piL30[tmpIndex_S30]; continue; }
					if (piL30[tmpIndex_S30 - 1] > 0) { piL30inter[iL] = piL30[tmpIndex_S30 - 1]; continue; }
					if (piL30[tmpIndex_S30 + 1] > 0) { piL30inter[iL] = piL30[tmpIndex_S30 + 1]; continue; }
				}

				//----------------------------------------------------------------------
				//-------- Extract the valid pairs of preS30 and L30 -------
				//-----------------------------------------------------------------------
				double *piX = new double[iFilesCount_S30];
				memset(piX, 0, iFilesCount_S30 * sizeof(double));
				double *piY = new double[iFilesCount_S30];
				memset(piY, 0, iFilesCount_S30 * sizeof(double));
				int iAllCount = 0;
				for (int iL = 0; iL < iFilesCount_S30; iL++)
				{
					if (piS30[iL] <= 0.0 || piL30inter[iL] <= 0.0) continue; //Invalid screening

					double tempX = (piS30[iL] - paraBandAdjust_intercept_V14A[iB]) / paraBandAdjust_slope_V14A[iB]; //Original 	
					double tempY = piL30inter[iL];
					if (tempX <= 0) continue; //Invalid original screening
					//Use the blue band reflectance to screen remaining clouds (Zhang et al., 2018)
					if (iB == 0) //Blue band
					{						
						if (fabs(tempX - tempY) > 0.5*(tempX + tempY)) piS30_QA[iL][iP] = 2;
					}
					if (fabs(tempX - tempY) > 0.5*(tempX + tempY)) continue; //Noise screening
					if (piS30_QA[iL][iP] == 2) continue; //QA screening

					piX[iAllCount] = tempX;
					piY[iAllCount] = tempY;  //fitting Y
					iAllCount++;  //fitting count
				}//iL

				if (iAllCount >= 4) pOutFittingQA[iP] = 1;  // indicate one-day matching
				if (iAllCount < 4)   // one-day matching < 4, try linear interpolation matching
				{					
					memset(piL30inter, 0, iFilesCount_S30 * sizeof(double));
					for (int iL = 0; iL < iFilesCount_S30; iL++)
					{
						int tmpIndex_S30 = (piDataYear_S30[iL] - 2015) * 365 + (piDataDate_S30[iL] - 1)
							+ int((piDataYear_S30[iL] - 2015 + 2) / 4);

						// Find valid L30 observation within 1 days of S30
						if (piL30[tmpIndex_S30] > 0) { piL30inter[iL] = piL30[tmpIndex_S30]; continue; }
						if (piL30[tmpIndex_S30 - 1] > 0) { piL30inter[iL] = piL30[tmpIndex_S30 - 1]; continue; }
						if (piL30[tmpIndex_S30 + 1] > 0) { piL30inter[iL] = piL30[tmpIndex_S30 + 1]; continue; }

						// Find valid L30 observation within 16 days of S30 in backward direction
						int iStep = 2;
						double interLeft = 0.0;
						int interLeft_pos = 0;
						while (iStep < 16)
						{
							if (piL30[tmpIndex_S30 - iStep] > 0)
							{
								interLeft = piL30[tmpIndex_S30 - iStep];
								interLeft_pos = tmpIndex_S30 - iStep;
								break;
							}
							iStep = iStep + 1;
						} // step left

						// Find valid L30 observation within 16 days of S30 in forward direction
						iStep = 2;
						double interRight = 0.0;
						int interRight_pos = 0;
						while (iStep < 16)
						{
							if (piL30[tmpIndex_S30 + iStep] > 0)
							{
								interRight = piL30[tmpIndex_S30 + iStep];
								interRight_pos = tmpIndex_S30 + iStep;
								break;
							}
							iStep = iStep + 1;
						} // step right

						// Linear interpolation
						if (interLeft != 0.0 && interRight != 0.0)
						{
							piL30inter[iL] = interLeft + (interRight - interLeft) * (tmpIndex_S30 - interLeft_pos) / (interRight_pos - interLeft_pos);
						} // linear interpolation						
					} // linear matching

					//-------- Extract the valid pairs of preS30 and L30
					memset(piX, 0, iFilesCount_S30 * sizeof(double));
					memset(piY, 0, iFilesCount_S30 * sizeof(double));
					iAllCount = 0;
					for (int iL = 0; iL < iFilesCount_S30; iL++)
					{
						if (piS30[iL] <= 0.0 || piL30inter[iL] <= 0.0) continue; //Invalid screening

						double tempX = (piS30[iL] - paraBandAdjust_intercept_V14A[iB]) / paraBandAdjust_slope_V14A[iB]; //Original 	
						double tempY = piL30inter[iL];
						if (tempX <= 0) continue; //Invalid original screening
						//Use the blue band reflectance to screen remaining clouds (Zhang et al., 2018)
						if (iB == 0) //Blue band
						{							
							if (fabs(tempX - tempY) > 0.5*(tempX + tempY)) piS30_QA[iL][iP] = 2;
						}
						if (piS30_QA[iL][iP] == 2) continue; //QA screening

						piX[iAllCount] = tempX;
						piY[iAllCount] = tempY;  //fitting Y
						iAllCount++;  //fitting count
					}//iL
				}

				if (iAllCount >= 4) pOutFittingQA[iP] = 2;  // indicate  linear interpolation matching
				if (iAllCount < 4)    // linear interpolation matching < 4, invalid point
				{
					pOutFittingQA[iP] = 3;
					delete piX;
					delete piY;
					delete piL30inter;
					delete piL30;
					delete piS30;
					continue;
				}


				/***** calculate the regression parameter *****/
				double R_square, slope, intercept;
				R_square = regression(piX, piY, iAllCount, &slope, &intercept);
				if (slope <= 0) { slope = 9999; intercept = 9999; }
				pOutSlope[iP] = slope;
				pOutIntercept[iP] = intercept;
				pOutFittingCount[iP] = int16(iAllCount);

				delete piX;
				delete piY;
				delete piL30inter;
				delete piL30;
				delete piS30;				
			} // iP
			
			//------------------------------------------------------------------------------
			//-------- Apply 3*3 Window Fitting for the failed pixels --------
			//------------------------------------------------------------------------------
			int winSize=1, winPixelCount = 9;
			for (int i = 0; i < iSrcDim[0]; i++)
			{
				for (int j = 0; j < iSrcDim[1]; j++)
				{
					int ijP = i * iSrcDim[1] + j;
					if (ijP % (iSrcDim[0] * iSrcDim[1] / 100 - 1) == 0) printf("*%d*", ijP / (iSrcDim[0] * iSrcDim[1] / 100 - 1));
					if (pOutSlope[ijP] != 9999) continue;

					double *piWFittingX = (double *)malloc(iTotalCount *winPixelCount * sizeof(double));
					memset(piWFittingX, 0, sizeof(double) * iFilesCount_S30*winPixelCount);
					double *piWFittingY = (double *)malloc(iTotalCount *winPixelCount * sizeof(double));
					memset(piWFittingY, 0, sizeof(double) * iFilesCount_S30*winPixelCount);
					int WFittingCount = 0;

					// Set the moving window, If the window extent is beyond the tile extent, clip it to the tile extent.
					int win_row_start = i - winSize;
					int win_row_end = i + winSize;
					int win_col_start = j - winSize;
					int win_col_end = j + winSize;
					if (i - winSize <= 0) win_row_start = 0;
					if (j - winSize <= 0) win_col_start = 0;
					if (i + winSize >= iSrcDim[0] - 1) win_row_end = iSrcDim[0] - 1;
					if (j + winSize >= iSrcDim[1] - 1) win_col_end = iSrcDim[1] - 1;

					//Select valid time series
					int iFittingCount = 0;
					for (int m = win_row_start; m < win_row_end; m++)
					{
						for (int n = win_col_start; n < win_col_end; n++)
						{
							int iP = m * iSrcDim[1] + n;
							//if (iP % (iSrcDim[0] * iSrcDim[1] / 100 - 1) == 0) printf("*%d*", iP / (iSrcDim[0] * iSrcDim[1] / 100 - 1));

							// Extract Landsat-8 time series of a pixel
				           // ------ Use a daily array to save Landsat-8 time series
							double *piL30_w = (double *)malloc(iTotalCount * sizeof(double));
							memset(piL30_w, 0, iTotalCount * sizeof(double));
							for (int iL = 0; iL < iFilesCount_L30; iL++)
							{
								//Invalid screening of Landsat-8
								if (piL30_QA[iL][iP] != 1 || piL30_Band[iL][iP] <= 0 || piL30_Band[iL][iP] > 10000) continue;
								int tmpIndex_L30 = (piDataYear[iL] - 2015) * 365 + (piDataDate[iL] - 1)
									+ int((piDataYear[iL] - 2015 + 2) / 4);
								piL30_w[tmpIndex_L30] = piL30_Band[iL][iP] * 0.0001;
							}

							// Extract Sentinel-2 time series of a pixel
							double *piS30_w = (double *)malloc(iFilesCount_S30 * sizeof(double));
							memset(piS30_w, 0, iFilesCount_S30 * sizeof(double));
							for (int iL = 0; iL < iFilesCount_S30; iL++)
							{
								//Invalid screening of Sentinel-2
								if (piS30_QA[iL][iP] != 1 || piS30_Band[iL][iP] <= 0 || piS30_Band[iL][iP] > 10000) continue;
								piS30_w[iL] = piS30_Band[iL][iP] * 0.0001;
							}

							//----------------------------------------------------------------------
							//-------- Linear matching: first 1D and then interpolation -------
							//-----------------------------------------------------------------------
							// Array used for linear matching (size = iFilesCount_S30)
							double *piL30inter_w = new double[iFilesCount_S30];
							memset(piL30inter_w, 0, iFilesCount_S30 * sizeof(double));
							for (int iL = 0; iL < iFilesCount_S30; iL++)
							{
								int tmpIndex_S30 = (piDataYear_S30[iL] - 2015) * 365 + (piDataDate_S30[iL] - 1)
									+ int((piDataYear_S30[iL] - 2015 + 2) / 4);
								//printf("%d\n%d\n", iL, tmpIndex_S30);

								// Find valid L30 observation within 1 days of S30
								if (piL30_w[tmpIndex_S30] > 0) { piL30inter_w[iL] = piL30_w[tmpIndex_S30]; continue; }
								if (piL30_w[tmpIndex_S30 - 1] > 0) { piL30inter_w[iL] = piL30_w[tmpIndex_S30 - 1]; continue; }
								if (piL30_w[tmpIndex_S30 + 1] > 0) { piL30inter_w[iL] = piL30_w[tmpIndex_S30 + 1]; continue; }

								// Find valid L30 observation within 16 days of S30 in backward direction
								int iStep = 2;
								double interLeft = 0.0;
								int interLeft_pos = 0;
								while (iStep < 16)
								{
									if (piL30_w[tmpIndex_S30 - iStep] > 0)
									{
										interLeft = piL30_w[tmpIndex_S30 - iStep];
										interLeft_pos = tmpIndex_S30 - iStep;
										break;
									}
									iStep = iStep + 1;
								} // step left

								// Find valid L30 observation within 16 days of S30 in forward direction
								iStep = 2;
								double interRight = 0.0;
								int interRight_pos = 0;
								while (iStep < 16)
								{
									if (piL30_w[tmpIndex_S30 + iStep] > 0)
									{
										interRight = piL30_w[tmpIndex_S30 + iStep];
										interRight_pos = tmpIndex_S30 + iStep;
										break;
									}
									iStep = iStep + 1;
								} // step right

								// Linear interpolation
								if (interLeft != 0.0 && interRight != 0.0)
								{
									piL30inter_w[iL] = interLeft + (interRight - interLeft) * (tmpIndex_S30 - interLeft_pos) / (interRight_pos - interLeft_pos);
								} // linear interpolation
							}

							for (int iL = 0; iL < iFilesCount_S30; iL++)
							{
								if (piS30_w[iL] <= 0.0 || piL30inter_w[iL] <= 0.0) continue; //Invalid screening

								double tempX = (piS30_w[iL] - paraBandAdjust_intercept_V14A[iB]) / paraBandAdjust_slope_V14A[iB]; //Original 	
								double tempY = piL30inter_w[iL];
								if (tempX <= 0) continue; //Invalid original screening
								//Use the blue band reflectance to screen remaining clouds (Zhang et al., 2018)
								if (iB == 0) //Blue band
								{
									if (fabs(tempX - tempY) > 0.5*(tempX + tempY)) piS30_QA[iL][iP] = 2;
								}
								if (piS30_QA[iL][iP] == 2) continue; //QA screening

								piWFittingX[WFittingCount] = tempX;
								piWFittingY[WFittingCount] = tempY;  //fitting Y
								WFittingCount++;  //fitting count
							}//iL
							delete piL30inter_w;
							delete piL30_w;
							delete piS30_w;
						}//n
					}//m

					 /***** calculate the regression parameter *****/
					double R_square_w, slope_w, intercept_w;
					R_square_w = regression(piWFittingX, piWFittingY, WFittingCount, &slope_w, &intercept_w);
					int tempPos = i * iSrcDim[1] + j;
					if (slope_w <= 0) {
						slope_w = 9999; intercept_w = 9999; pOutFittingQA[tempPos] = 4;
					}
					//The output of linear fitting parameters					
					pOutSlope[tempPos] = slope_w;
					pOutIntercept[tempPos] = intercept_w;
					pOutFittingCount[tempPos] = WFittingCount;

					delete piWFittingX;
					delete piWFittingY;
				}//j
			}//i


			for (int i = 0; i < iFilesCount_S30; i++) { delete piS30_Band[i]; }
			delete piS30_Band;
			for (int i = 0; i < iFilesCount_L30; i++) { delete piL30_Band[i]; }
			delete piL30_Band;

			//----------------------------------------------------------------
			//-------- Write the linear fitting parameters -------
			//----------------------------------------------------------------
			CString sOutputFileName2;
			sOutputFileName2.Format("%s.FittingPara%s.hdf", outputFileName, Landsat8_outputBands[iB]);
			printf("Writing: %s\n", sOutputFileName2);
			int32 sd_idw = SDstart(sOutputFileName2, DFACC_WRITE);
			int32 sds_idw = SDselect(sd_idw, SDnametoindex(sd_idw, "Slope"));
			SDwritedata(sds_idw, start_p, NULL, iSrcDim, pOutSlope);
			SDendaccess(sds_idw);
			sds_idw = SDselect(sd_idw, SDnametoindex(sd_idw, "Intercept"));
			SDwritedata(sds_idw, start_p, NULL, iSrcDim, pOutIntercept);
			SDendaccess(sds_idw);
			sds_idw = SDselect(sd_idw, SDnametoindex(sd_idw, "FittingCount"));
			SDwritedata(sds_idw, start_p, NULL, iSrcDim, pOutFittingCount);
			SDendaccess(sds_idw);
			sds_idw = SDselect(sd_idw, SDnametoindex(sd_idw, "FittingQA"));
			SDwritedata(sds_idw, start_p, NULL, iSrcDim, pOutFittingQA);
			SDendaccess(sds_idw);
			SDend(sd_idw);
		
			delete pOutSlope;
			delete pOutIntercept;
			delete pOutFittingCount;
		}//iB

		for (int i = 0; i < iFilesCount_S30; i++) { delete piS30_QA[i]; }
		delete piS30_QA; 
		for (int i = 0; i < iFilesCount_L30; i++) { delete piL30_QA[i]; }
		delete piL30_QA;
		printf("Part %d process ends!\n", iPart);
	}// iPart

	delete piDataYear_S30;
	delete piDataDate_S30;
	delete piDataYear;
	delete piDataDate;
	return;
}
