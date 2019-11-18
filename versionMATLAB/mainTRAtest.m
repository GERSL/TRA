%% This script shows an example of running TRA (By Rong Shang, 11/18/2019).

maindir = 'E:\TRA\'; % Change the dir here

% Read Landat 8 data
input_L30 = fullfile(maindir,'exampleData_L30.csv');
dataL30 = importdata(input_L30);
dataL30 = dataL30.data;
dataL30 = dataL30(find(dataL30(:,9)==1),:);

% Remove noise using a time series filter(threshold of 2000 can be changed).
idclr = TimeSeiresNoiseFilter(dataL30(:,3:8),2000);
dataL30 = dataL30(idclr,:);
srRef = dataL30(:,3:8);
dateRef = datenum(dataL30(:,1),1,0) + dataL30(:,2);


% Read Sentinel-2 data
input_S30 = fullfile(maindir,'exampleData_S30.csv');
dataS30 = importdata(input_S30);
dataS30 = dataS30.data;
dataS30 = dataS30(find(dataS30(:,9)==1),:);

% Remove noise using a time series filter(threshold of 2000 can be changed).
idclr = TimeSeiresNoiseFilter(dataS30(:,3:8),2000);
dataS30 = dataS30(idclr,:);
srTar = dataS30(:,3:8);
dateTar = datenum(dataS30(:,1),1,0) + dataS30(:,2);

% Match Landsat 8 and Sentinel-2
[matchRef, matchTar,interpFlag] = TRAmatch(srRef,dateRef,srTar,dateTar);

% Linear regression
[slope, intercept] = TRAregression(matchRef, matchTar);

% Adjust Sentinel-2
adjustTar = TRAadjust(srTar, slope, intercept);

