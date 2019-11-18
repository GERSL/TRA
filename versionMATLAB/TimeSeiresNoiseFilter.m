function idclr = TimeSeiresNoiseFilter(clry,threshold)
%% TimeSeiresNoiseFilter
%    This function returns the index of clear-sky observations. Notice that
% the first and last observations will be regarded as clear observations.
% Developed by Rong Shang, 11/18/2019.
%
% input: surface reflectance of six spectral bands
% output: the index of clear-sky observations


    new_clry = clry(:,:);    
    mean_clry = (clry(1:end-2,:)+clry(3:end,:))*0.5;    
    new_clry(2:end-1,:) = mean_clry;
    delta_clry = abs(new_clry - clry(:,:));    
    delta_clry = delta_clry';
    
    sumBands = sum(delta_clry(:,:));
    sumBands = sumBands';
    idclr = sumBands(:) < threshold;    
end