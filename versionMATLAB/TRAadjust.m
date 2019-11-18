function adjustTar = TRAadjust(srTar, slope, intercept)
%% TRAadjust
%     This function adjusts the target surface reflectance using the built
% linear relationship. Developed by Rong Shang, 10/23/2019.
%
% Funtion input:
%    srTar      Original surface reflectance to be adjusted;
%    slope      The slope of linear regression for the six bands;
%    intercept  The intercept of linear regression for the six bands.
%
% Funtion output:
%    adjustTar  The adjusted surface reflectance for the six bands.


    adjustTar = zeros(size(srTar,1),6);
    for i_B =1:6
        adjustTar(:,i_B) = srTar(:,i_B)*slope(i_B) + intercept(i_B);
    end
end