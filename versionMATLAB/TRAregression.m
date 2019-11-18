function [slope, intercept] = TRAregression(matchRef, matchTar)
%% TRAregression
%    This function conducts linear regression between the matched reference
% surface reflectance and the matched target surface reflectance. Developed
% by Rong Shang, 10/23/2019.
%
% Funtion input:
%    matchRef   Matched reference surface reflectance;
%    matchTar   Matched surface reflectance to be adjusted.
%
% Funtion output:
%    slope      The slope of linear regression for the six bands;
%    intercept  The intercept of linear regression for the six bands.

    slope = zeros(6,1);
    intercept = zeros(6,1);
    for i_B =1:6
        p = polyfit(matchTar(:,i_B),matchRef(:,i_B),1);
        slope(i_B) = p(1);
        intercept(i_B) = p(2);
    end
end