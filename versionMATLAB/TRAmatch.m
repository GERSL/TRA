function [matchRef, matchTar,interpFlag] = TRAmatch(srRef,dateRef,srTar,dateTar)
%% TRAmatch
%    This function matches the reference surface reflectance and the target
% surface reflectance using the one-day matching method. If there are less
% than 4 pairs of matched observations, the interpolation matching method
% will be used. Developed by Rong Shang, 10/23/2019.
%
% Funtion input:
%    srRef      Reference surface reflectance;
%    dateRef    The date of reference surface reflectance;
%    srTar      Surface reflectance to be adjusted;
%    dateTar    The date of surface reflectance to be adjusted.
%
% Funtion output:
%    matchRef   Matched reference surface reflectance;
%    matchTar   Matched surface reflectance to be adjusted;
%    interpFlag Flag of whether using linear interpolation(1) or not(0).

    
    % ------ One-day match ------ 
    interpFlag = 0;
    lengthTar = length(dateTar); 
    matchRef = zeros(lengthTar,6);
    matchTar = zeros(lengthTar,6);         
    numMatch = 0;
    for i = 1:lengthTar
        curDateTar = dateTar(i);
        % Same day match
        index = find(dateRef == curDateTar);
        if length(index) == 1
            numMatch = numMatch + 1;
            matchTar(numMatch,:) = srTar(i,:);
            matchRef(numMatch,:) = srRef(index,:);            
            continue;
        end
        % Previoud one day match
        index = find(dateRef == curDateTar - 1);
        if length(index) == 1
            numMatch = numMatch + 1;
            matchTar(numMatch,:) = srTar(i,:);
            matchRef(numMatch,:) = srRef(index,:);            
            continue;
        end
        % Next one day match
        index = find(dateRef == curDateTar + 1);
        if length(index) == 1
            numMatch = numMatch + 1;
            matchTar(numMatch,:) = srTar(i,:);
            matchRef(numMatch,:) = srRef(index,:);            
            continue;
        end   
    end
    
    % ------ Interpolation match ------
    if numMatch < 4
        interpFlag = 1;
        matchRef = zeros(lengthTar,6);
        matchTar = zeros(lengthTar,6);         
        numMatch = 0;        
        for i = 1:lengthTar
            curDateTar = dateTar(i);
            curInterp = zeros(1,6);
            
            % Left part
            index_left = find(dateRef < curDateTar);
            index_right = find(dateRef > curDateTar);
            preDate = dateRef(index_left);
            aftDate = dateRef(index_right);
            
            if length(preDate) < 1 || length(aftDate) < 1  
                continue;
            end
            preDate = preDate(end);           
            aftDate = aftDate(1);
            
            if abs(preDate - curDateTar) <= 16 && abs(aftDate - curDateTar) <= 16              
                preRef = srRef(index_left,:);               
                preRef = preRef(end,:);              
                aftRef = srRef(index_right,:);                
                aftRef = aftRef(1,:);
                
                for i_B = 1:6       
                    curInterp(1,i_B) = preRef(i_B) + ...
                        (aftRef(i_B) - preRef(i_B))*(curDateTar - preDate)/(aftDate - preDate);
                end
                
                numMatch = numMatch + 1;
                matchTar(numMatch,:) = srTar(i,:);
                matchRef(numMatch,:) = curInterp(1,:);            
            end
        end
    end
    
    matchTar = matchTar(find(matchTar(:,1)>0),:);
    matchRef = matchRef(find(matchRef(:,1)>0),:);
end