function TRAplot(srRef,dateRef,srTar,dateTar,adjustTar,PreFlag,pid)
%% TRAplot
%   This function plots the two kinds of surface reflectance (original or 
% adjusted). Developed by Rong Shang, 10/23/2019.
%
% Funtion input:
%    srRef      Reference surface reflectance;
%    dateRef    The date of reference surface reflectance;
%    srTar      Surface reflectance to be adjusted;
%    dateTar    The date of surface reflectance to be adjusted.
%    adjustTar  The adjusted surface reflectance for the six bands;
%    PreFlag    Plot original or adjusted;
%    pid        The ID of the point.

    % Set the figure
    fig_lts = figure('units','normalized','outerposition',[0 0 1 1]);
    ha = tight_subplot(3,2,[.06 .05],[.04 .04],[.05 .02]);
    bandNames = {'Blue', 'Green', 'Red', 'NIR', 'SWIR1', 'SWIR2'};
    % Plot figure
    for iiii =1:6
        axes(ha(iiii)); 
        maxFigureLength = 3000;
        tStep=500;
        sr_band = iiii;
        if sr_band > 3
            maxFigureLength = 5000;
            tStep = 1000;
        end

        hold on;
        ylim([0, maxFigureLength]);        
        plot(dateRef(:), srRef(:,sr_band),'.b');
        if PreFlag == 1
            plot(dateTar, srTar(:,sr_band),'.r');
        else
            plot(dateTar, adjustTar(:,sr_band),'.r');
        end
         

       
        set(gca,'YTick', [0:tStep:maxFigureLength]);
        set(gca,'YtickLabel',[0:tStep:maxFigureLength]);
        xlim([datenum(1984,1,0),datenum(2018,1,0)]);
        datetick('x','yyyy','keeplimits');           
        y_ti = bandNames(sr_band);
        ylabel(y_ti);
        
        if sr_band == 1
            if PreFlag == 1
                legend('PathRow1','PathRow2');
            else
                legend('PathRow1','AdjustPR2');
            end
        end
    end
    
    if PreFlag == 1
        filename = 'BRDF';
    else
        filename = 'TRA';
    end
    output_Path = char(strcat('E:\CurrentWork\RecursiveCCD\COLD_v13_03_FF\OverlapFigure\P',...
        int2str(pid),'_',filename,'.png'));
    saveas(fig_lts,output_Path,'png');
    close(fig_lts);
end