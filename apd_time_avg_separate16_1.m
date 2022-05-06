function apd_time_avg_separate16_1(method,num,backorder,backs,backmeth)
%Default function call:     apd_time_avg_separate16_1('amp',20)
%   background = backcor(time,averagepixel,backorder,backs,backmeth);

% possible method inputs: deriv(does derivative method)
% amp (does % of amplitude method)
% possible num inputs: number desired for either derivative (ie 1 or 2)
% or amplitude (ie 10, 20, etc) method.

%Developed by Nate Huebsch, Daniel Simmons, and David Schuftan based on code by Laughner
%et al. 2012

if exist('method', 'var') == 1
    
else
    method = 'amp';
end

if exist('num', 'var') == 1
    
else
    num = 20;
end

if exist('backorder','var') == 1

else
    backorder = 2;
end

if exist('backs','var') == 1

else 
    backs = 0.01;
end

if exist('backmeth','var')==1
    
else
    backmeth = 'atq';
end


[files,pathname] = uigetfile('*.*', 'Select One or More Files' , 'MultiSelect', 'on');

if isa(files, 'char') == 1
    total = 1;
else
    total = max(size(files));
end

if ispc
    delim = '\';
elseif ismac
    delim = '/';
end

count = 1;

while count <= total
    fprintf('\n %d of %d \n', count, total)
    if isa(files, 'char') == 1
        filename = files;
    else
        filename = files{1,count};
    end
    disp(filename);
    myVars = {'cmosData','frequency'};
    load([pathname filename],myVars{:});
    
    filename = filename(1:end-4);
    cmosData = cmosData(:,:,20:end-20);
    
    firstFrame = cmosData(:,:,1);
    [x,y] = find(firstFrame);
    top = min(x);
    bottom = max(x);
    left = min(y);
    right = max(y);
    
    if top == 1 && left == 1 && right == size(cmosData,2) && bottom == size (cmosData,1)
        cmosData = reshape(cmosData,[],size(cmosData,3));
    else
        dummy = zeros(length(x),size(cmosData,3));
        for i = 1:size(cmosData,3)
            for j = 1:length(x)
                dummy(j,i) = cmosData(x(j),y(j),i);
            end
        end
        cmosData = dummy;
    end
    
    
    
    
    [~,p] = size(cmosData);
    
    % Get average RFU for data
    averagepixel = zeros(1,p);
    
    for jj = 1:p
        
        dummy = cmosData(:,jj);
        averagepixel(jj) = mean(dummy(:));
        
    end

    averagepixel_orig = averagepixel; % Save unaltered original trace for exporting
    
    if frequency == 0 || isnan(frequency) || ~exist('frequency','var')
        framerate = str2double(cell2mat(inputdlg('Framerate was not found. Please input Framerate (frames/sec)')));
    else
        framerate = frequency;
    end
    
    % Median Filter
    averagepixel = medfilt1(averagepixel,5);
    
    %lowpass filter
%     buff1 = repmat(averagepixel(1),1,50);
%     buff2 = repmat(averagepixel(length(averagepixel)),1,50);
%     averagepixel1 = [buff1 averagepixel buff2];
%     
%     if framerate/2 <= 100
%         averagepixel1 = lowpass(averagepixel1,(framerate/2)-1,framerate);
%     else
%         averagepixel1 = lowpass(averagepixel1,100,framerate);
%     end
%     
%     averagepixel = averagepixel1(51:length(averagepixel1)-50);
    
    [~,p] = size(averagepixel);
    
    time = (1:p)/framerate;
    
    fig = figure;
    % Determine baseline
   % background = backcor(time,averagepixel,2,0.01,'atq');
    background = backcor(time,averagepixel,backorder,backs,backmeth);
    % plot raw BeRST signal
    subplot(3,2,(1:2));
    plot(time,averagepixel,'b')
    xlabel('Time (sec)'), ylabel({'BeRST','Intensity (RFU)'}), hold on
    % plot background
    plot(time, background,'r'), hold off
    title([filename ' ' method num2str(num)],'Interpreter','none');
    
    %% Try
    try
        % Remove background
        
        dummy = averagepixel - background';
        
        if numel(dummy) == numel(averagepixel)
            averagepixel = dummy;
        else
            averagepixel = averagepixel - background;
        end
        
        
        
        range = max(averagepixel(:)) - min(averagepixel(:));
        
        %plot berst with background removed
        subplot(3,2,(3:4)); plot(time,averagepixel,'b'),...
            xlabel('Time (sec)'), ylabel({'Corrected BeRST','Intensity (RFU)'}),...
            ylim([min(averagepixel(:)) (max(averagepixel(:))+0.1*range)]);
        hold on
        
        %%

        % Find individual waveforms based on extrema
        invertdata = abs((averagepixel/max(averagepixel))-1)*max(averagepixel);
        [~,maxs] = findpeaks(averagepixel,'MinPeakHeight',0.6*max(averagepixel)); % Find all local maxima
        [~,mins] = findpeaks(invertdata,'MinPeakHeight',0.4*max(invertdata)); % Find all local minima
        
        %Filter false maxima and minima
        
        extrema = sort([maxs mins]); %all unfiltered maximum and minimum local extrema
        
        maxi = zeros(length(maxs),1);
        for i= 1:length(maxs)
            maxi(i) = find(extrema==maxs(i)); % indices of local max indices within extrema array
        end
        
        mini = zeros(length(mins),1);
        for i = 1:length(mins)
            mini(i) = find(extrema==mins(i)); % indices of local min indices within extrema array
        end
        
        
        % Find true local maxima of each AP and use as peak
        keepmax = zeros(length(maxi));
        k = 1;
        
        
        if length(maxi) == 1 % If only one local maximum is found use as lone peak
            keepmax(k,1) = extrema(maxi);
            
        elseif isempty(maxi) % If no local maximum are found throw error
            figure
            plot(averagepixel)
            error('No Local Maxima Found')
            
        else % If multiple local maximum are found, make sure it is the true maximum with a local minimum the adjacent extrema
            i=1;
            while i<=length(maxi)
                good = zeros(length(maxi)+1,1);
                j = 1;
                
                if  i == 1 && maxi(i)+1 ~= maxi(i+1) % If local max is first extrema and next local max isn't second extrema
                    good(j) = maxi(i);
                elseif i == length(maxi) && maxi(i)-1 ~= maxi(i-1) % If local max is last extrema and previous extrema isn't a local max
                    good(j) = maxi(i);
                elseif i < length(maxi) &&  maxi(i)+1 == maxi(i+1) % If local max isn't last extrema and next extrema is also a local max
                    while i<length(maxi) &&  maxi(i)+1 == maxi(i+1)
                        good(j) = maxi(i);
                        good(j+1)=maxi(i+1);
                        j=j+1;
                        i=i+1;
                    end
                elseif i > 1 && i < length(maxi) && maxi(i)-1~=maxi(i-1) && maxi(i)+1 ~= maxi(i+1) % If local max isn't first or last extrema and there are no adjacent local max
                    good(j) = maxi(i);
                end
                
                good(good == 0) = [];
                [~,dummy] = max(averagepixel(extrema(good))); % Select maximum of the adjacent local max or use the local max if no adjacent local max
                keepmax(k,1) = extrema(good(dummy));
                k = k+1;
                i=i+1;
                
            end
            keepmax(keepmax==0) = [];
            keepmax = keepmax';
        end
        
        
        
        % Find true local minima of each AP peak
        keepmin = zeros(length(mini));
        k = 1;
        if length(mini) == 1% If only one local minimum is found use as lone peak
            keepmin(k,1) = extrema(mini);
            
        elseif isempty(mini) % If no local minimum are found throw error
            figure
            plot(averagepixel)
            error('No Local Minima Found')
            
        else % If multiple local minimum are found, make sure it is the true minimum with a local maximum the adjacent extrema
            i=1;
            while i<=length(mini)
                good = zeros(length(mini)+1,1);
                j = 1;
                
                if  i == 1 && mini(i)+1 ~= mini(i+1) % If local min is first extrema and next local min isn't second extrema
                    good(j) = mini(i);
                elseif i == length(mini) && mini(i)-1 ~= mini(i-1) % If local min is last extrema and previous extrema isn't a local min
                    good(j) = mini(i);
                elseif i < length(mini) &&  mini(i)+1 == mini(i+1)  % If local min isn't last extrema and next extrema is also a local min
                    while i<length(mini) &&  mini(i)+1 == mini(i+1)
                        good(j) = mini(i);
                        good(j+1)=mini(i+1);
                        j=j+1;
                        i=i+1;
                    end
                elseif i > 1 && i < length(mini) && mini(i)-1~=mini(i-1) && mini(i)+1 ~= mini(i+1) % If local min isn't first or last extrema and there are no adjacent local min
                    good(j) = mini(i);
                end
                good(good == 0) = [];
                [~,dummy] = min(averagepixel(extrema(good))); % Select maximum of the adjacent local min or use the local min if no adjacent local min
                keepmin(k,1) = extrema(good(dummy));
                k = k+1;
                i=i+1;
                
            end
            keepmin(keepmin == 0) = [];
            keepmin = keepmin';
        end
        
        keepextrema = [keepmax' keepmin'];
        keepextrema = sort(keepextrema)';
        if keepmax(1) == keepextrema(1)
            keepmax(1) = 0;
            keepextrema(1)= 0;
        end
        if keepmax(length(keepmax)) == keepextrema(length(keepextrema))
            keepmax(length(keepmax)) = 0;
        end
        
        keepmax((~keepmax)) = [];
        
        scatter(time(keepmax),averagepixel(keepmax),'k*');
        
       
        chopstart = zeros(length(keepmax),1);
        chopend = zeros(length(keepmax),1);
        totalsize = zeros(length(keepmax),1);
        for i = 1:length(keepmax)
            chopstart(i) = keepmax(i) - keepmin(i);
            chopend(i) = keepmin(i+1) - keepmax(i);
            totalsize(i) = chopstart(i) + chopend(i)+1;
        end
        
      
        averagedata = zeros(max(totalsize),length(keepmax));
        for i=1:length(keepmax)
            averagedata(1:totalsize(i),i) = averagepixel(1,keepmax(i)-chopstart(i):keepmax(i)+chopend(i));
        end
        
        % Find Beatrate
        AverageDist=0;
        for i=2:length(keepmax)
            AverageDist = AverageDist+ keepmax(i)-keepmax(i-1);
        end
        AverageDist = AverageDist/(length(keepmax)-1);
        Beatrate = 1/AverageDist * framerate * 60;
        RR = 60 ./ Beatrate; % 60 (sec/min)/bpm
        
        
        %Define constants
        requiredVal90 = 0.1;
        requiredVal80 = 0.2;
        requiredVal50 = 0.5;
        requiredVal30 = 0.7;
        plateaucutoff = 0.8;
        
        % pre-set average waveform info
        APD90 = zeros(length(keepmax),1);
        APD80 = zeros(length(keepmax),1);
        APD50 = zeros(length(keepmax),1);
        APD30 = zeros(length(keepmax),1);
        UPD = zeros(length(keepmax),1);
        upstarts = zeros(length(keepmax),1);
        APD90t = zeros(length(keepmax),1);
        APD80t = zeros(length(keepmax),1);
        APDMAX = zeros(length(keepmax),1);
        dVdtmax = zeros(length(keepmax),1);
        plateauDuration = zeros(length(keepmax),1);
       

        for i = 1:length(keepmax)
            
            % Get single waveform
            tempdata = averagedata(:,i);
            dummy = find(tempdata, 1, 'last');
            tempdata = tempdata(1:dummy);
            tempdata3 = tempdata;
            times = (1:totalsize(i))/frequency;

            %Find maximum of the signal
            [APDMAX(i),maxiter] = max(tempdata);

            % normalize
            minimum = min(tempdata);
            maximum = max(tempdata);
            difference = maximum-minimum;
            tempdata = (tempdata-minimum)./difference; 
            
            %%Determining activation time point
            if strcmp('deriv',method)    % Find chosen Derivative and time of maximum
                apd_Upstroke = tempdata(1:round(maxiter),1);
                apd_data2 = diff(apd_Upstroke,num,1); 
                [~, maxderiter] = max(apd_data2,[],1);

            elseif strcmp('amp',method)  % Use given percentage above baseline as activation point
                [~,dummypeak] = max(tempdata);
                tempdata4 = tempdata(1:dummypeak);
                [~,maxderiter] = min(abs((num/100)-tempdata4));                
            end         
            
            % Plot upstroke start time in panel 2
            upstarts(i) = time(keepmax(i)-chopstart(i)+maxderiter-1);
            subplot(3,2,(3:4)); hold on;
            scatter(upstarts(i),averagepixel(keepmax(i)-chopstart(i)+maxderiter-1),'filled','r');
            



            %%Find dV/dt Max%%
            upstarting = keepmin(i);
            upending = keepmax(i);
            
            if upending-upstarting < 8 % interp function requires 9 or more points
                upstarting = upstarting - (8-(upending-upstarting));
            end
            
            %Check if 9 or more points are included and expand end of
            %window if more points neeeded
            shortChopFlag = 0;
            if upstarting <= 0
                upstarting = 1;
                upending = 9;
                shortChopFlag = 1;
            end
            
            upstrokeonly_t = time(upstarting:upending) ; upstrokeonly_t = upstrokeonly_t';
            upstrokeonly_amp = averagepixel(upstarting:upending);
            
            
            r = 5; %artifically upsample to 5x actual sampling frequency for linear interpolation
            
            %Interpolate Data to find more accurate location of dV/dt max
            tempdata5 = interp(upstrokeonly_amp,r);          
            averagetime5 = interp(upstrokeonly_t,r);  averagetime5 =  averagetime5';

            % Find time closest to start and peak if points added to allow
            % 9 points for interpolation and remove points from expanded
            % window
            if shortChopFlag == 1
                [~,upstrokeonly_start] = min(abs(sqrt((tempdata5-upstrokeonly_amp(keepmin(i))).^2+(averagetime5-upstrokeonly_t(keepmin(i))).^2)));
                [~,upstrokeonly_end] = min(abs(sqrt((tempdata5-upstrokeonly_amp(keepmax(i))).^2+(averagetime5-upstrokeonly_t(keepmax(i))).^2)));
                tempdata5 = tempdata5(upstrokeonly_start:upstrokeonly_end);
                averagetime5 = averagetime5(upstrokeonly_start:upstrokeonly_end);
            end
            
            timeincrement = (upstrokeonly_t(2) - upstrokeonly_t(1))/r; %Set new timeinterval based on upsampling of interpolated data
            
            
            %Low Pass Filter
            n = 4; % Filter order
            Fs = framerate * r; % Sampling frequency from above
            Fc = 50;
            [bLow, aLow] = butter(n, Fc / (Fs/2)); % Creates a low pass Butterworth filter
            tempdata5 = filter(bLow,aLow,tempdata5);
            
            dVdt = diff(tempdata5);  % dV/dt of all data
            [dVdtmax(i),~] = max(dVdt(:)); % value of maximum dV/dt of data
             
            maxt = dVdt == dVdtmax(i); % Find index of max dV/dt within all dV/dt
            maxtime = averagetime5(maxt) ; % Find time of dV/dt max within single AP
            dVdtmax(i) = dVdtmax(i) ./ timeincrement; % Scale dV/dt based on upsampling
            [~, closestIndex] = min(abs(time-maxtime)); % Find time index of specific AP dV/dt max within entire trace
            
            
            % Plot end points that define max dV/dt
            subplot(3,2,(3:4)); hold on;
            scatter(time(closestIndex),averagepixel(closestIndex),'c'); % Last point of max dV/dt
            scatter(time(closestIndex-1),averagepixel(closestIndex-1),'c*'); % First point of max dV/dt
            
            
            clear tempdata5; clear averagetime5; clear dVdt;
            
            % starting from the peak of the signal, loop until we reach value for APD90
            for k = maxiter:max(size(tempdata))
                if tempdata(k) <= requiredVal90
                    apd90iter = k; % Save the index when the baseline is reached
                    % this is the repolarization time point
                    break;
                end
                
                if k == length(tempdata)
                    apd90iter = 0;
                end
                
            end
            
            % plot APD90 end time in panel 2
            APD90t(i) = time((keepmax(i)-chopstart(i))+apd90iter);
            
            if apd90iter ~= 0
                scatter(APD90t(i),averagedata(apd90iter,i),'g');
            end
            
            % starting from the peak of the signal, loop until we reach value for APD80
            for k = maxiter:max(size(tempdata))
                if tempdata(k) <= requiredVal80
                    apd80iter = k; % Save the index when the baseline is reached
                    % this is the repolarization time point
                    break;
                end
                
                if k == length(tempdata)
                    apd80iter = 0;
                end
                
            end
            
            if apd90iter == 0
                if apd80iter ~= 0
                    APD80t(i) = time((keepmax(i) - chopstart(i)) + apd80iter);
                    scatter(APD80t(i), averagedata(apd80iter,i),'c');
                end
            end
            
            % starting from the peak of the signal, loop until we reach value for APD50
            for k = maxiter:max(size(tempdata))
                if tempdata(k) <= requiredVal50
                    apd50iter = k; % Save the index when the baseline is reached
                    % this is the repolarization time point
                    break;
                end
                
                if k == length(tempdata)
                    apd50iter = 0;
                end
                
            end
            
            % starting from the peak of the signal, loop until we reach value for APD30
            for k = maxiter:max(size(tempdata))
                if tempdata(k) <= requiredVal30
                    apd30iter = k; % Save the index when the baseline is reached
                    % this is the repolarization time point
                    break;
                end
                
                if k == length(tempdata)
                    apd30iter = 0;
                end
                
            end
            if (apd50iter == 0 || apd90iter == 0)
                disp([i ' Did not find correct apd.']);
            end
            
            
            % Calculate APD in frames between upstroke initiation and the point where reached the desired amplitude percentage
            diffIndex90 = apd90iter - maxderiter;
            diffIndex80 = apd80iter - maxderiter;
            diffIndex50 = apd50iter - maxderiter;
            diffIndex30 = apd30iter - maxderiter;
            
            % Calculate APD values in ms based on frames and frequency
            apd90 = diffIndex90*(1000/framerate);
            apd80 = diffIndex80*(1000/framerate);
            apd50 = diffIndex50*(1000/framerate);
            apd30 = diffIndex30*(1000/framerate);
            
            if apd90iter == 0
                APD90(i) = NaN;
            else
                APD90(i) = apd90;
            end
            
            if apd80iter == 0
                APD80(i) = NaN;
            else
                APD80(i) = apd80;
            end
            
            if apd50iter == 0
                APD50(i) = NaN;
            else
                APD50(i) = apd50;
            end
            
            if apd30iter == 0
                APD30(i) = NaN;
            else
                APD30(i) = apd30;
            end
            
            % Calculate Upstroke Duration
            UPD(i) = (maxiter - maxderiter)*(1000/framerate);
            
                        
            % Determine width of peak plateau using time the AP is above 80% of the signal peak.
            
            for k = maxderiter:maxiter
                if tempdata(k) >= plateaucutoff
                    startplateau = k; % Save index when start of plateau is reached
                    break;
                end
            end
            
            for k = maxiter:apd90iter
                if tempdata(k) <=plateaucutoff
                    endplateau = k; % Save index when end of plateau is reached
                    break;
                end
            end
            
            % Calculate the duration of the peak plateau in frames
            plateauFrames = endplateau - startplateau;
            
            % Calculate the duration of the peak plateau in ms based on frames and frequency
            plateauDuration(i) = plateauFrames*(1000/framerate);
            
            




            % Plot individual waveform w/ location markers
            subplot(3,2,5)
            plot(times,tempdata3,'b');
            ylabel({'Corrected BeRST','Intensity (RFU)'}), xlabel('Time (sec)')
            ylim([min(averagedata(:)) (max(averagedata(:))+0.1*range)]);
            hold on
            % upstroke marker
            scatter(times(maxderiter),tempdata3(maxderiter),'r');
            % apd90 marker
            if apd90iter ~= 0
                scatter(times(apd90iter),tempdata3(apd90iter),'g')
            end
            % apd50 marker
            if apd50iter ~= 0
                scatter(times(apd50iter),tempdata3(apd50iter),'k')
            end
            % apd30 marker
            if apd30iter ~= 0
                scatter(times(apd30iter),tempdata3(apd30iter),'b')
            end
            % max marker
            scatter(times(maxiter),tempdata3(maxiter),'m')
        end % End of for loop for individual wave forms
        
        F0 = background(keepmin(1)); % Calcium Background Intensity
        
        %% Get Avg & Std Values
        avgMax = round(mean(APDMAX));
        stdMax = round(std(APDMAX));
        avg90 = round(mean(APD90,'omitnan'));
        std90 = round(std(APD90,'omitnan'));
        avg80 = round(mean(APD80,'omitnan'));
        std80 = round(std(APD80,'omitnan'));
        avg50 = round(mean(APD50,'omitnan'));
        std50 = round(std(APD50,'omitnan'));
        avg30 = round(mean(APD30,'omitnan'));
        std30 = round(std(APD30,'omitnan'));
        avgUp = round(mean(UPD));
        stdUp = round(std(UPD));
        avgSize = round(mean(totalsize));
        stdSize = round(std(totalsize));
        DelF_F0 = avgMax/F0;
        dV_dt = mean(dVdtmax);
        std_dV_dt = std(dVdtmax);
        peakDuration = round(mean(plateauDuration,'omitnan'));
        std_peakDuration = round(std(plateauDuration,'omitnan'));
        
        % Get Fridericia corrected output values
        coravg90 = round(avg90 ./ RR .^ (1/3));
        corstd90 = round(std90 ./ RR .^ (1/3));
        coravg80 = round(avg80 ./ RR .^ (1/3));
        corstd80 = round(std80 ./ RR .^ (1/3));
        coravg50 = round(avg50 ./ RR .^ (1/3));
        corstd50 = round(std50 ./ RR .^ (1/3));
        coravg30 = round(avg30 ./ RR .^ (1/3));
        corstd30 = round(std30 ./ RR .^ (1/3));
        coravgUp = round(avgUp ./ RR .^ (1/3));
        corstdUp = round(stdUp ./ RR .^ (1/3));
        coravgSize = round(avgSize ./ RR .^ (1/3));
        corstdSize = round(stdSize ./ RR .^ (1/3));
        corpeakDuration = round(peakDuration ./ RR .^ (1/3));
        corstd_peakDuration = round(std_peakDuration ./ RR .^ (1/3));
        
        
        % Plot output numbers
        subplot(3,2,6)
        axis([0, max(time(:)),0, max(averagedata(:))]);
        set(gca, 'XTick', []); set(gca,'Ytick',[]);
        
        text(double(0.05*max(time(:))),avgMax,sprintf('AP Amp: %.0f +- %.0f RFU',avgMax, stdMax), 'FontSize', 8);
        text(double(0.05*max(time(:))),0.85*avgMax,['APD_{30}:' num2str(coravg30) ' \pm ' num2str(corstd30) 'ms'], 'FontSize', 8);
        text(double(0.05*max(time(:))),0.7*avgMax,['APD_{80}:' num2str(coravg80) ' \pm ' num2str(corstd80) 'ms'], 'FontSize', 8);
        text(double(0.05*max(time(:))),0.55*avgMax,['APD_{90}:' num2str(coravg90) ' \pm ' num2str(corstd90) 'ms'], 'FontSize', 8);
        text(double(0.05*max(time(:))),0.4*avgMax,['UPD:' num2str(coravgUp) ' \pm ' num2str(corstdUp) 'ms'], 'FontSize', 8);
        text(double(0.05*max(time(:))),0.25*avgMax,sprintf('Beat Rate: %.0f bpm',Beatrate), 'FontSize', 8);
        text(double(0.05*max(time(:))),.1*avgMax,['\DeltaF/F_{0}: ' num2str(DelF_F0)], 'FontSize', 8);
        text(double(0.05*max(time(:))),-.05*avgMax,['dV/dt:' num2str(dV_dt) ' \pm ' num2str(std_dV_dt) ' RFU/s'], 'FontSize', 8);
        
        count = count + 1;
        
        
        dirname = pathname;
        
        if isequal(exist([dirname  filename], 'dir'),0)
            mkdir([dirname  filename])
            dirname = ([dirname filename delim]);
        else
            dirname = ([dirname filename delim]);
        end
        
        filenamemethnum = [filename '_' method num2str(num)];
        
        
        % Save background & average
        save([dirname delim filenamemethnum '_background.mat'],'background','time');
        save([dirname delim filenamemethnum '_originaltrace.mat'],'averagepixel_orig','time');
        save([dirname delim filenamemethnum '_correctedtrace.mat'],'averagepixel','time');
        
        % Save APD Data
        save([dirname delim filenamemethnum '_APD_Data.mat'],'avgMax','stdMax', 'Beatrate',...
            'avg90','std90', 'avg80', 'std80', 'avg50', 'std50','avg30',...
            'std30','avgUp','avgSize','stdSize','stdUp','DelF_F0','F0','dV_dt','std_dV_dt',...
            'peakDuration','std_peakDuration','coravg90','corstd90','coravg80','corstd80',...
            'coravg50','corstd50','coravg30','corstd30','coravgUp','corstdUp',...
            'coravgSize','corstdSize','corpeakDuration','corstd_peakDuration');        
        
        % Save Figure as Matlab File
        savefig(fig,[dirname delim filenamemethnum  '.fig']);
        
        % New Folder to Save APD Data Only
        
        if isequal(exist([pathname 'APD Data'], 'dir'),0)
            mkdir([pathname 'APD Data'])
        end
        
        save([pathname 'APD Data' delim filenamemethnum  '_APD_Data.mat'],'avgMax','stdMax', 'Beatrate',...
            'avg90','std90', 'avg80', 'std80', 'avg50', 'std50','avg30',...
            'std30','avgUp','avgSize','stdSize','stdUp','DelF_F0','F0','dV_dt','std_dV_dt',...
            'peakDuration','std_peakDuration','coravg90','corstd90','coravg80','corstd80',...
            'coravg50','corstd50','coravg30','corstd30','coravgUp','corstdUp',...
            'coravgSize','corstdSize','corpeakDuration','corstd_peakDuration');
        
        fprintf('Calculated Upstroke Duration: %.2f +- %.2f ms\n',avgUp,stdUp);
        fprintf('Calculated APD90: %.2f +- %.2f ms\n',avg90,std90);
        
        
        
        %% Catch      
    catch
        fprintf(['Unable to Fully Process File '   filename]);
        count = count+1;
        
        if exist('Beatrate','var') == 0
            Beatrate = NaN;
        end
        
        dV_dt = NaN;
        std_dV_dt = NaN;
        avgMax = NaN;
        stdMax = NaN;
        avg90 = NaN;
        std90 = NaN;
        avg80 = NaN;
        std80 = NaN;
        avg50 = NaN;
        std50 = NaN;
        avg30 = NaN;
        std30 = NaN;
        avgUp = NaN;
        stdUp = NaN;
        avgSize = NaN;
        stdSize = NaN;
        DelF_F0 = NaN;
        F0 = NaN;
        peakDuration = NaN;
        std_peakDuration = NaN;
        
        coravg90 = NaN;
        corstd90 = NaN;
        coravg80 = NaN;
        corstd80 = NaN;
        coravg50 = NaN;
        corstd50 = NaN;
        coravg30 = NaN;
        corstd30 = NaN;
        coravgUp = NaN;
        corstdUp = NaN;
        coravgSize = NaN;
        corstdSize = NaN;
        corpeakDuration = NaN;
        corstd_peakDuration = NaN;
        
        newdirname = pathname;
        filenamenul = ['nul_' filename];
        if isequal(exist([newdirname filenamenul], 'dir'),0)
            mkdir([newdirname filenamenul])
            newdirname = ([newdirname filenamenul delim]);
        else
            newdirname = ([newdirname filenamenul delim]);
        end
        
        filenamenulmethnum = [filenamenul '_' method num2str(num)];
        
        save([newdirname delim filenamenulmethnum '_APD_Data.mat'],'avgMax','stdMax', 'Beatrate',...
            'avg90','std90', 'avg80', 'std80', 'avg50', 'std50','avg30',...
            'std30','avgUp','avgSize','stdSize','stdUp','DelF_F0','F0','dV_dt','std_dV_dt',...
            'peakDuration','std_peakDuration','coravg90','corstd90','coravg80','corstd80',...
            'coravg50','corstd50','coravg30','corstd30','coravgUp','corstdUp',...
            'coravgSize','corstdSize','corpeakDuration','corstd_peakDuration');

        save([newdirname delim filenamenulmethnum '_background.mat'],'background','time');
        save([newdirname delim filenamenulmethnum '_originaltrace.mat'],'averagepixel_orig','time');
        save([newdirname delim filenamenulmethnum '_correctedtrace.mat'],'averagepixel','time');
        
        % Save Figure as Matlab File
        savefig(fig,[newdirname delim filenamenulmethnum '.fig']);
        
        
        % New Folder to Save APD Data Only
        if isequal(exist([pathname 'APD Data'], 'dir'),0)
            mkdir([pathname 'APD Data'])
        end
        
        save([pathname 'APD Data' delim filenamenulmethnum  '_APD_Data.mat'],'avgMax','stdMax', 'Beatrate',...
            'avg90','std90', 'avg80', 'std80', 'avg50', 'std50','avg30',...
            'std30','avgUp','avgSize','stdSize','stdUp','DelF_F0','F0','dV_dt','std_dV_dt',...
            'peakDuration','std_peakDuration','coravg90','corstd90','coravg80','corstd80',...
            'coravg50','corstd50','coravg30','corstd30','coravgUp','corstdUp',...
            'coravgSize','corstdSize','corpeakDuration','corstd_peakDuration');
        

    end % End of Try, Catch
    clearvars -except delim files method num pathname total count backorder backs backmeth
end



%% Create new excel sheet with compiled data

count2 = 1;
dirname2 = ([pathname 'APD Data' delim]);
addpath(dirname2)
files2 = dir(fullfile(dirname2, '*.mat'));
nfiles2 = cell(length(files2),1);
for i = 1:length(files2)
    nfiles2(i) = {files2(i).name};
end

total = length(nfiles2);

% create empty cells to be filled
filenames = cell(length(total),1);

zers = zeros(total,1);
newVars = struct('Amp',zers, 'Amp_std',zers, 'APD30',zers, 'APD30_std',zers,...
    'APD50',zers, 'APD50_std',zers, 'APD80',zers, 'APD80_std',zers, 'APD90',zers,...
    'APD90_std',zers,'UPD',zers, 'UPD_std',zers, 'TotalWidth',zers, 'TotalWidth_std',zers,...
    'BPM',zers, 'delF_f0',zers,'f0',zers,'DVDT',zers, 'DVDT_std',zers,'PeakTime',zers,...
    'PeakTime_std',zers, 'CORUPD',zers, 'CORUPD_std',zers,...
    'CORAPD30',zers, 'CORAPD30_std',zers, 'CORAPD50',zers, 'CORAPD50_std',zers,...
    'CORAPD80',zers, 'CORAPD80_std',zers, 'CORAPD90',zers, 'CORAPD90_std',zers,...
    'CORTotalWidth',zers, 'CORTotalWidth_std',zers,'CORPeakTime',zers,'CORPeakTime_std',zers);


% iterate through the number of files and pull out
while count2 <= total
    
    oldVars1 = struct('avgMax',NaN,'stdMax',NaN,'avg30',NaN,'std30',NaN,...
        'avg50',NaN,'std50',NaN,'avg80',NaN,'std80',NaN,'avg90',NaN,'std90',...
        NaN,'avgUp',NaN,'stdUp',NaN,'avgSize',NaN,'stdSize',NaN,'Beatrate',NaN,...
        'DelF_F0',NaN,'F0',NaN,'dV_dt',NaN,'std_dV_dt',NaN,'peakDuration',NaN,...
        'std_peakDuration',NaN,...
        'coravgUp',NaN,'corstdUp',NaN,'coravg30',NaN,'corstd30',NaN,'coravg50',...
        NaN,'corstd50',NaN,'coravg80',NaN,'corstd80',NaN,'coravg90',NaN,'corstd90',...
        NaN,'coravgSize',NaN,'corstdSize',NaN,'corpeakDuration',NaN,'corstd_peakDuration',NaN);

    oldVarsNames1 = fieldnames(oldVars1);
    filename2 = nfiles2{count2};
    
    
    oldVars2 = load([dirname2 filename2]);
    oldVarsNames2 = fieldnames(oldVars2);
    oldVarsCell2 = oldVarsNames2;
    oldVarsCell2(:,2) = struct2cell(oldVars2);
    
    
    
    missingVars = setdiff(oldVarsNames1,oldVarsNames2);
    if ~isempty(missingVars)
        missingVars(:,2) = num2cell(NaN(size(missingVars,1),1),size(missingVars,1));
        oldVars = cell2struct([oldVarsCell2(:,2);missingVars(:,2)],[oldVarsCell2(:,1);missingVars(:,1)]);
    else
        oldVars = oldVars2;
    end
    
    
    filename2 = filename2(1:end-4);
    f = filename2(1:end-9); % Remove '_APD_Data' from each file name so that you're just getting the video name.
   
    filenames{count2} = f;
    newVars.Amp(count2) = oldVars.avgMax;
    newVars.Amp_std(count2) = oldVars.stdMax;
    newVars.APD30(count2) = oldVars.avg30;
    newVars.APD30_std(count2) = oldVars.std30;
    newVars.APD50(count2) = oldVars.avg50;
    newVars.APD50_std(count2) = oldVars.std50;
    newVars.APD80(count2) = oldVars.avg80;
    newVars.APD80_std(count2) = oldVars.std80;
    newVars.APD90(count2) = oldVars.avg90;
    newVars.APD90_std(count2) = oldVars.std90;
    newVars.UPD(count2) = oldVars.avgUp;
    newVars.UPD_std(count2) = oldVars.stdUp;
    newVars.BPM(count2) = oldVars.Beatrate;
    newVars.TotalWidth(count2) = oldVars.avgSize;
    newVars.TotalWidth_std(count2) = oldVars.stdSize;
    newVars.delF_f0(count2) = oldVars.DelF_F0;
    newVars.f0(count2) = oldVars.F0;
    newVars.DVDT(count2) = oldVars.dV_dt;
    newVars.DVDT_std(count2) = oldVars.std_dV_dt;
    newVars.PeakTime(count2) = oldVars.peakDuration;
    newVars.PeakTime_std(count2) = oldVars.std_peakDuration;
    
    
    newVars.CORAPD30(count2) = oldVars.coravg30;
    newVars.CORAPD30_std(count2) = oldVars.corstd30;
    newVars.CORAPD50(count2) = oldVars.coravg50;
    newVars.CORAPD50_std(count2) = oldVars.corstd50;
    newVars.CORAPD80(count2) = oldVars.coravg80;
    newVars.CORAPD80_std(count2) = oldVars.corstd80;
    newVars.CORAPD90(count2) = oldVars.coravg90;
    newVars.CORAPD90_std(count2) = oldVars.corstd90;
    newVars.CORUPD(count2) = oldVars.coravgUp;
    newVars.CORUPD_std(count2) = oldVars.corstdUp;
    newVars.CORTotalWidth(count2) = oldVars.coravgSize;
    newVars.CORTotalWidth_std(count2) = oldVars.corstdSize;
    newVars.CORPeakTime(count2) = oldVars.corpeakDuration;
    newVars.CORPeakTime_std(count2) = oldVars.corstd_peakDuration;
    
    
    
    count2 = count2 + 1;
end



% Create and Save Table into Folder where Data is Taken From
T = struct2table(newVars,'RowNames', filenames);


writetable(T, [dirname2 'Compiled APD Data.xlsx'],'WriteRowName',true);




end