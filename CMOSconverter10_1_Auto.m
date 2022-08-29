function cmosData = CMOSconverter10_1_Auto(fdim1,fdim2,sigma,ngauss,scale,upperframelimit,inputfreq,drawbounds)

%If code throws an error when correct input variables are used, please run
%the "setJavaPath.m" code using the instructions in that script

% To pull framerate from metadata, set inputfreq to 'auto'. To manually
% input framerate, set inputfreq to integer value of framerate


%cmosData = CMOSconverter10_1_Auto(3,3,0.5,2,1,5e3,'auto',0);


warning off
%parameters that need to be encoded
% Modified 3/20/2020 by Nate Huebsch
%
% If the "drawbounds" is set to 1, the user can input to only analyze a
% selected ROI. The program includes a message box that allows a delay, so
% that the user can zoom in and pan to select a very precise ROI.

% Notes - 8/13/18 Nate Huebsch: this was created to allow importing files
% from microscopy formats not supported by the original rhythm toolbox, like the proprietary
% files from Zeiss and Nikon systems
% INPUTS
% Filetype =  microscopy filename extension (e.g. 'czi', 'cxd', etc)
% NOTE USED --> olddir = directory where file is located (use uigetdir instead)
% NOT USED --> oldfilename = filename (use uigetdir in combination with
% filetype search instead)
%
%fdim1 and fdim2 are dimensions of the median and gaussian filters,
%respectively, used to pre-smooth data; sigma is the standard deviation
%applied in Gaussian filtering; ngauss is the number of (iterative)
%gaussian filtration steps applied
%scale = size to scale up or down (e.g. 0.25 does a 4x4 bin)
%upperlimit = maximim time step allowed (do this to optimize file size)cl
%
%inputfreq = input frequency (add in case bioformats incorrectly reads in
%acquisition frequency) Set 'auto' to automatically pull from 
% OUTPUT
% cmosData = pertinent data as saved in this structure including:
%                + intensity = cmosData and cmosData2c
%                + analog channels = channel1 and channel2
%                + time per frame = acqFrequency (msec)
%                + rate of acquisition = frequency (Hz)
%                + background image = bgimage
%                + dual camera setting = dual(1 - yes / 0 - no)

dual = 0;

[files,pathname] = uigetfile('*.*', 'Select One or More Files' , 'MultiSelect', 'on');

if isa(files, 'char') == 1
    total = 1;
else
    total = max(size(files));
end

count = 1;


while count <= total
    if isa(files, 'char') == 1
        filename = files;
    else
        filename = files{1,count};
    end
    disp(filename);


    nfilename = [filename(1:end-4) '.mat'];
    tic

    % Read-In File and MetaData

    warning off

    %     reader = loci.formats.Memoizer(bfGetReader(), 0);
    % reader.setId([pathname filename]);
    % reader.close();
    % loci.common.DebugTools.enableLogging('INFO');

    reader = bfGetReader([pathname filename]);
    omeMeta = reader.getMetadataStore();
    PlaneCount = omeMeta.getPlaneCount(0);



    if PlaneCount > 1
        % Get Framerate


        if (isa(inputfreq, 'char') || isa(inputfreq,'string')) && matches(inputfreq,'auto')
            omeXML = char(omeMeta.dumpXML());
            omeXML = omeXML(1:1500);
            k=strfind(omeXML, "ExposureTime=");
            k = k(1);
            k1 = strfind(omeXML, 'Plane DeltaT=');
            k1 = k1(2);

            xpotime = omeXML(k+14:end);
            ftime = omeXML(k1+14:end);
            xpoquoteend = find(xpotime == '"',1,'first');
            fquoteend = find(ftime == '"',1,'first');
            xpotime = str2double(xpotime(1:xpoquoteend-1));
            ftime = str2double(ftime(1:fquoteend-1));

            if ftime ~= 0 && ~isnan(ftime) && contains(omeXML, "http://www.openmicroscopy.org/Schemas/OME/2016-06")
                frequency = 1/ftime; %Aquisition Frequency
                acqFreq = ftime; % Acquisition time in msec
            elseif xpotime ~= 0 && ~isnan(xpotime) &&contains(omeXML, "http://www.openmicroscopy.org/Schemas/OME/2015-01")
                frequency = 1000/xpotime; %Aquisition Frequency
                acqFreq = xpotime/1000;
                disp('Updating Bioformats Recommended. Aquired Framerate may be incorrect.')
            elseif xpotime == 0 || isnan(xpotime)
                fprintf(['Unable to Find Exposure Time for '   filename]);
                error('Please check bioformats OME schema version or input framerate')
            end

        else
            try
                frequency = inputfreq;
            catch
                disp('Inputfreq must be integer or auto')
            end
        end
        %else
        
        % end


        % Read-In the spatially-averaged signal as the average over ALL pixels
        % at each time-slice

        %     bgimage(1:m,1:n) = bgimage1;


        %limitsize = 1e4;
        % upperframelimit = 1e4;
        %maxcount = min(PlaneCount,limitsize);

        maxcount = min(PlaneCount,upperframelimit);

        testslice = int32((bfGetPlane(reader,1)));
        %   imageslice = imrotate(imageslice,90,'bilinear','crop');
        %cmosData(:,:,jj) = imageslice;





        testslice2 = testslice; %medfilt2(imageslice,[fdim1 fdim1]);

        h = fspecial('gaussian',fdim2,sigma);
        %     imageslice2 = imfilter(imageslice2,h);
        %     %imageslice2 = imgaussfilt(imageslice2,fdim2); %Median filter with a fdim x fdim box; gaussian filter with a
        %     %standard deviation of 2
        testslice2 = imresize(testslice2,scale);
        [m,n] = size(testslice2);
        cmosData = int32(zeros(m,n,maxcount));

        if drawbounds

            ChooseFigure = figure;
            set(gcf,'numbertitle','off','name','Choose Area of Interest');
            imagesc(testslice2), colormap(gray)
            axis image
            set(gca, 'YTick', []);
            set(gca, 'XTick', []);

            inputcheck = 0;
            while inputcheck ~= 1
                try
                    numROI = str2double(cell2mat(inputdlg('Number of ROI to be drawn?','ROI Number',1,{'1'})));
                    if numROI >= 0
                        inputcheck = 1;
                    else
                    end
                catch
                end

            end

            if numROI == 0
                close(ChooseFigure);
                %can we initialize bioformats to parfor this loop????
                for jj = 1:maxcount
                    %     reader = javaObject('loci.formats.Memoizer', bfGetReader(), 0);
                    %     reader.setId(char(filenames(ii)));
                    imageslice = int32((bfGetPlane(reader,jj))); %squeeze(tiffStackOriginal(:,:,i));
                    %   imageslice = imrotate(imageslice,90,'bilinear','crop');
                    %cmosData(:,:,jj) = imageslice;
                    imageslice = imresize(imageslice,scale);
                    imageslice2 = medfilt2(imageslice,[fdim1 fdim1]);

                    for ii = 1:ngauss
                        imageslice2 = imfilter(imageslice2,h);
                    end
                    %imageslice2 = imgaussfilt(imageslice2,fdim2); %Median filter with a fdim x fdim box; gaussian filter with a
                    %standard deviation of 2

                    %cmosData(1:m,1:n,jj) = imageslice2;
                    cmosData(1:m,1:n,jj) = imageslice2;


                end



                bgimage = cmosData(:,:,1);

                %NOTE: Matlab fails to properly save files > 2GB, which would include zero-padded time series
                nfilename2 = [pathname nfilename];
                %nfilename2 = [pathname 'test.mat'];
                %dummyvar = 0;

                save(nfilename2,'cmosData','acqFreq','frequency', 'bgimage','dual');
                %save(nfilename2,'dummyvar','cmosData');



            else
                mask = zeros(m,n,numROI);
                h2 = gobjects(numROI,1);

                %  What is this? ask Nate
                for i = 1:numROI
                    %d = msgbox('Use zoom and pan to find desired ROI. Click "OK" when done');
                    uiwait(msgbox('Use zoom and pan to find desired ROI. Click "OK" when done'));

                    h2(i) = drawfreehand; % user outlines ROI on the image
                    mask(:,:,i) = createMask(h2(i));
                    % drawnow;
                end
                close(ChooseFigure);
                mask = int32(mask);
                % Create New Filename for each new cropped Code - Change max of
                % i if you plan to do more than 10 different regions
                for i = 1:numROI
                    nfilename =  [filename(1:end-4) '_ROI_' num2str(i) '.mat'];
                    if isequal(exist([pathname nfilename], 'file'),0)
                        break
                    end
                end
                %can we initialize bioformats to parfor this loop????
                for jj = 1:maxcount
                    %     reader = javaObject('loci.formats.Memoizer', bfGetReader(), 0);
                    %     reader.setId(char(filenames(ii)));
                    imageslice = int32((bfGetPlane(reader,jj))); %squeeze(tiffStackOriginal(:,:,i));
                    %   imageslice = imrotate(imageslice,90,'bilinear','crop');
                    %cmosData(:,:,jj) = imageslice;
                    imageslice = imresize(imageslice,scale);
                    imageslice2 = medfilt2(imageslice,[fdim1 fdim1]);
                    for ii = 1:ngauss
                        imageslice2 = imfilter(imageslice2,h);
                    end
                    %imageslice2 = imgaussfilt(imageslice2,fdim2); %Median filter with a fdim x fdim box; gaussian filter with a
                    %standard deviation of 2

                    %imageslice2 = mask(:,:,i).*imresize(imageslice2,scale);
                    %cmosData(1:m,1:n,jj) = imageslice2;
                    cmosData(1:m,1:n,jj) = imageslice2;


                end

                cmosData1 = cmosData;
                for i = 1:numROI
                    cmosData = mask(:,:,i).*cmosData1;

                    bgimage = cmosData(:,:,1);

                    %NOTE: Matlab fails to properly save files > 2GB, which would include zero-padded time series
                    if i == 1
                        nfilename2 = [pathname nfilename];
                    end
                    %nfilename2 = [pathname 'test.mat'];
                    %dummyvar = 0;



                    name = strcat(nfilename2(1:end-5),num2str(i),'.mat');
                    save(name,'cmosData','acqFreq','frequency', 'bgimage','dual','-v7.3');
                    %save(nfilename2,'dummyvar','cmosData');

                    %         if PlateCount < 1
                    %                 sprintf('Snapshot, No Video Detected  For')
                    %                 disp(char([pathname filename]));
                end
            end
        else
            mask = int32(ones(m,n));
            %can we initialize bioformats to parfor this loop????
            for jj = 1:maxcount
                %     reader = javaObject('loci.formats.Memoizer', bfGetReader(), 0);
                %     reader.setId(char(filenames(ii)));
                imageslice = int32((bfGetPlane(reader,jj))); %squeeze(tiffStackOriginal(:,:,i));
                %   imageslice = imrotate(imageslice,90,'bilinear','crop');
                %cmosData(:,:,jj) = imageslice;
                imageslice = imresize(imageslice,scale);
                imageslice2 = medfilt2(imageslice,[fdim1 fdim1]);
                for ii = 1:ngauss
                    imageslice2 = imfilter(imageslice2,h);
                end
                %imageslice2 = imgaussfilt(imageslice2,fdim2); %Median filter with a fdim x fdim box; gaussian filter with a
                %standard deviation of 2

                imageslice2 = mask.*imresize(imageslice2,scale);
                %cmosData(1:m,1:n,jj) = imageslice2;
                cmosData(1:m,1:n,jj) = imageslice2;


            end

            bgimage = cmosData(:,:,1);

            %NOTE: Matlab fails to properly save files > 2GB, which would include zero-padded time series
            nfilename2 = [pathname nfilename];
            %nfilename2 = [pathname 'test.mat'];
            %dummyvar = 0;

            save(nfilename2,'cmosData','acqFreq','frequency', 'bgimage','dual');
            %save(nfilename2,'dummyvar','cmosData');

        end



    else %PlateCount < 1
        sprintf('Snapshot, No Video Detected  For')
        disp(char([pathname filename]));

    end % Planecount ("is it a snapshot or a video?")

    toc

    fprintf('Remaining videos to convert: %d\n', total - count);

    count = count+1;
end
