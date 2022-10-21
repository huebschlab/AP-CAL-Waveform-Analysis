function CMOSconverter11_Auto(fdim1,fdim2,sigma,ngauss,scale,upperframelimit,inputfreq,drawbounds)

%If code throws an error when correct input variables are used, please run
%the "setJavaPath.m" code using the instructions in that script

% To pull framerate from metadata, set inputfreq to 'auto'. To manually
% input framerate, set inputfreq to integer value of framerate


%CMOSconverter10_3_Auto(3,3,0.5,2,1,5e3,'auto',0);


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
    [~,~,filetype] = fileparts([pathname filename]);


    nfilename = [filename(1:end-4) '.mat'];
    tic


    %% Read in File and MetaData

    if matches(filetype,'.tif')
        tifinfo = imfinfo([pathname filename]);
        PlaneCount = length(tifinfo);


        if PlaneCount > 1

            if (isa(inputfreq, 'char') || isa(inputfreq,'string')) && matches(inputfreq,'auto') % If frequnecy is set to auto, pull FPS from metadata
                try
                    imageDescript = tifinfo.ImageDescription;
                    splitDescript = splitlines(imageDescript);
                    splitDescript(~contains(splitDescript,'fps')) = [];
                    fpsStr = cell2mat(splitDescript);
                    fpsQuoteStart = find(fpsStr == '=');
                    frequency = str2double(fpsStr(fpsQuoteStart+1:end));
                    acqFreq = 1/frequency;
                catch % If FPS cannot be found in metadata, query user
                    try
                        frequency = str2double(cell2mat(inputdlg('Frame rate could not be found. Please input Frame rate (Frames/Second).','No Frame rate')));
                        acqFreq = 1/frequency;
                    catch
                        disp('Frame rate must be integer')
                    end
                end


            else % If Frequency is not set to auto
                try
                    frequency = inputfreq;
                    acqFreq = 1/frequency;
                catch
                    disp('Inputfreq must be integer or auto')
                end
            end
        end

    else %file type is not a .tif

        reader = bfGetReader([pathname filename]);
        omeMeta = reader.getMetadataStore();
        PlaneCount = omeMeta.getPlaneCount(0);



        if PlaneCount > 1
            % Get Framerate


            if (isa(inputfreq, 'char') || isa(inputfreq,'string')) && matches(inputfreq,'auto')
                try
                    omeXML = char(omeMeta.dumpXML());
                    omeXML = omeXML(1:1500);

                    k=strfind(omeXML, "ExposureTime=");
                    if ~isempty(k)
                        k = k(1);
                        xpotime = omeXML(k+14:end);
                        xpoquoteend = find(xpotime == '"',1,'first');
                        xpotime = str2double(xpotime(1:xpoquoteend-1));
                    end


                    k1 = strfind(omeXML, 'Plane DeltaT=');
                    if ~isempty(k1)
                        k1 = k1(2);
                        ftime = omeXML(k1+14:end);
                        fquoteend = find(ftime == '"',1,'first');
                        ftime = str2double(ftime(1:fquoteend-1));
                    end

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
                catch
                    try
                        frequency = str2double(cell2mat(inputdlg('Frame rate could not be found. Please input Frame rate (Frames/Second).','No Frame rate')));
                        acqFreq = 1/frequency;
                    catch
                        disp('Frame rate must be integer')
                    end
                end
            else
                try
                    frequency = inputfreq;
                    acqFreq = 1/frequency;
                catch
                    disp('Inputfreq must be integer or auto')
                end
            end

        else %PlateCount < 1
            sprintf('Snapshot, No Video Detected  For')
            disp(char([pathname filename]));

        end % Planecount ("is it a snapshot or a video?")
    end

    %% Read-In the spatially-averaged signal as the average over ALL pixels at each time-slice



    maxcount = min(PlaneCount,upperframelimit);


    if matches(filetype,'.tif') % Read in Tif
        testslice = int32(imread([pathname filename],1));
        testslice = multiRGB(testslice);
    else % Read in microscopy file types (.nd2, .cxd, etc.)
        testslice = int32((bfGetPlane(reader,1)));
    end



    % Filter Test Slice
    testslice2 = testslice;
    h = fspecial('gaussian',fdim2,sigma);
    testslice2 = imresize(testslice2,scale);
    [m,n] = size(testslice2);
    cmosData = int32(zeros(m,n,maxcount));


    %% Find ROIs
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
            mask = int32(ones(m,n));


                for jj = 1:maxcount
                    if matches(filetype,'.tif')
                        imageslice = int32(imread([pathname filename],jj));
                        imageslice = multiRGB(imageslice);
                    else
                        imageslice = int32((bfGetPlane(reader,jj)));
                    end
                    imageslice = imresize(imageslice,scale);
                    imageslice2 = medfilt2(imageslice,[fdim1 fdim1]);
                    for ii = 1:ngauss
                        imageslice2 = imfilter(imageslice2,h);
                    end

                    imageslice2 = mask.*imresize(imageslice2,scale);
                    cmosData(1:m,1:n,jj) = imageslice2;
                end

            

            bgimage = cmosData(:,:,1);

            %NOTE: Matlab fails to properly save files > 2GB, which would include zero-padded time series
            nfilename2 = [pathname nfilename];


            save(nfilename2,'cmosData','acqFreq','frequency', 'bgimage','dual');


        else
            mask = zeros(m,n,numROI);
            h2 = gobjects(numROI,1);

            for i = 1:numROI
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


            %%
            for jj = 1:maxcount
                if matches(filetype,'.tif')
                    imageslice = int32(imread([pathname filename],jj));
                    imageslice = multiRGB(imageslice);
                else
                    imageslice = int32((bfGetPlane(reader,jj)));
                end
                imageslice = imresize(imageslice,scale);
                imageslice2 = medfilt2(imageslice,[fdim1 fdim1]);
                for ii = 1:ngauss
                    imageslice2 = imfilter(imageslice2,h);
                end

                cmosData(1:m,1:n,jj) = imageslice2;


            end


            %%




            cmosData1 = cmosData;
            for i = 1:numROI
                cmosData = mask(:,:,i).*cmosData1;

                bgimage = cmosData(:,:,1);

                %NOTE: Matlab fails to properly save files > 2GB, which would include zero-padded time series
                if i == 1
                    nfilename2 = [pathname nfilename];
                end

                name = strcat(nfilename2(1:end-5),num2str(i),'.mat');
                save(name,'cmosData','acqFreq','frequency', 'bgimage','dual','-v7.3');
            end
        end


    else
        mask = int32(ones(m,n));


        for jj = 1:maxcount
            if matches(filetype,'.tif')
                imageslice = int32(imread([pathname filename],jj));
                imageslice = multiRGB(imageslice);
            else
                imageslice = int32((bfGetPlane(reader,jj)));
            end

                imageslice = imresize(imageslice,scale);
                imageslice2 = medfilt2(imageslice,[fdim1 fdim1]);
                for ii = 1:ngauss
                    imageslice2 = imfilter(imageslice2,h);
                end
                imageslice2 = mask.*imresize(imageslice2,scale);
                cmosData(1:m,1:n,jj) = imageslice2;
        end

        bgimage = cmosData(:,:,1);

        %NOTE: Matlab fails to properly save files > 2GB, which would include zero-padded time series
        nfilename2 = [pathname nfilename];


        save(nfilename2,'cmosData','acqFreq','frequency', 'bgimage','dual');

    end





    toc

    fprintf('Remaining videos to convert: %d\n', total - count);

    count = count+1;
end
end


%% if tif data type is not a matrix (RGB, transparancy, etc.)
function data = multiRGB(data)
if ~ismatrix(data)
    if size(data,3)>3
        data = data(:,:,3);
    end
    if length(find(all(data,[1 2]))) ~= 1
        data = rgb2gray(data);
    else
        data(:,:,all(~data,[1 2])) = [];
    end
end


end
