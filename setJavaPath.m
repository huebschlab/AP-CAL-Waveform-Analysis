% 1.Before running go to https://www.openmicroscopy.org/bio-formats/downloads/
% and download the Bio-Formats Package. 
% 
% 2. Place the downloaded file, titled "bioformats_package.jar"
% in the bfmatlab folder to replace the old version of the file and reset
% the matlab path.
% 
%3.Run this function and select the new "bioformats_package.jar" file
%
%4. This should solve the issue. If it does not, download the the matlab
%toolbox from the above link and replace the existing bfmatlab folder with
%the new downloaded version, reseting the matlab path, and then run this code, selecting the
%"bioformats_package.jar" in the new bfmatlab folder.
%
%5. If this still does not work, try to restart matlab prior to downloading
%the new files and running this code.


function setJavaPath


[file,pathname] = uigetfile('*.*', 'Select One or More Files' , 'MultiSelect', 'off');
javaaddpath([pathname file]);
end
