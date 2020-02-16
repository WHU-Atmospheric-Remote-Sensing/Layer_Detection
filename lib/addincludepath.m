% add all the subdirectories in the folder of '../include/'
% History:
%   2019-08-14 Add the comments by Zhenping Yin

includePath = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'include');
addpath(includePath);

%% find subdirectories in lib 
subdirs = listdir(includePath);

for iSubdir = 1:length(subdirs)
    addpath(subdirs{iSubdir});
end