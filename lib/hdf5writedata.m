function [] = hdf5writedata(filename, location, data, dataAttr, deflate, chunksize)
%hdf5writedata Write data and data attributes to HDF5 file.
%   Example:
%       [] = hdf5writedata(filename, location, data, dataAttr, deflate, chunksize)
%   Inputs:
%       filename: char
%           HDF5 filename. (with absolute path)
%       location: char
%           location of the dataset.
%       data: array
%           exported data.
%       dataAttr: struct
%           dataset attributes.
%       deflate: integer
%           compression level.
%       chunksize: 2-element array
%           chunking size.
%   Outputs:
%       
%   History:
%       2019-11-11. First Edition by Zhenping
%   Contact:
%       zp.yin@whu.edu.cn

if ~ exist('chunksize', 'var')
    chunksize = ceil(size(data)/8);
end

if ~ exist('deflate', 'var')
    deflate = 6;
end

if iscell(data) || ischar(data)
    slashPos = strfind(location, '/');

    if isempty(slashPos)
        error('Not a valid position: %s', location);
    end

    dset_details.Location = location(1:slashPos(end));
    dset_details.Name = location(slashPos(end) + 1:end);

    hdf5write(filename, dset_details, data, 'WriteMode', 'append');

    keys = fieldnames(dataAttr);
    for iKey = 1:length(keys)
        attr_details.Name = keys{iKey};
        attr_details.AttachedTo = location;
        attr_details.AttachType = 'dataset';

        hdf5write(filename, attr_details, dataAttr.(keys{iKey}), 'WriteMode', 'append');
    end

elseif isnumeric(data)

    h5create(filename, location, size(data), 'DataType', class(data), 'Deflate', deflate, 'Chunksize', chunksize);
    h5write(filename, location, data);
    
    keys = fieldnames(dataAttr);
    for iKey = 1:length(keys)
        attName = keys{iKey};
        h5writeatt(filename, location, attName, dataAttr.(attName));
    end

else
    error('Unsupported dataset for HDF5 file.');
end

end