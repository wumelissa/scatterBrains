
function data=getPhotonHistory(history_filename)

% loads .mch file, assumes that momentum transfer information is available

% input: .mch filename

% output: data array with shape [ndets 2*ntissue_layers+1]
% columns for data array are [det_idx pl_tissue1 pl_tissue2 ... pl_tissuen
% mt_tissue1 mt_tissue2 ...mt_tissuen]

% author: Melissa Wu (wu.melissa.m <at> gmail.com)
% contributing author: Stefan Carp (stefan.carp@mgh.harvard.edu)

% this file is part of scatterBrains
% License: GPLv3

[mch_data,mch_header]=loadmch(history_filename);

if (mch_header.recordnum-2)<(2*mch_header.medianum),
    error('History file does not contain momentum transfer information \n');
end

pl_index_start=2+mch_header.medianum;
pl_index_end=pl_index_start+mch_header.medianum-1;
mt_index_start=pl_index_end+1;
mt_index_end=mt_index_start+mch_header.medianum-1;
data=mch_data(:,[1 pl_index_start:pl_index_end mt_index_start:mt_index_end]);
