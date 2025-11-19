
% segments scalp, skull, CSF, grey matter, and white matter from an MRI structural scan
% Uses Freesurfer brain surfaces and freesurfer utilities
% Also expects that iso2mesh is in your file path

% author: Katherine Perdue, PhD

% Please cite the following work if using this code:
% Perdue KL, Diamond SG; 
% T1 magnetic resonance imaging head segmentation for diffuse optical tomography and electroencephalography. 
% J. Biomed. Opt. 19(2):026011.  doi:10.1117/1.JBO.19.2.026011.

%% Load Freesurfer volumes and create a rough brain mask

t1=MRIread([mr_dir filesep 'T1.mgz']);

brain_t1=MRIread([mr_dir filesep 'brainmask.mgz']);
wm=MRIread([mr_dir filesep 'wm.seg.mgz']);

% create a rough mask of the brain
aseg=MRIread([mr_dir filesep 'aseg.mgz']);
bmask=aseg.vol>0.5;
bmask=i2m_close(bmask, 6); 
bmask_d=thickenbinvol(bmask, 1);

%% Segment outer skull

% surface parameters
optskull.radbound=4;
optskull.distbound=1;
dx=0.5; % size of voxels for meshing purposes
mindist=1;
smooth_param=mindist/dx;
volinds=-1:254;

t1_mask=t1.vol> .15*range(t1.vol(:)); 
[bin, regionnum]=bwlabeln(t1_mask, 6);

%%

[histn, ~]=hist(bin(bin>0), 1:regionnum);
[~, k]=max(histn); % find max nonzero connected element

scalp=bin==k;

scalp2=zeros(size(scalp));
scalp3=zeros(size(scalp));
scalp4=zeros(size(scalp));
for ind=1:size(scalp, 2)
    fprintf(['calculating outer skull in slice %f \n' ind])
    [x, y]=find(squeeze(scalp(:, :, ind))'>0);
    if ~isempty(x) && size(x, 1)>2
    [V, S]=alphavol([x y], 100);
    loops=extractloops(S.bnd);
    largest_loop=min(find(isnan(loops)));
    scalp2(:, :, ind)=poly2mask(x(loops([1:largest_loop-1 1]), :)', y(loops([1:largest_loop-1 1]), :)', 256, 256);
    else
        scalp2(:, :, ind)=zeros(256, 256);
    end
    
 
    [x, y]=find(squeeze(scalp(:, ind, :))'>0);
    if ~isempty(x) && size(x, 1)>2
    [V, S]=alphavol([x y], 100);
    loops=extractloops(S.bnd);
    largest_loop=min(find(isnan(loops)));
    scalp3(:, ind, :)=poly2mask(x(loops([1:largest_loop-1 1]), :)', y(loops([1:largest_loop-1 1]), :)', 256, 256);
    else
        scalp3(:, ind, :)=zeros(256, 256);
    end
    
    [x, y]=find(squeeze(scalp(ind, :, :))'>0);
    if ~isempty(x) && size(x, 1)>2
    [V, S]=alphavol([x y], 100);
    loops=extractloops(S.bnd);
    largest_loop=min(find(isnan(loops)));
    scalp4(ind, :, :)=poly2mask(x(loops([1:largest_loop-1 1]), :)', y(loops([1:largest_loop-1 1]), :)', 256, 256);
    else
        scalp4(ind, :, :)=zeros(256, 256);
    end
    
end

scalp=(scalp2+scalp3+scalp4)==3;
scalp=thickenbinvol(scalp, 1);
scalp_e=thinbinvol(scalp, 2); %changed this to 2 to better reflect anatomy

headonly=t1.vol.*(~bmask).*scalp; % also get rid of stuff outside scalp
% ventricles
vents=[4 43 5 14 15 24 44 72];
ventvol=ismember(aseg.vol, vents);
sort_csf=sort(t1.vol(ventvol));
sort_t1=sort(t1.vol(scalp>0));
skull_thresh2=sort_t1(round(.25*length(sort_t1))) ;
skull_thresh=50;

fprintf('calculating outer skull boundary \n')
os=(((headonly<skull_thresh)+(bmask>0)).*scalp_e)>0;
[bin, regionnum]=bwlabeln(os, 6);
[histn, ~]=hist(bin(bin>0), 1:regionnum);
[~, k]=max(histn); % find max nonzero connected element
os=bin==k;

% smooth surface

skull_fill_smooth=smoothbinvol(os, smooth_param);
[skullns, skullfs]=v2s(skull_fill_smooth, 0.75, optskull); 
os=surf2vol(skullns, skullfs, volinds, volinds, volinds);
os=fillholes3d(os, 0);
os=os.*scalp_e;

%% Segment inner skull
fprintf('calculating inner skull boundary \n')
optskull.radbound=2; % allow more wiggles in inner skull
optskull.distbound=1;
max_skull=15;
os_e=thinbinvol(os, 1);
headonly=os_e.*brain_t1.vol;
is=((headonly>skull_thresh2)+(bmask_d>0))>0;
skull_max=thinbinvol(os, max_skull);
is=(is+skull_max)>0;
[bin, regionnum]=bwlabeln(is, 6);
[histn, histx]=hist(bin(bin>0), 1:regionnum);
[j, k]=max(histn); % find max nonzero connected element
is=bin==k;
skull_fill_smooth=smoothbinvol(is, smooth_param);
[skullns, skullfs]=v2s(skull_fill_smooth, 0.48, optskull); 
is=surf2vol(skullns, skullfs, volinds, volinds, volinds);
is=fillholes3d(is, 0);
is=i2m_close(is, 4); 

%%
% cerebellum white matter
cb_white=[7 46];
cbwhitevol=ismember(aseg.vol, cb_white);

% subcortical gray matter;
% 10 and 49 are thalamus proper 
sc_gray=[10 11 12 13 17 18 49 50 51 52 53 54];
scgrayvol=ismember(aseg.vol, sc_gray);

fprintf('writing out segmentation \n')
m_seg=t1;
if include_fiducial
    m_seg.vol=is+os+scalp+(aseg.vol>0.5)+(wm.vol>0.5)+vitE1+vitE2;
else
    m_seg.vol=is+os+scalp+(aseg.vol>0.5)+(wm.vol>0.5);
end
m_seg.vol(ventvol>0)=3; % set ventricles to CSF val
m_seg.vol(cbwhitevol>0)=5; % set cerebellum white to white val
m_seg.vol(scgrayvol>0)=4; % subcortical gray matter to gray val

permuted_vol=permute(m_seg.vol,[2 1 3]);
writevol=permuted_vol(:);
%%

fwrite(fileID,writevol,'uint8');
fclose(fileID)
