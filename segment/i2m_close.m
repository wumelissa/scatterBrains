function [ img ] = i2m_close( img, arg )
% Close using iso2mesh functions, fill holes as middle step

 img=thickenbinvol(img, arg);
 img=fillholes3d(img, 0);
 img=thinbinvol(img, arg);

end

