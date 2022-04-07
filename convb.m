function c=convb(g,r)
% This function convolutes the image g with an out of focus blurring kernel 
% with reflection boundary. function c=convb(g,r);
%r is the radius of the point spread function.

PSF = fspecial('disk', r);
% to do: try other PSF's, Gaussian, motion blur etc.

p = (size(PSF, 1) - 1) / 2;
g = padarray(g,[p p], 'symmetric');
c = conv2(g, PSF, 'valid');