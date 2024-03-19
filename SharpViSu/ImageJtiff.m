function [] = ImageJtiff(MultiDimImg, FileName, pixsize)
%saves a multi-channel tiff image for ImageJ
%MultiDimImg = zeros(300,400,4,5,6,'uint16'); - X, Y, Color Channels, Z-Slices,
%T-Frames;
%example for a two-color image: MultiDimImg = cat(3, R, G);
%https://www.mathworks.com/matlabcentral/answers/389765-how-can-i-save-an-image-with-four-channels-or-more-into-an-imagej-compatible-tiff-format
if ~exist("pixsize", "var")
    pixsize = 1;
end
fiji_descr = convertCharsToStrings(['ImageJ=1.52p' newline ...
            'images=' num2str(size(MultiDimImg,3)*...
                              size(MultiDimImg,4)*...
                              size(MultiDimImg,5)) newline... 
            'channels=' num2str(size(MultiDimImg,3)) newline...
            'slices=' num2str(size(MultiDimImg,4)) newline...
            'frames=' num2str(size(MultiDimImg,5)) newline... 
            'hyperstack=true' newline...
            'mode=composite' newline...  
            'loop=false' newline...  
            ['min=' num2str(min(MultiDimImg,[],'all'),'%.1f')] newline...      
            ['max=' num2str(max(MultiDimImg,[],'all'),'%.1f')] newline...
            'unit=nm']);            

if isa(MultiDimImg,'double')
bits = 64;
elseif isa(MultiDimImg,'single')
bits = 32;
elseif isa(MultiDimImg,'uint16')
bits = 16;
elseif isa(MultiDimImg,'uint8')
bits = 8;
end

t = Tiff(FileName,'w'); % file name 
tagstruct.ImageLength = size(MultiDimImg,1);
tagstruct.ImageWidth = size(MultiDimImg,2);
tagstruct.Photometric = Tiff.Photometric.MinIsBlack;
tagstruct.BitsPerSample = bits;
tagstruct.SamplesPerPixel = 1;
tagstruct.Compression = Tiff.Compression.LZW;
tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
if bits < 20
tagstruct.SampleFormat = Tiff.SampleFormat.UInt;
else
tagstruct.SampleFormat = Tiff.SampleFormat.IEEEFP;
end
tagstruct.ImageDescription = fiji_descr;
tagstruct.XResolution = 1/pixsize;
tagstruct.YResolution = 1/pixsize;

for frame = 1:size(MultiDimImg,5)
    for slice = 1:size(MultiDimImg,4)
        for channel = 1:size(MultiDimImg,3)
            t.setTag(tagstruct)
            %t.write(im2uint16(MultiDimImg(:,:,channel,slice,frame)));
            t.write(MultiDimImg(:,:,channel,slice,frame));
            t.writeDirectory(); % saves a new page in the tiff file
        end
    end
end
t.close() 
end