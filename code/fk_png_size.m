function [w, h] = fk_png_size(path)
%FK_PNG_SIZE  Pixel width and height from a PNG header.
%
%  The IHDR chunk is always the first chunk and always at byte 17, so the
%  dimensions are two big-endian 32-bit integers at offsets 17 and 21.  Read
%  directly rather than with imread, which would decode the whole image.

fid = fopen(path, 'r');
if fid < 0, error('fk_png_size:open', 'cannot read %s', path); end
c = onCleanup(@() fclose(fid)); %#ok<NASGU>

sig = fread(fid, 8, '*uint8')';
if ~isequal(sig, uint8([137 80 78 71 13 10 26 10]))
    error('fk_png_size:notpng', '%s is not a PNG', path);
end
fseek(fid, 16, 'bof');
w = double(fread(fid, 1, 'uint32', 0, 'b'));
h = double(fread(fid, 1, 'uint32', 0, 'b'));
if w <= 0 || h <= 0
    error('fk_png_size:bad', 'implausible size in %s', path);
end
end
