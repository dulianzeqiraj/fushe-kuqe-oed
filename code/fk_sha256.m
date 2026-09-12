function h = fk_sha256(path, nchar)
%FK_SHA256  Hex SHA-256 of a file, truncated to nchar characters (default 16).
%
%  Through the JVM that ships with MATLAB, so no toolbox and no shell call.

if nargin < 2, nchar = 16; end
fid = fopen(path, 'r');
if fid < 0, error('fk_sha256:open', 'cannot read %s', path); end
bytes = fread(fid, inf, '*uint8');
fclose(fid);

md = java.security.MessageDigest.getInstance('SHA-256');
md.update(typecast(bytes, 'int8'));
d = typecast(md.digest(), 'uint8');
h = lower(reshape(dec2hex(d, 2)', 1, []));
h = h(1:min(nchar, numel(h)));
end
