function fk_zip_write(template, outfile, names, contents)
%FK_ZIP_WRITE  Copy a zip package, replacing some members and adding others.
%
%  fk_zip_write(template, outfile, names, contents)
%
%  names    cellstr of part names, e.g. {'word/document.xml', 'word/media/image1.png'}
%  contents cell of matching payloads: char for XML, uint8 for binary
%
%  A name already in the template replaces that part; a name not in it is
%  appended.  Everything else is carried across untouched, which is how this
%  project edits a Word file rather than rebuilding it.
%
%  Two details, both learned the hard way.
%
%  Writing goes through the JVM that ships with MATLAB rather than through
%  zip(), because the package contains [Content_Types].xml and the brackets
%  are wildcard characters to MATLAB's file matching.
%
%  Reading does NOT go through java.io.InputStream.read(byte[]).  MATLAB
%  passes arrays to Java by value, so the JVM fills its own copy, the caller
%  gets its buffer back unchanged, and every copied part comes out as a run of
%  zero bytes.  The document then looks structurally sound and Word refuses it
%  as corrupt.  The parts are therefore extracted to a temporary folder and
%  read with fread, which is the direction that works.

if ischar(names), names = {names}; end
if ~iscell(contents), contents = {contents}; end
assert(numel(names) == numel(contents), 'names and contents must match');

outfile  = local_abspath(outfile);
template = local_abspath(template);

tmp = tempname;
mkdir(tmp);
c1 = onCleanup(@() rmdir(tmp, 's')); %#ok<NASGU>
unzip(template, tmp);

existing = local_entry_names(template);
payload = containers.Map('KeyType', 'char', 'ValueType', 'any');
for k = 1:numel(names)
    v = contents{k};
    if ischar(v) || isstring(v)
        v = unicode2native(char(v), 'UTF-8');
    end
    payload(names{k}) = uint8(v(:))';
end

extra = setdiff(names, existing, 'stable');
order = [existing(:); extra(:)];

if exist(outfile, 'file'), delete(outfile); end
fos = java.io.FileOutputStream(outfile);
zos = java.util.zip.ZipOutputStream(fos);
c2 = onCleanup(@() local_close(zos, fos)); %#ok<NASGU>

for k = 1:numel(order)
    name = order{k};
    zos.putNextEntry(java.util.zip.ZipEntry(name));
    if isKey(payload, name)
        bytes = payload(name);
    else
        bytes = local_readbytes(fullfile(tmp, strrep(name, '/', filesep)));
    end
    if ~isempty(bytes)
        zos.write(typecast(uint8(bytes(:)'), 'int8'), 0, numel(bytes));
    end
    zos.closeEntry();
end
end

% =========================================================================
function names = local_entry_names(zipfile)
zf = java.util.zip.ZipFile(zipfile);
c = onCleanup(@() zf.close()); %#ok<NASGU>
names = {};
en = zf.entries();
while en.hasMoreElements()
    e = en.nextElement();
    if ~e.isDirectory()
        names{end+1} = char(e.getName()); %#ok<AGROW>
    end
end
end

function b = local_readbytes(path)
fid = fopen(path, 'r');
if fid < 0
    error('fk_zip_write:read', 'cannot read extracted part %s', path);
end
b = fread(fid, inf, '*uint8');
fclose(fid);
b = b(:)';
end

function local_close(zos, fos)
try, zos.close(); catch, end %#ok<CTCH>
try, fos.close(); catch, end %#ok<CTCH>
end

function p = local_abspath(p)
if isempty(regexp(p, '^([A-Za-z]:|\\\\|/)', 'once'))
    p = fullfile(pwd, p);
end
end
