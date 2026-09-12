function fk_zip_replace(template, outfile, partname, content)
%FK_ZIP_REPLACE  Copy a zip package entry by entry, swapping one member.
%
%  Used to write a .docx by replacing word/document.xml in a template while
%  every other part, styles, numbering, theme, fonts, is carried across
%  untouched.  That is the same discipline this project uses when editing a
%  Word file rather than rebuilding it.
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
%  zero bytes.  The document then looks structurally sound and Word refuses
%  it as corrupt.  The parts are therefore extracted to a temporary folder and
%  read with fread, which is the direction that works.

outfile = local_abspath(outfile);
template = local_abspath(template);

tmp = tempname;
mkdir(tmp);
c1 = onCleanup(@() rmdir(tmp, 's')); %#ok<NASGU>
unzip(template, tmp);

names = local_entry_names(template);
newbytes = unicode2native(content, 'UTF-8');

if exist(outfile, 'file'), delete(outfile); end
fos = java.io.FileOutputStream(outfile);
zos = java.util.zip.ZipOutputStream(fos);
c2 = onCleanup(@() local_close(zos, fos)); %#ok<NASGU>

seen = false;
for k = 1:numel(names)
    name = names{k};
    zos.putNextEntry(java.util.zip.ZipEntry(name));
    if strcmp(name, partname)
        bytes = newbytes;
        seen = true;
    else
        bytes = local_readbytes(fullfile(tmp, strrep(name, '/', filesep)));
    end
    if ~isempty(bytes)
        zos.write(typecast(bytes(:)', 'int8'), 0, numel(bytes));
    end
    zos.closeEntry();
end
if ~seen
    error('fk_zip_replace:missing', '%s is not in %s', partname, template);
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
    error('fk_zip_replace:read', 'cannot read extracted part %s', path);
end
b = fread(fid, inf, '*uint8');
fclose(fid);
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
