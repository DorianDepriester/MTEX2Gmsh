function helptoc2yml(xmlFile, outFile)
%HELPTOC2YML Convert helptoc.xml into the Jekyll data file _data/toc.yml.
%   HELPTOC2YML() reads docs/html/helptoc.xml and writes docs/_data/toc.yml so
%   that the sidebar of the GitHub Pages site stays in sync with the table of
%   contents displayed by the MATLAB help browser. Run it after each change
%   of helptoc.xml.
%
%   HELPTOC2YML(XMLFILE, OUTFILE) uses the given input and output files.

if nargin < 1
    here = fileparts(mfilename('fullpath'));
    xmlFile = fullfile(here, 'html', 'helptoc.xml');
end
if nargin < 2
    here = fileparts(mfilename('fullpath'));
    outFile = fullfile(here, '_data', 'toc.yml');
end

doc = xmlread(xmlFile);
root = doc.getDocumentElement();
items = directChildren(root, 'tocitem');

fid = fopen(outFile, 'w', 'n', 'UTF-8');
if fid < 0
    error('helptoc2yml:cannotOpen', 'Cannot open %s for writing.', outFile);
end
cleanup = onCleanup(@() fclose(fid)); %#ok<NASGU>
writeItems(fid, items, 0);
end

function children = directChildren(node, tag)
children = {};
cl = node.getChildNodes();
for k = 0:cl.getLength-1
    c = cl.item(k);
    if c.getNodeType() == c.ELEMENT_NODE && strcmp(char(c.getNodeName()), tag)
        children{end+1} = c; %#ok<AGROW>
    end
end
end

function writeItems(fid, items, indent)
pad = repmat(' ', 1, indent);
for k = 1:numel(items)
    node = items{k};
    target = char(node.getAttribute('target'));
    txt = '';
    cl = node.getChildNodes();
    for j = 0:cl.getLength-1
        c = cl.item(j);
        if c.getNodeType() == c.TEXT_NODE || c.getNodeType() == c.CDATA_SECTION_NODE
            txt = [txt char(c.getNodeValue())]; %#ok<AGROW>
        end
    end
    label = strtrim(txt);
    fprintf(fid, '%s- title: %s\n', pad, yq(label));
    fprintf(fid, '%s  target: %s\n', pad, yq(target));
    children = directChildren(node, 'tocitem');
    if ~isempty(children)
        fprintf(fid, '%s  children:\n', pad);
        writeItems(fid, children, indent + 4);
    end
end
end

function s = yq(s)
s = ['"' strrep(strrep(s, '\', '\\'), '"', '\"') '"'];
end
