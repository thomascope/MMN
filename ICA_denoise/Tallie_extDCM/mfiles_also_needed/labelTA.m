function str = labelTA(str,varargin)

if nargin==2
    delim = varargin{1};
else delim = '_';
end
if ~iscell(str), str = {str}; end
str = cellfun(@(x) strsplit(x,'_'),str,'Uni',0);
str = cellfun(@(x) strcat(x{:}),str,'Uni',0);

end

