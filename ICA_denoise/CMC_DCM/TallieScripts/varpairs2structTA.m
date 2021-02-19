function [V,flds] = varpairs2structTA(varargin)


if ~isempty(varargin)
    names = varargin(1:2:length(varargin));
    values = varargin(2:2:length(varargin));
    for k = 1:length(names)
        V.(names{k}) = values{k};
    end
else V = struct();
end
flds = fieldnames(V);

end

