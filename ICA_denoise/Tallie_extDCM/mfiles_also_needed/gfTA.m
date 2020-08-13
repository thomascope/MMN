function fn = gfTA(varargin)

fn = dir;
fn = arrayfun(@(x) x.name,fn,'Uni',0);
fn(1:2) = [];
fn = sort_nat(fn);
if nargin==1
    idx = cellfun(@(x) regexp(x,varargin{1}),fn,'Uni',0);
    logi = cellfun(@isempty,idx);
    fn(logi) = [];
end

end

