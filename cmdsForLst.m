function cmdsForLst(lst,varargin)
disp(nargin);
for i=1:nargin-1
    disp(varargin{i});
    runLst(lst,varargin{i});
end
end