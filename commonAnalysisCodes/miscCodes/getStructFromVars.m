function rawAnalysis = getStructFromVars(varargin)
    for iArg = 1:nargin
        varName = inputname(iArg);
        rawAnalysis.(varName) = varargin{iArg};
    end
end