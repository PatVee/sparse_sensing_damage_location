function labels=kw51_normaliseLabels(value)
    if ischar(value), labels=string(cellstr(value)); else, labels=string(value); end
    labels=strtrim(labels(:));
end
