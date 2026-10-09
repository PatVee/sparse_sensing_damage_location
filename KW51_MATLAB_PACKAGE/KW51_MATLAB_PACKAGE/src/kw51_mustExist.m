function kw51_mustExist(path,label)
    if ~isfile(path), error('KW51:MissingInput','Missing %s: %s',label,path); end
end
