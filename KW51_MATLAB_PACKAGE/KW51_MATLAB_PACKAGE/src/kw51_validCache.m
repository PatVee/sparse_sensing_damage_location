function yes=kw51_validCache(path,signature)
    yes=false;
    if isfile(path)
        info=whos('-file',path);
        if ismember('runSignature',{info.name})
            saved=load(path,'runSignature'); yes=isequaln(saved.runSignature,signature);
        end
    end
end

%% ------------------------------------------------------------------------
