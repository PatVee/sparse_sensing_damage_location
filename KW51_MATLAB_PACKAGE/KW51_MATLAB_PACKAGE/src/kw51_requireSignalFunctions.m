function kw51_requireSignalFunctions()
    for name={'resample','butter','filtfilt','dct'}
        if exist(name{1},'file')==0, error('KW51:Toolbox','Signal Processing Toolbox required: %s',name{1}); end
    end
end
