function value=kw51_sha256File(path)
    md=java.security.MessageDigest.getInstance('SHA-256');
    f=fopen(path,'rb'); assert(f>=0,'Cannot open file: %s',path); cleanup=onCleanup(@() fclose(f)); %#ok<NASGU>
    while true
        bytes=fread(f,1024*1024,'*uint8'); if isempty(bytes), break; end
        md.update(typecast(bytes(:),'int8'));
    end
    bytes=typecast(md.digest(),'uint8'); value=lower(reshape(dec2hex(bytes,2).',1,[]));
end
