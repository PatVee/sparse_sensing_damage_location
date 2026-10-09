function kw51_writeText(path,text)
    f=fopen(path,'w'); assert(f>=0); cleanup=onCleanup(@() fclose(f)); %#ok<NASGU>
    fprintf(f,'%s\n',text);
end
