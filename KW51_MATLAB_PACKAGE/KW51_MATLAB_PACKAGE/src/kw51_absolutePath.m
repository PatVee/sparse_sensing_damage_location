function path=kw51_absolutePath(path)
    jf=javaObject('java.io.File',char(path)); path=char(jf.getCanonicalPath());
end
