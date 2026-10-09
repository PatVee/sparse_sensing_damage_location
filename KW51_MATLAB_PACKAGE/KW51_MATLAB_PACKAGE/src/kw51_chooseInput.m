function path=kw51_chooseInput(pattern,prompt)
    [name,folder]=uigetfile(pattern,prompt);
    if isequal(name,0), error('KW51:Cancelled','Input selection was cancelled.'); end
    path=fullfile(folder,name);
end
