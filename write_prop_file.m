function write_prop_file(mua,mus,g,n,filepath,session_id)

% this function writes the optical property .dat file
% input:
% mua: absorption coefficient for each tissue layer, mm-1
% mus: scattering coefficient for each tissue layer, mm-1 (not reduced)
% g: anisotropy for each tissue layer
% n: refractive index for each tissue layer
% session_id: name of session id

% output:
% none, file gets written out

% author: Melissa Wu (wu.melissa.m <at> gmail.com)

% this file is part of scatterBrains

file_content=cell(length(mua)+1,1);

file_content{1}=['1 ' num2str(length(mua))];

for I=1:length(mua)
    file_content{I+1}=[num2str(I) ' ' num2str(mua(I)) ' ' num2str(mus(I)) ' ' num2str(g(I)) ' ' num2str(n(I))];
end

fid=fopen([filepath filesep 'prop_' session_id '.dat'],'w');

for line=1:length(file_content)
    fprintf(fid,[file_content{line} '\n']);
end

fclose(fid);