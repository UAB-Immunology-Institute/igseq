function igseq_generateFASTA(db,l,info)

% Chris Fucile
% Alex Rosenberg
% August 7 2017
% University of Alabama at Birmingham, Informatics Institute
% Copyright (c) 2017. All rights reserved.
% This software is offered with no guarantees of any kind.
db_file = strsplit(db, '/');
db_file = strsplit(db_file{end}, '.');
if l == 0
    conn = igseq_db(db);
    q = 'SELECT distinct(lineageID) from sequence';
    l = cell2mat(get(fetch(exec(conn, q)), 'data'));
    close(conn);
end
for i = 1:length(l)
 % open database and get populations, lineages IDs and isotypes
    if info == 1
       igseq_generateINFO(db,l,0); 
    end
    conn = igseq_db(db);
    q = ['SELECT sequenceID, population, Vgene, Dgene, Jgene, NT_VDJ from sequence where lineageID = ' num2str(l(i)) ''];
    a = get(fetch(exec(conn, q)), 'data');
    close(conn);
   
    fid = fopen([db_file{1} '_lineage_' num2str(l(i)) '.fasta'],'w');
    for j = 1:size(a)
        seq = regexprep(a(j,6),'\.','');
        
        fprintf(fid,'%1s\n',['>' a{j,1} '|' a{j,2} '|' a{j,3} '|' a{j,4} '|' a{j,5}]);
        fprintf(fid,'%1s\n',seq{1});
    end
    fclose(fid);
end