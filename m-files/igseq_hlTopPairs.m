function n = igseq_hlTopPairs(db)

% IGSEQ_HLTopPairs lineages that having matching heavy and light chains
%
%    t = igseq_hlTopPairs(dbname, idx) returns a cell array that
%    summarizes 
%
%    Requires Matlab Database Toolbox
%
%    See also IGSEQ_MHOVERLAP, IGSEQ_GROUPEDLIN, IGSEQ_COMMONLIN

% Chris Fucile
% Alex Rosenberg
% 9 Feb 2018
% University of Alabama at Birmingham, Informatics Institute
% Copyright (c) 2018. All rights reserved.
% This software is offered with no guarantees of any kind.


    conn = igseq_db(db);
    q = 'select pairedChain, count(distinct(chainType)) as n from sequence group by pairedChain having n = 2';
    a = get(fetch(exec(conn, q)), 'data');
    
    q2 = char(strcat('select lineageID, VGene, pairedChain from sequence where pairedChain != ''null'''));
    a2 = get(fetch(exec(conn, q2)), 'data');
    
    nseq = size(a,1);
    n = {};
    nindex = 1;
    for i = 1:nseq
        %q2 = char(strcat('select lineageID, VGene, count(*) from sequence where pairedChain = ''', a(i,1), ''' group by lineageID'));
        %a2 = get(fetch(exec(conn, q2)), 'data');
        tmpa2 = a2(find(strcmp(a2(:,3),a(i,1))),:);
        t  = tabulate(cell2mat(tmpa2(:,1)));
        t2 = tabulate(tmpa2(:,2));
        t = t(find(t(:,2) ~=0),1:2);
        
        for j = 1:size(t,1)
            n{nindex,1} = a(i,1);
            n(nindex,2:3) = tmpa2(j,:);
            n(nindex,4) = num2cell(t(j,2));
            
            idx = find([tmpa2{:}] == t(j,1));
            
            n(nindex,5) = tmpa2(idx(1),2);
            nindex = nindex + 1;
        end
    end
    
   return