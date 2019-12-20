function out = igseq_runMultiDB( directory )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

%get list of db files from directory
list = dir([directory '/' '*.db']);

%output variable that tempOut gets vertcat into
out = {};
%tempory variable for unique isotype usage
totalIso = {};

%loop through list of DB files
  for j = 1:length(list)
     db = [list(j).folder '/' list(j).name];
     pop = igseq_upop(db);
     iso = repmat({''}, length(pop), 1);
     col = mat2cell(parula(length(pop)), ones(1, length(pop)),3);
     tempOut = {};
     tempOut(:,2) = pop;
     tempOut(:,1) = {list(j).name};
     
     %vj plot
     igseq_vjplot(db,'sequence')
     print(strcat(db,'_vjplot_sequence.eps'),'-depsc');
     
     %clonality plots
%      for i = 1:size(pop)
%          n = pop{i};
%          igseq_clonalityplot(db,n,n,'b',0,[])
%          print(strcat(n,'_clonality.eps'),'-depsc');
%      end
     
     %get number of lineages
     tmp = igseq_nlin(db);
     tempOut(:,3) = tmp(:,2);
     
     %get number of sequences
     tmp = igseq_nseq(db);
     tempOut(:,4) = tmp(:,2);
     
     %get number of unique sequences
     tmp = igseq_nUniqueSeq(db);
     tempOut(:,5) = tmp(:,2);
     
     %get D values for 20 and 50
     tmp = igseq_dx(db,[20 50]);
     t = tmp(strcmp(tmp(:,2),''),:);
     s = size(t);
     tempOut(:,6:(6 + (s(2)-4))) = t(:,4:s(2));
     
     %calculate Hill Diversity
     tmp = igseq_hillplot(db,0:1:2,0,pop,iso,col);
     r = num2cell(tmp);
     t = r';
     s = size(t);
     tempOut(:,8:(8 + (s(2)-2))) = t(:,2:s(2));
     
     %get total mutation average
     tmp = igseq_mutavg(db,'total');
     tempOut(:,10) = tmp(:,2);
     
     %get silent mutation average
     tmp = igseq_mutavg(db,'silent');
     tempOut(:,11) = tmp(:,2);
     
     %get non-silent mutation average
     tmp = igseq_mutavg(db,'nonsilent');
     tempOut(:,12) = tmp(:,2);
     
     %get isotypes
     tmp = igseq_isotypes(db);
     totalIso = [totalIso tmp(1,2:end)];
     
     out = vertcat(out, tempOut);
     
  end
  
  totalIso = unique(totalIso);
  outIso = {};
  outPercentIso = {};
  
  %aggregate Isotype information across databases
  for j = 1:length(list)
    tmpOut = {};
    tmpPercent = {};
    db = [list(j).folder '/' list(j).name];
    tmp = igseq_isotypes(db);
    s = igseq_nseq(db);
    
    for i = 1:length(totalIso)
        if ~isempty(find(strcmp(tmp(1,:),totalIso(i))))
            values = tmp(:,strcmp(tmp(1,:),totalIso(i)));
            tmpOut(:,i) = values(2:end);
            tmpPercent(:,i)  = num2cell(cell2mat(tmpOut(:,i)) ./ cell2mat(s(:,2)) * 100);
            
        else
            tmpOut(:,i) = num2cell(zeros(length(tmp)-1,1));
            tmpPercent(:,i)  = num2cell(zeros(length(tmp)-1,1));
        end
    end
    
    outIso = vertcat(outIso, tmpOut);
    outPercentIso = vertcat(outPercentIso, tmpPercent);
    
  end
  
  % add header and print out file
  out = [out outIso outPercentIso];
  header = ['DB' 'Pop' 'Number of Lineages' 'Number of Sequences' 'Number Unique Sequences' 'D20' 'D50' 'q1' 'q2' 'Avg total mutation rate' 'Avg silent mutation rate' 'Avg nonsilent mutation rate' totalIso strcat(totalIso, ' Percentage')];
  out = vertcat(header,out);
  cell2csv('db_stats.csv',out,',','2016','.');
end
