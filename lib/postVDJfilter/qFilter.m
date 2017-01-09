function [ikeep] = qFilter( QualityScoreAlign, name, q1start, q1stop,q1, q2start,q2stop,q2 )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    q1start = str2num(q1start);
    q1stop = str2num(q1stop);
    q1 = str2num(q1);
    q2start = str2num(q2start);
    %q2stop = str2num(q2stop);
    q2 = str2num(q2);
    
    
    n = size(QualityScoreAlign,1);
    n2 = size(QualityScoreAlign,2);
    
    r1 = zeros(n,1);
    r2 = zeros(n,1);
    
    m1 = nanmean(QualityScoreAlign(:,q1start:q1stop),2);
    m2 = nanmean(QualityScoreAlign(:,q2start:q2stop),2);
    
    r1(m1 >= q1) = 1;
    r2(m2 >= q2) = 1;
    
    disp('Number of sequences passing threshold 1')
    disp(length(find(m1 >= q1)))
    
    disp('Number of sequences passing threshold 2')
    disp(length(find(m2 >= q2)))
    
    ikeep=find(prod([r1 r2],2)>0);
    
    disp('Number of overlapping sequences')
    disp(length(ikeep))
    
     baseviz_r2(QualityScoreAlign(ikeep,:));
     print(strcat(name,'qualityScoreDistribution_filtered.eps'),'-depsc');
     
     baseqc_r2(QualityScoreAlign(ikeep,:));
     print(strcat(name,'qualityHistogram_filtered.eps'),'-depsc');

end

