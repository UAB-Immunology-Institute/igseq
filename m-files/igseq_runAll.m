function iso = igseq_runAll( db )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
     pop = igseq_upop(db);
     igseq_nlin(db)
     igseq_nseq(db)
     igseq_nUniqueSeq(db)
     iso = igseq_isotypes(db);
     igseq_vjplot(db,'sequence')
     print('vjplot_sequence.eps','-depsc');
     igseq_vjplot(db,'lineage')
     print('vjplot_lineage.eps','-depsc');
     igseq_dx(db,20)
     igseq_dx(db,50)
     for i = 1:size(pop)
         n = pop{i};
         igseq_clonalityplot(db,n,n,'b',0,[])
         print(strcat(n,'_clonality.eps'),'-depsc');
     end
    
     igseq_muthistplot(db,'total',1,40,pop,iso,col);

end

