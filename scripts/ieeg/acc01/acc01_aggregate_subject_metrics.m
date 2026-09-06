function out = acc01_aggregate_subject_metrics(rows)
subject=string({rows.subject})';run=string({rows.run})';keys=unique(table(subject,run),'rows','sorted');
x=nan(height(keys),2);nPairs=zeros(height(keys),1);
for i=1:height(keys)
    keep=subject==keys.subject(i)&run==keys.run(i);r=rows(keep);
    x(i,:)=[median([r.power],'omitnan'),median([r.phase],'omitnan')];nPairs(i)=nnz(keep);
end
t=table(keys.subject,keys.run,x(:,1),x(:,2),nPairs,'VariableNames',["subject" "run" "Mpower" "Mphase" "n_pairs"]);
out=struct('subjects',keys.subject,'runs',keys.run,'matrix',x,'table',t);
end
