function dump_p0_3_subject_vectors(root, outCsv)
%DUMP_P0_3_SUBJECT_VECTORS Minimal data bridge for the independent verifier.
% Dumps ONLY the persisted subject_metrics tables of all 640 Mpower
% checkpoints to a single TSV. Contains zero inference logic; the
% independent exact sign-flip recompute happens in Python
% (verify_p0_3_ablation_independent.py). Mphase subject vectors are plain
% float64 arrays in the checkpoint and are read by the verifier directly
% via h5py, so they are not dumped here.
arguments
    root (1,1) string
    outCsv (1,1) string
end
paper=jsondecode(fileread(fullfile(root,"scripts","ieeg","paper","acc00_sim_p0_3_ablation_contract.json")));
base=fullfile(root,string(paper.outputs.root));
rows=cell(0,8);
for i=1:numel(paper.cells.Mpower)
    cellId=string(paper.cells.Mpower(i).id);
    for s=1:double(paper.cells.Mpower(i).worlds)
        p=fullfile(base,"mpower","checkpoints",cellId,sprintf("world_%03d.mat",s));
        q=load(p,"checkpoint");c=q.checkpoint;
        assert(string(c.cell_id)==cellId&&c.seed_index==s,"identity mismatch at %s",p);
        t=c.subject_metrics;
        for j=1:height(t)
            rows(end+1,:)={char(cellId),s,j,char(string(t.subject(j))),char(string(t.run(j))),t.Mpower(j),t.Mphase(j),t.n_pairs(j)}; %#ok<AGROW>
        end
    end
end
T=cell2table(rows,"VariableNames",["cell_id","seed_index","row_index","subject","run","Mpower","Mphase","n_pairs"]);
writetable(T,outCsv,"FileType","text","Delimiter","\t");
fprintf("DUMPED %d rows\n",height(T));
end
