function [] = savePKLtoMAT(fn, scc, ascc)

plot_flag = 1;

MaxFeat = 30;

% read the numpy MD array in pickle files and save to a MD matrix
% I used Table data container, you can use something else too.
% Main idea is to save it to a MAT file so that we can use
% powerful matrix manipulation for further analysis.

for i = 1:MaxFeat
 
    fid = py.open([fn num2str(scc) '_' num2str(ascc) '_' num2str(i) '.pickle'],'rb');
    VT = py.pickle.load(fid);
    VTraces(i, :, :, :) = VT.double;
end

%Save Traces to MAT file
save(['/Users/soga1/Downloads/assemble/post_syn_mat/PostSyn_' num2str(scc) '_' num2str(ascc) '.mat'], "VTraces")