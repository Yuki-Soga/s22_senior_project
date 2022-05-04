function [] = saveMATtoPKL(fn, scc, ascc)

load(fn)

dim = size(spike_times);

% Stimulus onset latency
stimOnset = 100; % ms

for nrnID = 1:dim(2)
    for featID = 1:dim(1)
        trials = spike_times(featID,nrnID,:);
        trl = py.list();
        for t = 1:length(trials)
            spt = spike_times{featID,nrnID,t} + stimOnset;
            trl.append(py.numpy.asarray(spt, 'int').tolist());
        end

        fid = py.open(['Feat' num2str(featID) '_SynGrp' num2str(nrnID) '_' num2str(scc) '_' num2str(ascc) '.pickle'],'wb');
        py.pickle.dump(trl,fid)
    end
end