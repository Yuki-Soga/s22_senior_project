ascc_list = [0.2];
scc_list = [0.12, 0.14, 0.16, 0.18];

for i = 1:length(scc_list)

    for j = 1:length(ascc_list)

        if scc_list(i) < ascc_list(j)

            scc = scc_list(i);
            ascc = ascc_list(j);
   
            fn = ['/Users/soga1/Downloads/assemble/post_syn_pickle/' num2str(scc) '_' num2str(ascc) '/PostSyn_'];
            
            savePKLtoMAT(fn, scc, ascc)
            
            load(['/Users/soga1/Downloads/assemble/post_syn_mat/PostSyn_' num2str(scc) '_' num2str(ascc) '.mat'])
            
            dim = size(VTraces);
            
            SOMA = 2; % site ID
            dt = 0.1; % sim time step. check with NEURON params
            THRESHOLD = 0;
            CorrWin = 400; % ms
            stimOnset = 100; % ms
            
            % Init
            PostSynSpikeCount = zeros(dim(1), dim(3));
            PostSynSpikeTimes = cell(dim(1), dim(3));
            
            % Calc
            for fID = 1:dim(1) % features
            
                for tID = 1:dim(3) % trials per feature
                    tmp(1,:) = VTraces(fID,SOMA,tID,:);
                    [a,b] = count_spike_times(tmp,THRESHOLD);
                    PostSynSpikeCount(fID,tID) = sum(b*dt>stimOnset & b*dt<=CorrWin);
                    PostSynSpikeTimes{fID,tID}(:) = b*dt';
                end
            end
            
            % Save
            save(['/Users/soga1/Downloads/assemble/post_syn_spikes/PostSynSpikes_' num2str(scc) '_' num2str(ascc) '.mat'], 'PostSynSpikeCount', 'PostSynSpikeTimes');
            
            % Num of trials to plot
            trial_num = 100;
            
            % Test resp var
            clear tmp
            % clr = {'r','g','b','c','m','y','k'};
            % for i = 1:7, figure; tmp(1,:) = PostCellSpikeCount(i,1:trial_num); plot(tmp,'o','color', clr{i}); xlim([0 trial_num]); ylim([0 10]); %hold on;
            % end
        end
    end
end