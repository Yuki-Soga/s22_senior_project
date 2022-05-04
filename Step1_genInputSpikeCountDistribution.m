%%
function Step1_genInputSpikeCountDistribution(scc, ascc, MaxTrials)
% Frame work for generating spiking patterns from multiple
% presynaptic neurons to a post synaptic neuron. 
% Let us say, there are three presynaptic neurons with distinct 
% tuning functions for a set of features (MaxFeat).
% This code does the following:
% 1. Generate mean spike counts (SC) for each neuron,   
%    with different levels of correlations between pairs (scc),
%    and within the neuron (ascc), for a number of trials (MaxTrials).
% 2. Simulate spike trains with these mean SCs.
% 3. Estimate SC for all simulated neuronal spike trains
%    and verify levels of correlations are mainatined between them.
%    This is for verification only. Step 2 generates the output for
%    next stage of the project.
% 
% To be done: Save these input patterns and use them to study
% in the post synpatic neuron the consequences of the different 
% SC correlation patterns in the input. 
% 
% Jadi, March 2021, Yale
%----------------------------------------------------------------
plotFig = 1;

if(plotFig)
    figure()
    set(gca, 'fontsize',14)
end

% Color associated with each neuron for plotting
clr = {'r', 'b', 'g','c','m','y','k'};

%% Define Tunning functions
MaxFeat = 30; %number of features
features = [1:MaxFeat];
MaxNrn = 3; %number of neurons in the population
tunDiff = 4; % Overlap: circular tunning distance of neurons in the population.
for nrnID = 1:MaxNrn
    prefFeat = 1+tunDiff*(nrnID-1); % peak of tunning function
    tunFun = 10*(2 + sin(2*pi*(features-prefFeat)/max(features)));
    tunFun(tunFun<0) = 0;
    TF(nrnID,:)= tunFun;
end

if(plotFig)
        subplot(2,3,1:2),
        p1=plot(features, TF(1,:), 'Color', clr{1}); hold on;
        p2=plot(features, TF(2,:), 'Color', clr{2}); hold on;
        p3=plot(features, TF(3,:), 'Color', clr{3}); hold on;
        xlabel('Stimulus \theta');ylabel('Spike Count')
        set(gca, 'fontsize',20)
end

%%Define the spike time data struct
% MaxTrials = 10;
for feat = 1:MaxFeat
    for nrn = 1:MaxNrn
        for trl = 1:MaxTrials
            spike_times{feat,nrn,trl} = [];
        end
    end
end

%%
% 1. Generate mean spike counts (SC) for all neurons with 
% different levels of correlations between them.
% 2. Generate spike trains with these mean SCs.
%------------------------------------------------------
for featID = [1:MaxFeat] % Loop through all features
    
%     scc = 0.015; %cross-corr coeff
%     ascc = .02; %auto-corr coeff
    corrType = 0; %0 = +ve, 1 = -ve
    mu = [TF(1,featID) TF(2,featID) TF(3,featID)];
    
    sigma = [ascc scc/2 scc;
             scc/2 ascc scc/3;
             scc scc/3 ascc];
    
    if(corrType)
        sigma = sigma.*[1 -1 1;
            -1 1 -1;
            1 -1 1];
    end
   
    % 1. Generate mean SC with the above stats over many experimental trials
    trials = round(mvnrnd(mu,sigma,MaxTrials)); %MaxTrialsxMaxNrn
    
    % Visualize the scatter plot of SC over all "recorded" neurons
    if(plotFig)
        % Plot SCs
        %for nrnID = 1:MaxNrn
        %    subplot(2,3,1:2), plot(featID, trials(:,nrnID), '.', 'Color',clr{nrnID}), hold on
        %end

        
        legend([p1(1), p2(1), p3(1)], {'Neuron1','Neuron2','Neuron3'});
        
        subplot(2,3,3), plot3(trials(:,1), trials(:,2),trials(:,3),'.');
        xlabel('Neuron1');ylabel('Neuron2'); zlabel('Neuron3');
        set(gca, 'fontsize',20)
        hold on;
        subplot(2,3,4), plot(trials(:,3), trials(:,2),'.');
        xlabel('N3');ylabel('N2');
        set(gca, 'fontsize',14)
        hold on;
        subplot(2,3,5), plot(trials(:,1), trials(:,2),'.');
        xlabel('N1');ylabel('N2');
        set(gca, 'fontsize',14)
        hold on;
        subplot(2,3,6), plot(trials(:,1), trials(:,3),'.');
        xlabel('N1');ylabel('N3');
        set(gca, 'fontsize',14)
        hold on;
    end
    
    % 2. Simulate spike trains
    corWin = 200; %ms window of spike-count correlations
    for nrnID = 1:MaxNrn
        for trialID = 1:MaxTrials
            st = sort(randi(corWin, [trials(trialID,nrnID),1]));%us
            spike_times{featID, nrnID, trialID} = st(1:end);%ms

            %spk_tim = zeros(1,corWin);
            %spk_tim(st) = 1;
            %spike_times(featID, trialID, nrnID,:) = spk_tim;
        end
    end
end



% figure()
% %%
% % 3. Estimate SC for all simulated neuronal spike trains
% % and verify levels of correlations are mainatined between them.
% %----------------------------------------------------------------
% corWin = 200; %ms window of spike-counts
% for featID = [1:1:MaxFeat]
%     
%     % Estimate SC
%     for nrnID = 1:MaxNrn
%         for trialID = 1:MaxTrials
%             st = spike_times{featID, nrnID, trialID};
%             trials(trialID, nrnID) = length(st);
%         end
%     end
%     
%     % Visualize the scatter plot of *estimated* SC over all "recorded" neurons
%     if(plotFig)
%         % Plot SCs
%         subplot(2,3,3), plot3(trials(:,1), trials(:,2),trials(:,3),'.');
%         xlabel('N1');ylabel('N2'); zlabel('N3');
%         set(gca, 'fontsize',14);xlim([0 4]); ylim([0 4]); zlim([0 4]); axis tight
%         hold on;
%         subplot(2,3,4), plot(trials(:,3), trials(:,2),'.'); axis tight
%         xlabel('N3');ylabel('N2');
%         set(gca, 'fontsize',14);xlim([0 4]); ylim([0 4]); axis tight
%         hold on;
%         subplot(2,3,5), plot(trials(:,1), trials(:,2),'.'); axis tight
%         xlabel('N1');ylabel('N2');
%         set(gca, 'fontsize',14);xlim([0 4]); ylim([0 4]); axis tight
%         hold on;
%         subplot(2,3,6), plot(trials(:,1), trials(:,3),'.'); axis tight
%         xlabel('N1');ylabel('N3');
%         set(gca, 'fontsize',14);xlim([0 4]); ylim([0 4]); axis tight
%         hold on;
%     end
% end

% Save the simulated spike trains for the specified spike count
% distribution
fn = 'presynSpikingData';
spikeTimeDSformat = 'spike_times{feature,neuron,trial}';
save([fn '.mat'], 'spike_times', 'scc', 'ascc', 'spikeTimeDSformat')
saveMATtoPKL(fn)
