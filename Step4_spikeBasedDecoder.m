ascc_list = [0.2];
scc_list = [0.0002, 0.0004, 0.0006, 0.0008, 0.001, 0.002, 0.004, 0.006, 0.008, 0.01, 0.02, 0.04, 0.06, 0.08, 0.1, 0.12, 0.14, 0.16, 0.18];

% Number of trees we want in the forest
nTrees = 20;

for i = 1:length(scc_list)

    for j = 1:length(ascc_list)

        if scc_list(i) < ascc_list(j)

            scc = scc_list(i);
            ascc = ascc_list(j);
            
            fn = '/Users/soga1/Downloads/assemble/post_syn_spikes/PostSynSpikes';
            
            % Load post syn spikes data
            data = load([fn '_' num2str(scc) '_' num2str(ascc) '.mat']);
            
            % Extract post syn spike count data
            spike_count = data.PostSynSpikeCount;
            
            % Extract size of the post cell spike count data
            dim = size(spike_count);
            
            rng default
            
            % For how many folds we calculate decoding accuracy (different train-test 
            % split every time but note that this is NOT a k-fold cross validation)
            fold_num = 100;
            
            % How many trials to use for testing for each fold
            test_num = 20;
            
            % Store decoding performance of all folds
            performance = zeros(1, fold_num);
            
            % For each fold, the assignment of columns to training set and testing set 
            % will be decided randomly
            for k = 1:fold_num
            
                % Choose which columns to use for testing
                test_index = randperm(dim(2), test_num);
            
                % Concatenate them to create a testing matrix (30 x test_num)
                test_spike_count = spike_count(:, test_index);
            
                % Reshape the testing matrix into an array
                test_feature = transpose(reshape(transpose(test_spike_count), [1, dim(1)*test_num]));
                
                train_spike_count = spike_count;
            
                % Remove the columns used for testing to create a training matrix
                train_spike_count(:, test_index) = [];
                
                % Reshape the training matrix into an array (30 x (number of trials - test_num))
                train_feature = transpose(reshape(transpose(train_spike_count), [1, dim(1)*(dim(2)-test_num)]));
                
                % Training and testing labels
                test_label = [];
                train_label = [];
            
                for l = 1:dim(1)
            
                    test_label = [test_label, l.* ones(1, test_num)];
                    train_label = [train_label, l.* ones(1, dim(2)-test_num)];
                end
            
                % Train the TreeBagger 
                classifier = TreeBagger(nTrees, train_feature, transpose(train_label), 'Method', 'classification');
            
                % Use the trained decision forest to make predictions
                prediction = predict(classifier, test_feature);
                prediction = str2double(transpose(prediction));
            
                % Calculate decoding performance
                performance(k) = 1 - mean(test_label~=prediction);
            
            end
            
            ave_performance = mean(performance);
            
            display(ave_performance)
            
            % Save mean decoding performance as well as decoding performance for all trials 
            save(['/Users/soga1/Downloads/assemble/decoding_performance/decoding_performance' '_' num2str(scc) '_' num2str(ascc) '_' num2str(nTrees) '.mat'], 'ave_performance', 'performance')
        end
    end
end