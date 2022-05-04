ascc_list = [0.0002, 0.0004, 0.0006, 0.0008, 0.001, 0.002, 0.004, 0.006, 0.008, 0.01, 0.02, 0.04, 0.06, 0.08, 0.1, 0.2];
scc_set = 0.0001;

nTrees = 20;

scc_1 = zeros(1, length(ascc_list));
ascc_1 = zeros(1, length(ascc_list));

performance_1 = zeros(1, length(ascc_list));
counter = 1;

for i = 1:length(ascc_list)

    if scc_set < ascc_list(i)

        data = load(['/home/ys589/project/assemble/decoding_performance/decoding_performance' '_' num2str(scc_set) '_' num2str(ascc_list(i)) '_' num2str(nTrees) '.mat']);

        scc_1(counter) = scc_set;
        ascc_1(counter) = ascc_list(i);
        performance_1(counter) = data.ave_performance;

        counter = counter + 1;
    end
end

scc_list = [0.0001, 0.0002, 0.0004, 0.0006, 0.0008, 0.001, 0.002, 0.004, 0.006, 0.008, 0.01, 0.02, 0.04, 0.06, 0.08, 0.1];

ascc_set = 0.2;

scc_2 = zeros(1, length(scc_list));
ascc_2 = zeros(1, length(scc_list));

performance_2 = zeros(1, length(scc_list));
counter = 1;

for i = 1:length(scc_list)

    if scc_list(i) < ascc_set 

        data = load(['/home/ys589/project/assemble/decoding_performance/decoding_performance' '_' num2str(scc_list(i)) '_' num2str(ascc_set) '_' num2str(nTrees) '.mat']);

        scc_2(counter) = scc_list(i);
        ascc_2(counter) = ascc_set;
        performance_2(counter) = data.ave_performance;

        counter = counter + 1;
    end
end

%tbl_1 = table(transpose(scc_1), transpose(ascc_1), transpose(performance_1), 'VariableNames', {'scc','ascc','performance'});
%tbl_1

lm_ascc = fitlm(log10(transpose(ascc_1)),performance_1);
lm_ascc

%tbl_2 = table(transpose(scc_2), transpose(ascc_2), transpose(performance_2), 'VariableNames', {'scc','ascc','performance'});
%tbl_2

lm_scc = fitlm(log10(transpose(scc_2)),performance_2);
lm_scc

figure;
plot(lm_ascc)
title('decoding performance vs ascc (scc = ' num2str(scc_set) ')')
xlabel('ascc (log scale)')
ylabel('decoding performance')

figure;
plot(lm_scc)
title('decoding performance vs scc (ascc = ' num2str(ascc_set) ')')
xlabel('ascc (log scale)')
ylabel('decoding performance')

figure;

%coeffs = lm_scc.Coefficients.Estimate;

%scatter(ax1, ascc, performance);
%hold on
%h1 = refline(coeffs(2), coeffs(1));
%h1 = lsline(ax1);
%h1.Color = 'r';

fit_1 = polyfit(log10(ascc_1), performance_1, 1);
fitline_1 = polyval(fit_1, log10(ascc_1));

subplot(1,2,1);
plot(log10(ascc_1), performance_1, '.')
xlim([-5 0])
xticks([-4, -3, -2, -1])
xticklabels({'10^{-4}','10^{-3}','10^{-2}', '10^{-1}', '10^{0}'})
hold on
plot(log10(ascc_1), fitline_1, '-r')
hold off
title('decoding performance vs ascc (scc = ' num2str(scc_set) ')')
xlabel('ascc (log scale)')
ylabel('decoding performance')

fit_2 = polyfit(log10(scc_2), performance_2, 1);
fitline_2 = polyval(fit_2, log10(scc_2));

subplot(1,2,2);
plot(log10(scc_2), performance_2, '.')
xlim([-5 0])
xticks([-4, -3, -2, -1])
xticklabels({'10^{-4}','10^{-3}','10^{-2}', '10^{-1}', '10^{0}'})
hold on
plot(log10(scc_2), fitline_2, '-r')
hold off
title('decoding performance vs scc (ascc =' num2str(ascc_set) ')');
xlabel('scc (log scale)');
ylabel('decoding performance');