%% plot Fig 2C + D
% freq_all = AllRegions_AllCellTypeFreq.csv
% freq_infect = AllRegions_InfectedCellTypeFreq.csv
% freq_ifn = AllRegions_IFNCellTypeFreq.csv

% find mean of cell type frequency per region 
mean_all = mean(freq_all, 2);
mean_infect = mean(freq_infect, 2);
mean_ifn = mean(freq_ifn, 2);
std_infect = std(freq_infect, [], 2); 
std_ifn = std(freq_ifn, [], 2); 
std_all = std(freq_all, [], 2); 

custom_hex = ["#1F77B4","#AEC7E8","#FF7F0E","#FFBB78","#2CA02C", "#98DF8A","#D62728",...
    "#FF9896","#9467BD","#8C564B","#C49C94","#E377C2","#F7B6D2","#7F7F7F",...
    "#C7C7C7", "#BCBD22","#DBDB8D","#17BECF","#9EDAE5","#393B79","#637939",... 
    "#B5CF6B","#843C39"];
figure
errorbar(mean_all(:),mean_infect(:),std_infect/sqrt(10),'k')
hold on
for k = 1:length(mean_infect)
    plot(mean_all(k,:),  mean_infect(k,:),'.', 'Color',custom_hex(k), 'MarkerSize', 45)
    hold on
end


figure
errorbar(mean_infect(:),mean_ifn(:),std_ifn/sqrt(10),'k')
hold on
for k = 1:length(mean_infect)
    plot(mean_infect(k,:),  mean_ifn(k,:),'.', 'Color',custom_hex(k), 'MarkerSize', 45)
    hold on
end

%% Plot Fig 2F
% comp_foci = AllFoci_CellTypeFreq.csv
% frac = AllFoci_FracIFN.csv

% separate airway and alveolar foci based on the frequency of airway epi 
% first two columns are airway
alv_ind = find(comp_foci(:,1) + comp_foci(:,2) == 0);
epi_ind = find(comp_foci(:,1) + comp_foci(:,2) ~=0); 

%plot all data using scatter plot
alv_foci = frac(alv_ind,:);
figure; scatter([1:23], frac(alv_ind,:), 'ko')
hold on
scatter([1:23], frac(epi_ind,:), 'kx')

% plot error bars
mean_frac = nanmean(frac);
stderr = nanstd(frac)/sqrt(85);
hold on; errorbar( mean_frac, stderr)

%% intra-celltype statistics on fraction of IFN+ 
epi_stat = kruskalwallis(frac(:, 1:4));
fib_stat = kruskalwallis(frac(:, 5:9));
endo_stat = kruskalwallis(frac(:, 10:13));
immune_stat = kruskalwallis(frac(:, 14:23));

%% inter-celltype statistics on fraction of IFN+ 
epi = frac(:,1:4);
fib = frac(:, 5:9);
endo = frac(:, 10:13);
immune = frac(:, 14:23);

alpha = 0.05;
numTests = 7;
p_values = zeros(numTests, 1);
[~, p_values(1)] = ttest2(epi(:), fib(:), 'Tail', 'both'); 
[~, p_values(2)] = ttest2(fib(:), endo(:), 'Tail', 'both');
[~, p_values(3)] = ttest2(endo(:), immune(:), 'Tail', 'both');
[~, p_values(4)] = ttest2(fib(:), immune(:), 'Tail', 'both');
[~, p_values(5)] = ttest2(epi(:), immune(:), 'Tail', 'both');
[~, p_values(6)] = ttest2(epi(:), endo(:), 'Tail', 'both');
[~, p_values(7)] = ttest2(fib(:), endo(:), 'Tail', 'both');

% Apply Bonferroni Correction ---
bonferroni_alpha = alpha / numTests;
fprintf('Original Alpha: %.4f, Bonferroni Alpha: %.4f\n', alpha, bonferroni_alpha);

% Determine Significant Results ---
significant_results = p_values < bonferroni_alpha;
