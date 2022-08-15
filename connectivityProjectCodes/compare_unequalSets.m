function out_pval = compare_unequalSets(set1,set2)
    % set1 - N1 X 1
    % set2 - N2 X 1
    N1 = length(set1); N2 = length(set2);
    group = [zeros(N1,1); ones(N2,1)];
    d = [set1; set2];
    valid = ~isnan(d);
    out_pval = kruskalwallis(d(valid),group(valid),'off');
end

