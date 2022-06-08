% implementing bootstrapping for comparing two unequal size data matrices
% from my understanding:
% psuedo code:
% 1. Get anova F-value (thres)
% 2. Get a distribution of null F-value
% 3. The proportion of combinations greater than thres F-val is new "p-value"


% there's "unbalanced ANOVA" as well

function [out_pval,ref_pval] = compare_unequalSets(set1,set2,method)
    % set1 - N1 X 1
    % set2 - N2 X 1
    % Null distribution
    nboot = 2000;
    global group;
    N1 = length(set1); N2 = length(set2);
    group = [zeros(N1,1); ones(N2,1)];
    if(method==1 || method==3)
        d = [set1; set2];
        [refF,ref_pval] = f_stat(d);
        bootstat = bootstrp(nboot,@f_stat,d); %,group);
        out_pval = sum(bootstat > refF)/nboot;
        if(method==3)
            out_pval=ref_pval;
        end
    elseif(method==2)
        b_set1m = bootstrp(nboot,@nanmean,set1);
        b_set2m = bootstrp(nboot,@nanmean,set2);
        nboot_nnan = sum(~isnan(b_set2m));
        refDiff = nanmean(set1)-nanmean(set2);
        bootDiff = b_set1m - b_set2m;
        out_pval = nansum(bootDiff > refDiff)/nboot_nnan;
        ref_pval = []; % if unassigned, shoots error
    elseif(method==4) % KW test
        d = [set1; set2];
        valid = ~isnan(d);
        out_pval = kruskalwallis(d(valid),group(valid),'off');
        ref_pval = [];
    end
end

% another implementation with θ' = z - y, where z and y are, respectively,
% the means of two samples drawn from two distributions
% through this empirical distribution, it is possible to estimate the probability of occurrence 
% of the observed 8' under the null hypothesis-that is, how many 8'*s are equal or larger than
% 8' by chance-and therefore its level ofsignificance
% Ref: Di Nocera, F., & Ferlazzo, F. (2000). 
% Resampling approach to statistical inference: bootstrapping from event-related potentials data. 
% Behavior research methods, instruments, & computers, 32(1), 111-119.

function [out,p_val] = f_stat(d)
global group;
[~,tbl] = anova1(d,group,'off');
out = tbl{2,5};
p_val = tbl{2,6};
end

% Comment by Murty:
% By bootstrap, Supi might be referring to just reporting mean and std of
% data1 and data2 without any testing through bootstrapping (N~20000). This
% is to just support the p-val from unbalanced ANOVA. Your method is not
% widely known.
% Comment by Vinay:
% I know how to do it, but not exactly. Try referring to papers that have
% implemented that and following them.