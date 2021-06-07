function [list,numbers] = list_subband_flips(X_filt,Y_filt,Z_filt)
% X_filt --> either 'L' or 'H'
% Y_filt --> either 'L' or 'H'
% Z_filt --> either 'L' or 'H'

list = cell(24,1);
list{1} = [X_filt,Y_filt,Z_filt]; % g_0_0
list{2} = ['j',Z_filt,Y_filt,X_filt]; % g_0_pi/2_0
list{3} = ['j',X_filt,Y_filt,'j',Z_filt]; % g_0_pi_0
list{4} = [Z_filt,Y_filt,'j',X_filt]; % g_0_3pi/2_0
list{5} = [Y_filt,Z_filt,X_filt]; % g_pi/2_0_pi/2
list{6} = [Y_filt,'j',Z_filt,'j',X_filt]; % g_pi/2_0_3pi/2
list{7} = [Y_filt,'j',X_filt,Z_filt]; % g_pi/2_0_0
list{8} = ['j',X_filt,'j',Y_filt,Z_filt]; % g_pi_0_0
list{9} = ['j',Y_filt,X_filt,Z_filt];  % g_3pi/2_0_0
list{10} = ['j',Z_filt,'j',X_filt,Y_filt]; % g_0_pi/2_3pi/2
list{11} = ['j',Z_filt,'j',Y_filt,'j',X_filt]; % g_0_pi/2_pi
list{12} = ['j',Z_filt,X_filt,'j',Y_filt]; % g_0_pi/2_pi/2
list{13} = ['j',Y_filt,'j',X_filt,'j',Z_filt]; % g_pi/2_pi_0
list{14} = [X_filt,'j',Y_filt,'j',Z_filt]; % g_pi_pi_0
list{15} = [Y_filt,X_filt,'j',Z_filt]; % g_3pi/2_pi_0
list{16} = [Z_filt,'j',X_filt,'j',Y_filt]; % g_0_3pi/2_pi/2
list{17} = [Z_filt,'j',Y_filt,X_filt]; % g_0_3pi/2_pi
list{18} = [Z_filt,X_filt,Y_filt]; % g_0_3pi/2_3pi/2
list{19} = ['j',X_filt,Z_filt,Y_filt]; % g_pi_0_pi/2
list{20} = ['j',Y_filt,Z_filt,'j',X_filt]; % g_3pi/2_0_pi/2
list{21} = [X_filt,Z_filt,'j',Y_filt]; % g_0_0_pi/2
list{22} = ['j',X_filt,'j',Z_filt,'j',Y_filt];  % g_pi_0_3pi/2
list{23} = ['j',Y_filt,'j',Z_filt,X_filt]; % g_3pi/2_0_3pi/2
list{24} = [X_filt,'j',Z_filt,Y_filt]; % g_0_0_3pi/2
list = unique(list,'rows','stable');

list_new = strrep(list,'j','');
numbers = zeros(24,1);
for n = 1:24
    if strcmp(list_new{n},'LLL')
        numbers(n) = 1;
    elseif strcmp(list_new{n},'LLH')
        numbers(n) = 2;
    elseif strcmp(list_new{n},'LHL')
        numbers(n) = 3;
    elseif strcmp(list_new{n},'LHH')
        numbers(n) = 4;
    elseif strcmp(list_new{n},'HLL')
        numbers(n) = 5;
    elseif strcmp(list_new{n},'HLH')
        numbers(n) = 6;
    elseif strcmp(list_new{n},'HHL')
        numbers(n) = 7;
    elseif strcmp(list_new{n},'HHH')
        numbers(n) = 8;
    end
end

end