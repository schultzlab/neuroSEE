function [ spk_trials, phi_trials ] = div_spikes_intro_trials( activespk, activet, activephi )

p = activephi;
dthr = 20; % threshold for finding laps
Ncells = size(activespk,1);
for ii = 1:Ncells
    % find the delineations for the video: find t = 0
    idx_file = find(diff(activet) < 0);
    idx_file = [0; idx_file; numel(activet)] +1; 
    s = activespk(ii,:);
    Ntrial = 1;
    
    for jj = 1:numel(idx_file)-1
        % find the delineations per trial (i.e. loop)
        p_tr = p(idx_file(jj):idx_file(jj+1)-1);
        s_tr = s(idx_file(jj):idx_file(jj+1)-1);
        idx_tr = find( abs(diff(p_tr)) > dthr );
        for k = numel(idx_tr):-1:2
            if (idx_tr(k) - idx_tr(k-1)) <= 20 
                idx_tr(k) = 0;
            end
        end
        idx_tr = idx_tr( idx_tr > 0 );
        if numel(idx_tr)==1, idx_tr = [idx_tr; numel(p_tr)]; end
        
        for k = 1:numel(idx_tr)-1
            phi_trials{Ntrial} = p_tr(idx_tr(k)+1:idx_tr(k+1));
            spk_trials{ii}{Ntrial} = s_tr(idx_tr(k)+1:idx_tr(k+1));
            Ntrial = Ntrial + 1;
        end
    end
end