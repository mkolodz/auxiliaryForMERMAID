%calculate the normalization coefficients:

function [c, g, history] = estimate_normalization_components(m_ij, s_ij, S_ij, i_crys, j_crys, n_crystals, n_crystals_one_side, max_iter, threshold)
% m_ij: Measured counts [LORs x 1]
% s_ij: Simulated counts [LORs x 1]
% S_ij: geometry lookup table (geometrically equivalent sets of LORs)
% i_crys, j_crys: crystal indices [LORs x 1]
% n_crystals: total number of crystals
% max_iter: maximum number of iterations
% threshold: NRMS convergence threshold (e.g., 1e-3)



n_lors = length(m_ij);

% Initialize
c = ones(n_crystals, 1);
g = ones(n_lors, 1);

% For convergence tracking
history = struct('delta_c', [], 'delta_g', []);

for iter = 1:max_iter
    c_old = c;
    g_old = g;

    
    % --- Update crystal efficiency c_i ---
    for i = 1:n_crystals
        idx_i = (i_crys == i);
        if all(~isnan(s_ij(idx_i))) && all(~isnan(m_ij(idx_i)))
            j_set = j_crys(idx_i);
            s_vals = s_ij(idx_i);
            g_vals = g(idx_i);
            c_j = c(j_set);
            m_vals = m_ij(idx_i);
    
            numerator = sum(m_vals .* c_j .* g_vals .* s_vals);
            denominator = sum((c_j.^2) .* (g_vals.^2) .* s_vals) + eps;
    
            c(i) = numerator / denominator;
        else 
            c(i) = -100;
        end
    end

    % --- Update crystal efficiency c_j ---
    for j = 1:n_crystals
        idx_j = (j_crys == j);
        if all(~isnan(s_ij(idx_j))) && all(~isnan(m_ij(idx_j)))
            i_set = i_crys(idx_j);
            s_vals = s_ij(idx_j);
            g_vals = g(idx_j);
            c_i = c(i_set);
            m_vals = m_ij(idx_j);
    
            numerator = sum(m_vals .* c_i .* g_vals .* s_vals);
            denominator = sum((c_i.^2) .* (g_vals.^2) .* s_vals) + eps;
    
            c(j) = numerator / denominator;
        else
            c(j) = -100;
        end
    end

    % % --- Update geometric factors g_ij ---
    % %if all(~isnan(c(i_crys))) && all(~isnan(c(j_crys)))
    % g = m_ij ./ (c(i_crys) .* c(j_crys) .* s_ij + eps);
    % %else
    % %    g = -100;
    % %end
    %g=-100;
    l_indices_1 = S_ij{2,257}; % S_ij: cell array with list of LOR indices for each i,j pair
    disp(l_indices_1);
    disp("aaa"+ l_indices_1(1,2));
    for i = 1:n_crystals_one_side
        for j = (n_crystals_one_side+1):n_crystals
            % Find indices l in S_ij for LOR-ij (geometrically equivalent)
            l_indices = S_ij{i,j-256}; % S_ij: cell array with list of LOR indices for each i,j pair
    
            num_sum = 0; % numerator
            den_sum = 0; % denominator
            
            for l = l_indices

                l_singleID = l(i,1)*256+l(i,2);
                %l(i,:) - one index pair; translate this to 1D: l(i,1)-x,
                %l(i,2)-y, actuallyNeededID = x*256+y --->edit this and
                %correct and validate
                num_sum = num_sum + m_ij(l_singleID) * c_old(i) * c_old(j) * s_ij(l_singleID);
                den_sum = den_sum + (c_old(i) * c_old(j) * s_ij(l_singleID))^2;
            end
    
            if den_sum ~= 0
                g(i,j) = num_sum / den_sum;
            else
                g(i,j) = 1; % fallback (or NaN)
            end
        end
    end




    % --- Convergence criteria ---
    delta_c = sqrt(sum((c - c_old).^2) / n_crystals);
    delta_g = sqrt(sum((g - g_old).^2) / n_lors);

    history.delta_c(end+1) = delta_c;
    history.delta_g(end+1) = delta_g;

    fprintf('Iter %d: Δc = %.3e, Δg = %.3e\n', iter, delta_c, delta_g);

    if delta_c < threshold && delta_g < threshold
        disp('Converged.');
        break;
    end
end
end
