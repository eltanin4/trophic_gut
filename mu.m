function [m2b_total, m2m_total, m2m_layer] = mu(f, m2b, b2m, numLayer_max)
% This mu function is used to calculate the m2b_total which corresponds 
% to the conversion matrix from nutrient intake to bacterial abundance and m2m_total 
% which corresponds to the conversion matrix from nutrient intake to final metabolome.
% m2m_layer transform the nutrient intake into nutrient at each layer

    % matrices save metabolites' contribution to bacteria in each trophic layer
    m2m_layer = zeros(2244,2244,numLayer_max);  
    % matricx save metabolites' contribution to bacteria in all trophic layer
    m2b_total = zeros(2244,2244);  
    
    s_step =  b2m * m2b';  % s_step is the conversion matrix of each layer
    s_step_ii = eye(2244,2244);
    f_mul = repmat(f, 1, 2244);
    for ii = 1:numLayer_max
        % m2b_total is a series made of s_step of each layer
        m2b_total = m2b_total + f_mul.^(ii-1) .* s_step_ii; 
        m2m_layer(:,:,ii) = f_mul.^(ii-1) .* s_step_ii;
        s_step_ii = s_step_ii * s_step;
    end    
    m2m_total = m2b_total;
    m2b_total = (1-f_mul) .* m2b' * m2b_total;  % m2b_total has an extra multiplication of m2b and (1-f).
end