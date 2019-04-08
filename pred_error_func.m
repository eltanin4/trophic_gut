function [pred_error] = pred_error_func(x, f, b_real, mu_matrix, c2b_real, bp_out,i_fit_diet)
% pred_error_func is used to calculate the bacterial abundance error for
% each species existing in the patients data.
    % x is the nutrients existing in the diet we wants to predict. x_full
    % is used to construct the full size array of 2244 * 1 to make the
    % multiplication by mu_matrix and x_full be possible. Thus mu_matrix * 
    % x_full can give us the result of bacterial abundance.
    x_full = zeros(2244,1);
    x_full(i_fit_diet) = x;
    %x = x./sum(x);
    
    % pred_error is a vector rather than scalar because "lsqnonlin"
    % optimization function requires the vector output from the function rather
    % than scalar output. For other linear/nonlinear optimization method like 
    % fmincon, fminsearch, the scalar output like square error is required.
    ba_pred = mu_matrix * x_full ;
    %pred_error = (log10(ba_pred + 1e-6) - log10(b_real * sum(ba_pred)+1e-6)) ./ log10(b_real * sum(ba_pred)+1e-6);
    pred_error = (log10(ba_pred + 1e-6) - log10(b_real +1e-6)) ./ log10(b_real +1e-6);
    %pred_error = ((mu_matrix * x_full + 1e-6) - (b_real+1e-6)) ./ (b_real+1e-6);
    %pred_error = sum(pred_error .* pred_error);
    
end