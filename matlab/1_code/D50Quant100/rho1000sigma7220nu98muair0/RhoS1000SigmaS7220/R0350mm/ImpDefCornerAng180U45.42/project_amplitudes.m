function coefs = project_amplitudes(ff, harmonics_qtt, endpoints, PROBLEM_CONSTANTS, flag) % PROBLEM CONSTANTS used to be here 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% project_amplitudes.m - Projects function onto legendre polynomials
%
% Calculates the nth amplitude coefficient of function ff, where n <=
% harmonics_qtt. Note that the zeroth mode is not calculated.
% THe following formula applies
% coefs(n) = (2n+1)/2 \int_{0}^{pi} ff(w) * sin(w) * P_l(cos(w)) dw
% where Pl is the l-th legendre Polynomial. (see
% https://en.wikipedia.org/wiki/Legendre_polynomials#Rodrigues'_formula_and_other_explicit_formulas) 

% Arguments:
% - ff: Function handle that is going to be projected.
% - harmonics_qtt: number of harmonics to project the function to
% - endpoints: The function may be set to zero outside of the endpoints
% - PROBLEM_CONSTANTS: Struct that has a field with the precomputed weights
% for the Curtis-Clensahw Quadrature.
% - flag: Boolean that flags whether to use Curtis-CLenshaw. If set to
% false, it uses Curtis Clenshaw.
%
% Outputs:
% - coefs: Array of coefficients that correspond to the amplitudes.

% EXAMPLES
%
% Written by: Elvis Aguero- 02/01/2023
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if flag

            coefs = ((1:harmonics_qtt) + 0.5)' .* ...
                integral(@(theta) (ff(theta) .* sin(theta)) .* collectPl(harmonics_qtt, cos(theta)), ...
                    endpoints(1), endpoints(2), 'RelTol', 1e-5, 'AbsTol', 1e-5, 'ArrayValued', 1);
            assert(min(size(coefs)) == 1);
            %fprintf("The norm is: %.4f", norm(coefs
        %end
    else
        warning("Im deprecated!");
        if endpoints(2) > PROBLEM_CONSTANTS.nodes(1) || endpoints(1) < PROBLEM_CONSTANTS.nodes(end)
            [PROBLEM_CONSTANTS.nodes, PROBLEM_CONSTANTS.weights] = fclencurt(2^19+1, endpoints(1), endpoints(2));
            warning("Integration vector recalculated");
        end
        %loc = 1; n = length(PROBLEM_CONSTANTS.nodes);
        
        loc = knnsearch(PROBLEM_CONSTANTS.nodes, endpoints(1));
        disp(loc);
%         while PROBLEM_CONSTANTS.nodes(loc) > endpoints(1) && 2 * loc < n
%             loc = 2 * loc;
%         end
%         if PROBLEM_CONSTANTS.nodes(loc) > endpoints(1)
%             locmax = n;
%             locmin = loc;
%         else
%             locmin = loc/2;
%             locmax = loc;
%         end
% 
%         while locmax - locmin > 1
%             idxmed = round((locmax + locmin)/2);
%             if PROBLEM_CONSTANTS.nodes(idxmed) > endpoints(1)
%                 locmin = idxmed;
%             else
%                 locmax = idxmed;
%             end
%         end
        nodes   = PROBLEM_CONSTANTS.nodes(1:loc);
        weights = PROBLEM_CONSTANTS.weights(1:loc);
        if loc < 1000; warning(srpintf("Too few evaluations for numerical integration: %d", loc)); end
        
%         coefs = zeros(harmonics_qtt, 1)
%         for idx = 1:harmonics_qtt
%             funcs = my_legendre(idx, cos(theta))
%         end
        leg_matrix = collectPl(harmonics_qtt, cos(nodes));
        integral_results =  ff(nodes) .* sin(nodes) .* weights;
        integral_results = leg_matrix .* (integral_results');
        coefs = sum(((1:harmonics_qtt + 1/2))' .* integral_results, 2);
        
%         coefs = arrayfun(@(idx) ...
%             (2*idx+1)/2 * dot(integral(idx, :), weights), 1:harmonics_qtt);
    end
end