function [fn, fn_prime, fn_primeprime] = create_function_handle(LEGENDRE_POLYNOMIALS, ...
            LEGENDRE_DERIVATIVES, LEGENDRE_SECOND_DERIVATIVES)
   % Returns an array of handles that correspond to generators of legendre
   % polynomials 
   fn = @(amplitudes) combine_functions(amplitudes, LEGENDRE_POLYNOMIALS);
   fn_prime = @(amplitudes) combine_functions(amplitudes, LEGENDRE_DERIVATIVES);
   fn_primeprime = @(amplitudes) combine_functions(amplitudes, LEGENDRE_SECOND_DERIVATIVES);

end
function f = combine_functions(scalars, functionHandles)
    % This function returns a linear combination of functions with weights
    % given by scalars(ii). 
    if size(scalars, 2) > 1; scalars = scalars'; end
    f = @(x) myGeneralFunction(x, scalars, functionHandles);
    function y = myGeneralFunction(x, scalars, funHandles)
        if size(x, 1) > 1; x = x'; end
        scalars(1) = 0;
        y = cell2mat(cellfun(@(fcn) fcn(x), funHandles, 'UniformOutput', false));
        y = sum(scalars .* y, 1);
    end

end