function tests = zs_from_spherical_unittest
    tests = functiontests(localfunctions);
end

function setupOnce(testCase)
    testCase.TestData.N = 3;
    N = testCase.TestData.N;
    syms ang;
    LEGENDRE_POLYNOMIALS = cell(N, 1);
    LEGENDRE_DERIVATIVES = cell(N, 1);
    LEGENDRE_SECOND_DERIVATIVES = cell(N, 1);
    
    for ii = 1:N
        P = legendreP(ii, cos(ang));
        LEGENDRE_POLYNOMIALS{ii} = matlabFunction(P);
        LEGENDRE_DERIVATIVES{ii} = matlabFunction(diff(P));
        LEGENDRE_SECOND_DERIVATIVES{ii} = matlabFunction(diff(P, 2));
    end
    
    [oscillation_handle, oscillation_handle_prime, oscillation_handle_prime_prime] = ...
        create_function_handle(LEGENDRE_POLYNOMIALS, LEGENDRE_DERIVATIVES, LEGENDRE_SECOND_DERIVATIVES);
    clear ang;
    
    testCase.TestData.oscillation_handle             = oscillation_handle;
    testCase.TestData.oscillation_handle_prime       = oscillation_handle_prime;
    testCase.TestData.oscillation_handle_prime_prime = oscillation_handle_prime_prime;
end

function test_both_frameworks(testCase)
    theta = pi * rand();
    amplitudes = rand(testCase.TestData.N, 1);
    zeta = testCase.TestData.oscillation_handle(amplitudes);
    actSolution = zs_from_spherical(theta, zeta);
    expSolution = zsoftheta(theta, amplitudes(2), amplitudes(3));
    
    verifyEqual(testCase, actSolution, expSolution, ...
        "RelTol", 1e-5, sprintf("Input: %.3f", theta));

end