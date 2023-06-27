function tests = project_amplitudes_test
    tests = functiontests(localfunctions);
end

function setupOnce(testCase)
    import matlab.unittest.TestCase
    import matlab.unittest.constraints.IsEqualTo
    import matlab.unittest.constraints.RelativeTolerance
    
    testCase.TestData.N = 10;
    %N = testCase.TestData.N;
    %syms ang;
    %LEGENDRE_POLYNOMIALS = cell(N, 1);
    %LEGENDRE_DERIVATIVES = cell(N, 1);
%     LEGENDRE_SECOND_DERIVATIVES = cell(N, 1);
%     
%     for ii = 1:N
%         P = legendreP(ii, cos(ang));
%         LEGENDRE_POLYNOMIALS{ii} = matlabFunction(P);
%         LEGENDRE_DERIVATIVES{ii} = matlabFunction(diff(P));
%         LEGENDRE_SECOND_DERIVATIVES{ii} = matlabFunction(diff(P, 2));
%     end
%     
%     clear ang;
%     
%     testCase.TestData.polynomials = LEGENDRE_POLYNOMIALS;
%     %testCase.TestData.oscillation_handle_prime       = oscillation_handle_prime;
    %testCase.TestData.oscillation_handle_prime_prime = oscillation_handle_prime_prime;
end

function testLegendrePolynomials(testCase)
    ff = @(theta) cos(theta);
    endpoints = [0, pi];
    
    %amplitudes = zeros(testCase.TestData.N, 1);
    %amplitudes(2) = 1;
    
    %[handles, ~, ~] = create_function_handle(testCase.TestData.polynomials, {}, {});
    
    actual_value = project_amplitudes(ff, testCase.TestData.N, endpoints, NaN, true);
    expected_value = [1; zeros(testCase.TestData.N - 1, 1)];
    
    %input = rand(); % Random but deterministic value
    
    %testCase.verifyThat(actual_handle(input), IsEqualTo(expected_handle(input), "Within", RelativeTolerance(1e-5)));
    verifyEqual(testCase, actual_value, expected_value, ...
         "AbsTol", 1e-5); %, sprintf("Input value: %.3f", input));
end

% function testsimplePolynomialsTwo(testCase)
%     amplitudes = zeros(testCase.TestData.N, 1);
%     amplitudes(2:3) = [1, 1];
%     
%     [handles, ~, ~] = create_function_handle(testCase.TestData.polynomials, {}, {});
%     
%     actual_handle = handles(amplitudes);
%     expected_handle = @(ang) legendreP(2, cos(ang)) + legendreP(3, cos(ang));
%     
%     input = rand(); % Random but deterministic value
%     
%     %testCase.verifyThat(actual_handle(input), IsEqualTo(expected_handle(input), "Within", RelativeTolerance(1e-5)));
%     verifyEqual(testCase, actual_handle(input), expected_handle(input), ...
%         "RelTol", 1e-5, sprintf("Input value: %.3f", input));
% end
