function tests = custom_project_amplitudes_test
    tests = functiontests(localfunctions);
    %disp("lol")
end

function setupOnce(testCase)
    import matlab.unittest.TestCase
    import matlab.unittest.constraints.IsEqualTo
    import matlab.unittest.constraints.RelativeTolerance
    
    testCase.TestData.N = 50;
end

function testLegendrePolynomials(testCase)
    ff = @(theta) cos(theta);
    % endpoints = [0, pi];
    angles = linspace(0, pi, 100);
    values = ff(angles);
    %amplitudes = zeros(testCase.TestData.N, 1);
    %amplitudes(2) = 1;
    
    %[handles, ~, ~] = create_function_handle(testCase.TestData.polynomials, {}, {});
    
    actual_value = custom_project_amplitudes(angles, values, testCase.TestData.N, NaN, true);
    %fprintf("Test 1: Time %.3e seconds\n", timeit(@() project_amplitudes(ff, testCase.TestData.N, endpoints, NaN, true)));
    expected_value = [1; zeros(testCase.TestData.N - 1, 1)];
    
    %input = rand(); % Random but deterministic value
    
    %testCase.verifyThat(actual_handle(input), IsEqualTo(expected_handle(input), "Within", RelativeTolerance(1e-5)));
    verifyEqual(testCase, actual_value, expected_value, ...
         "AbsTol", 1e-5); %, sprintf("Input value: %.3f", input));
end

function testLegendrePolynomialsTwo(testCase)
    ff = @(theta) legendreP(5, cos(theta)) + legendreP(12, cos(theta));
    % endpoints = [0, pi];
    angles = linspace(0, pi, 2000);
    values = ff(angles);
    %amplitudes = zeros(testCase.TestData.N, 1);
    %amplitudes(2) = 1;
    
    %[handles, ~, ~] = create_function_handle(testCase.TestData.polynomials, {}, {});
    
    actual_value = custom_project_amplitudes(angles, values, testCase.TestData.N, NaN, true);
    %fprintf("Test 1: Time %.3e seconds\n", timeit(@() project_amplitudes(ff, testCase.TestData.N, endpoints, NaN, true)));
    expected_value = zeros(testCase.TestData.N, 1); expected_value([5, 12]) = 1;
    
    %input = rand(); % Random but deterministic value
    
    %testCase.verifyThat(actual_handle(input), IsEqualTo(expected_handle(input), "Within", RelativeTolerance(1e-5)));
    verifyEqual(testCase, actual_value, expected_value, ...
         "AbsTol", 1e-5, "RelTol", 1e-3); %, sprintf("Input value: %.3f", input));
end

function testSimpleLinearRandomFunctions(testCase)
    for ii = 1:25
        a = rand();
        b = a + rand();
        ff = @(theta) (cos(theta) - cos(a))/(cos(b) - cos(a));
        custom = custom_project_amplitudes([a, b], [0, 1], testCase.TestData.N, NaN, NaN);
        old    = project_amplitudes(ff, testCase.TestData.N, [a, b], NaN, 1);
        
        verifyEqual(testCase, custom, old, "AbsTol", 1e-5);
    end
end

function testCompoundLinearRandomFunctions(testCase)
    for ii = 1:25
        a = rand();
        b = a + rand();
        c = b + rand();
        angles = [a, b, c];
        values = [0, 1, 2];
        ff = @(theta) interp1(cos(angles), values, cos(theta), 'linear', 0);
        custom = custom_project_amplitudes(angles, values, testCase.TestData.N, NaN, NaN);
        old    = project_amplitudes(ff, testCase.TestData.N, [a, c], NaN, 1);
        
        verifyEqual(testCase, norm(custom-old), 0, "AbsTol", 1e-4);
        
    end
end