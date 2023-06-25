syms x;

N = 20;

f = matlabFunction(legendreP(N, x));

for ii = 1:100
   a = rand();
   f(a);
   legendreP(N, a);
end

% N = 20;
% 
% syms ang;
% LEGENDRE_POLYNOMIALS = cell(N, 1);
% LEGENDRE_DERIVATIVES = cell(N, 1);
% LEGENDRE_SECOND_DERIVATIVES = cell(N, 1);
% 
% for ii = 1:N
%     P = legendreP(ii, cos(ang));
%     LEGENDRE_POLYNOMIALS{ii} = matlabFunction(P);
%     LEGENDRE_DERIVATIVES{ii} = matlabFunction(diff(P));
%     LEGENDRE_SECOND_DERIVATIVES{ii} = matlabFunction(diff(P, 2));
% end
% [oscillation_handle, oscillation_handle_prime, ~] = create_function_handle(LEGENDRE_POLYNOMIALS, LEGENDRE_DERIVATIVES, LEGENDRE_SECOND_DERIVATIVES);
% clear ang;
% 
% 
% % To test whether this code below  is independent of the matlab symbolic
% % toolbox. Number of calls to mupadmex should be independent of the number
% % of iterations in the for loop below.
% for ii = 1:10000
%    f = oscillation_handle(rand(1, N));
%    
%    disp(f(rand()));
% 
% end