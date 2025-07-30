function DTNnew345 = DTNVectorized(nr, D)    

rn = 0:nr+1;
dr = D/(2*nr);

refp = 10;
num_batches = 20; % Number of batches to reduce memory footprint

%Finding the DTN operator
drp = dr/refp;
numer = ceil(pi*D/drp);
if mod(numer,2)==1
    numer = numer+1;
end
dtheta = 2*pi/numer;%use pi/even number
DTNnew345=zeros(nr,nr);


%% Integrating away from the singularity
k = 1;
% Define the ii range as a vector and compute idx1 for all values of ii at once
i_valsk1 = 2:(rn(k) + nr + 1);
idx1_vals = round(i_valsk1);

% Preallocate arrays for computed values for each term
term_0 = (-i_valsk1.^2 / 2 - i_valsk1 - 1/3) .* log((i_valsk1 + 1) ./ i_valsk1) ...
         + (i_valsk1.^2 / 6 + i_valsk1 / 2 + 1/3) .* (1 - i_valsk1 ./ (i_valsk1 + 1)) ...
         - (i_valsk1 + 0.5) / 6 + i_valsk1 / 2 + 1/2;

term_1 = (3 * i_valsk1.^2 / 2 + 2 * i_valsk1 - 1/2) .* log((i_valsk1 + 1) ./ i_valsk1) ...
         + (-i_valsk1.^2 / 2 - i_valsk1 + 1/2 + 1 ./ i_valsk1) .* (1 - i_valsk1 ./ (i_valsk1 + 1)) ...
         + (i_valsk1 + 0.5) / 2 - 3 * i_valsk1 / 2 - 1;

term_2 = (-3 * i_valsk1.^2 / 2 - i_valsk1 + 1) .* log((i_valsk1 + 1) ./ i_valsk1) ...
         + (i_valsk1.^2 / 2 + i_valsk1 / 2 - 1) .* (1 - i_valsk1 ./ (i_valsk1 + 1)) ...
         - (i_valsk1 + 0.5) / 2 + 3 * i_valsk1 / 2 + 1/2;

term_3 = (i_valsk1.^2 / 2 - 1/6) .* log((i_valsk1 + 1) ./ i_valsk1) ...
         + (-i_valsk1.^2 + 1) / 6 .* (1 - i_valsk1 ./ (i_valsk1 + 1)) ...
         + (i_valsk1 + 0.5) / 6 - i_valsk1 / 2;


% Collect indices and values for accumarray
indices_0 = idx1_vals(idx1_vals < nr + 1);
values_0 = term_0(idx1_vals < nr + 1);

indices_1 = idx1_vals(idx1_vals < nr) + 1;
values_1 = term_1(idx1_vals < nr);

indices_2 = idx1_vals(idx1_vals < nr - 1) + 2;
values_2 = term_2(idx1_vals < nr - 1);

indices_3 = idx1_vals(idx1_vals < nr - 2) + 3;
values_3 = term_3(idx1_vals < nr - 2);

% Apply the update to DTNnew345 for row k
DTNnew345(k, :) = DTNnew345(k, :) ...
     - accumarray(indices_0.', values_0.', [nr, 1]).' ...
     - accumarray(indices_1.', values_1.', [nr, 1]).' ...
     - accumarray(indices_2.', values_2.', [nr, 1]).' ...
     - accumarray(indices_3.', values_3.', [nr, 1]).';

DTNnew345(1,1) = DTNnew345(1,1) + 1/2;
DTNnew345(1,:) = DTNnew345(1,:)/dr;
%Integrating the vincinity of the origin (= sing)
DTNnew345(1,1)   = DTNnew345(k,1) +209/(54*dr);%- (-15+    4)/(6*dr);
DTNnew345(1,2)   = DTNnew345(k,2) -29/(6*dr);%+ (-16+4/3*4)/(6*dr);
DTNnew345(1,3)   = DTNnew345(k,3) + 7/(6*dr);%+ (  1-1/3*4)/(6*dr);
DTNnew345(1,4)   = DTNnew345(k,4) - 11/(54*dr);%- (-15+    4)/(6*dr);


%% Integrating the vincinity of the origin
% Define vector of i values and compute Kern for all values
% Define i values and compute Kern
i_vals21 = (1:2 * refp)';  % Column vector for i values
i_vals22 = ((2 * refp + 1) : ((rn(2) + nr) * refp))';  % Column vector for ii
i_vals31 = (1:2 * refp)';  % Column vector for i values
i_vals32 = (2 * refp + 1 : (rn(3) + nr) * refp)';  % Column vector for ii
i_valsk1 = (1:2*refp).';    % Column vector for ii (transposed to column)

% Computing kernels
Kern21 = 2 * (1 ./ (i_vals21 - 1/2) - 1 ./ (i_vals21 + 1/2));  
Kern22 = 2 * (1 ./ (i_vals22 - 1/2) - 1 ./ (i_vals22 + 1/2));  
Kern31 = 2 * (1 ./ (i_vals31 - 1/2) - 1 ./ (i_vals31 + 1/2));  
Kern32 = 2 * (1 ./ (i_vals32 - 1/2) - 1 ./ (i_vals32 + 1/2));  
Kernk1 = 2 * (1 ./ (i_valsk1 - 1/2) - 1 ./ (i_valsk1 + 1/2));  % Column vector

% Define l values and trigonometric terms
l_vals = dtheta / 2 : dtheta : pi - dtheta / 4;
total_l_vals = length(l_vals);
batch_size = ceil(total_l_vals / num_batches);  % Batch size for l_vals

% Initialize cosine and sine vectors for l_vals
cos_l = cos(l_vals);
sin_l = sin(l_vals);
% Initialize Line update vector
LineUpdate2 = zeros(1, nr);
LineUpdate3 = zeros(1, nr);

% Define @sliceAndSum function for accumulating updates
sliceAndSum = @(expression, condition, idx, nr) accumarray(idx(condition), expression(condition), [nr, 1])';
% Loop over batches of l_vals
checkpoints = (1:9)/10;
for pp = 1:length(checkpoints)
[~, checkpoints(pp)] = min(abs((4:nr) - nr*(checkpoints(pp) + .01))); 
end

for batch = 1:num_batches
    k = 2;
    % Determine the range for this batch in l_vals
    idx_start = (batch - 1) * batch_size + 1;
    idx_end = min(batch * batch_size, total_l_vals);
    cos_l_batch = cos_l(idx_start:idx_end);
    sin_l_batch = sin_l(idx_start:idx_end);

    % Compute radn, x1, and posr for this batch
    radn = abs(sqrt((rn(k) + i_vals21 * cos_l_batch / refp).^2 + (i_vals21 * sin_l_batch / refp).^2));
    x1 = i_vals21 * cos_l_batch / refp;
    posr = radn - rn(k);

    % Vectorized updates for DTNnew345 for this batch
    DTNnew345(k, k) = DTNnew345(k, k) - sum((( -6 - 3 * posr + 9 * posr.^2 + 3 * posr.^3 - 3 * posr.^4) .* posr + 6 * x1) / 72 .* Kern21, 'all');
    DTNnew345(k, k-1) = DTNnew345(k, k-1) - sum(((-88 + 48 * posr + 62 * posr.^2 - 12 * posr.^3 - 10 * posr.^4) .* posr + 88 * x1) / 72 .* Kern21, 'all');
    DTNnew345(k, k) = DTNnew345(k, k) - sum((( 72 - 90 * posr - 90 * posr.^2 + 18 * posr.^3 + 18 * posr.^4) .* posr - 72 * x1) / 72 .* Kern21, 'all');

    % Condition k < nr - 0.5
    if k < nr - 0.5
        DTNnew345(k, k+1) = DTNnew345(k, k+1) - sum((( 24 + 48 * posr + 18 * posr.^2 - 12 * posr.^3 - 6 * posr.^4) .* posr - 24 * x1) / 72 .* Kern21, 'all');
    end
    
    % Condition k < nr - 1.5
    if k < nr - 1.5
        DTNnew345(k, k+2) = DTNnew345(k, k+2) - sum(((-2 - 3 * posr + posr.^2 + 3 * posr.^3 + posr.^4) .* posr + 2 * x1) / 72 .* Kern21, 'all');
    end

    % Now the second part
    % Compute radn, idxs, and weights for this batch
    
    radn = abs(sqrt((rn(k) + i_vals22 .* cos_l_batch / refp).^2 + (i_vals22 .* sin_l_batch / refp).^2));
    idxs = floor(radn);   % Floor values of radn
    %full_counter = full_counter + numel(radn);
    w1 = min(max(0, radn - idxs), 1);  % Clipped w1 values
    
    w2 = w1.^2;  % Square of w1 values
    w3 = w1.^3;  % Cube of w1 values

    % Apply conditions and accumulate updates for LineUpdate for this batch
    cond1 = idxs < 0.5; 
    LineUpdate2 = LineUpdate2 + sliceAndSum(-(3 * w3 / 4 - 7 * w2 / 4 + 1) .* Kern22, cond1, idxs + 1, nr);
    LineUpdate2 = LineUpdate2 + sliceAndSum(-(-w3 + 2 * w2) .* Kern22, cond1, idxs + 2, nr);
    LineUpdate2 = LineUpdate2 + sliceAndSum(-(w3 - w2) / 4 .* Kern22, cond1, idxs + 3, nr);

    cond2 = (idxs >= 0.5) & (idxs < nr); 
    LineUpdate2 = LineUpdate2 + sliceAndSum(-(-w3 / 6 + w2 / 2 - w1 / 3) .* Kern22, cond2, idxs, nr);
    LineUpdate2 = LineUpdate2 + sliceAndSum(-(w3 / 2 - w2 - w1 / 2 + 1) .* Kern22, cond2, idxs + 1, nr);
    %full_counter = full_counter + sum(sliceAndSum(-(w3 / 2 - w2 - w1 / 2 + 1) .* Kern22, cond2, idxs + 1, nr));
    
    cond3 = (idxs < nr - 1) & (idxs >= 0.5); 
    cond4 = (idxs < nr - 2) & (idxs >= 0.5); 
    LineUpdate2 = LineUpdate2 + sliceAndSum(-(-w3 / 2 + w2 / 2 + w1) .* Kern22, cond3, idxs + 2, nr);
    LineUpdate2 = LineUpdate2 + sliceAndSum(-(w3 - w1) / 6 .* Kern22, cond4, idxs + 3, nr);
    %DTNnew345(k, :) = (DTNnew345(k, :) + LineUpdate2);
    %full_counter = full_counter + LineUpdate2(end);

    k = 3;
    % Compute radn, x1, and posr for this batch
    radn_batch = abs(sqrt((rn(k) + i_vals31 .* cos_l_batch / refp).^2 + (i_vals31 .* sin_l_batch / refp).^2));
    x1_batch = i_vals31 .* cos_l_batch / refp;
    posr_batch = radn_batch - rn(k);

    % Update DTNnew345 elements based on the batch calculations
    DTNnew345(k, k-2) = DTNnew345(k, k-2) - sum(((124 - 12 * posr_batch - 149 * posr_batch.^2 + 12 * posr_batch.^3 + 25 * posr_batch.^4) .* posr_batch - 124 * x1_batch) / 288 .* Kern31, 'all');
    DTNnew345(k, k-1) = DTNnew345(k, k-1) - sum(((-384 + 192 * posr_batch + 288 * posr_batch.^2 - 48 * posr_batch.^3 - 48 * posr_batch.^4) .* posr_batch + 384 * x1_batch) / 288 .* Kern31, 'all');
    DTNnew345(k, k)   = DTNnew345(k, k)   + sum(((-4 + 10 * posr_batch + 5 * posr_batch.^2 - 2 * posr_batch.^3 - posr_batch.^4) .* posr_batch + 4 * x1_batch) / 8 .* Kern31, 'all');
    
    % Condition k < nr - 0.5
    if k < nr - 0.5
        DTNnew345(k, k+1) = DTNnew345(k, k+1) - sum(((128 + 192 * posr_batch + 32 * posr_batch.^2 - 48 * posr_batch.^3 - 16 * posr_batch.^4) .* posr_batch - 128 * x1_batch) / 288 .* Kern31, 'all');
    end
    
    % Condition k < nr - 1.5
    if k < nr - 1.5
        DTNnew345(k, k+2) = DTNnew345(k, k+2) - sum(((-12 - 12 * posr_batch + 9 * posr_batch.^2 + 12 * posr_batch.^3 + 3 * posr_batch.^4) .* posr_batch + 12 * x1_batch) / 288 .* Kern31, 'all');
    end


    % Second part
    % Calculate radn, indices, and weights
    radn = abs(sqrt((rn(k) + i_vals32 .* cos_l_batch / refp).^2 + (i_vals32 .* sin_l_batch / refp).^2));
    idxs = floor(radn);   % Floor values of radn
    w1 = min(max(0, radn - idxs), 1);  % Clipped w1 values
    w2 = w1.^2;  % Square of w1 values
    w3 = w1.^3;  % Cube of w1 values
    
    % Apply conditions and accumulate updates for Line
    cond1 = idxs < 0.5;
    LineUpdate3 = LineUpdate3 + sliceAndSum(-(3 * w3 / 4 - 7 * w2 / 4 + 1) .* Kern32, cond1, idxs + 1, nr);
    LineUpdate3 = LineUpdate3 + sliceAndSum(-(-w3 + 2 * w2) .* Kern32, cond1, idxs + 2, nr);
    LineUpdate3 = LineUpdate3 + sliceAndSum(-(w3 - w2) / 4 .* Kern32, cond1, idxs + 3, nr);
    
    cond2 = (idxs >= 0.5) & (idxs < nr);
    LineUpdate3 = LineUpdate3 + sliceAndSum(-(-w3 / 6 + w2 / 2 - w1 / 3) .* Kern32, cond2, idxs, nr);
    LineUpdate3 = LineUpdate3 + sliceAndSum(-(w3 / 2 - w2 - w1 / 2 + 1) .* Kern32, cond2, idxs + 1, nr);
    
    cond3 = (idxs < nr - 1) & (idxs >= 0.5);
    cond4 = (idxs < nr - 2) & (idxs >= 0.5);
    LineUpdate3 = LineUpdate3 + sliceAndSum(-(-w3 / 2 + w2 / 2 + w1) .* Kern32, cond3, idxs + 2, nr);
    LineUpdate3 = LineUpdate3 + sliceAndSum(-(w3 - w1) / 6 .* Kern32, cond4, idxs + 3, nr);



    %% Now k >= 4
    parfor k = 4:nr
        if ismember(k, checkpoints); disp([batch, k/nr]); end
        Line = zeros(1,nr);
        
        radn_vals = abs(sqrt((rn(k) + i_valsk1 .* cos_l_batch / refp).^2 + (i_valsk1 .* sin_l_batch / refp).^2));
        x1_vals = i_valsk1 .* cos_l_batch / refp;  % Matrix of x1 values
        posr_vals = radn_vals - rn(k);     % Matrix of posr values

        % Vectorized update for Line(k-2)
        Line(k-2) = Line(k-2) - sum(((2 - posr_vals - 2 * posr_vals.^2 + posr_vals.^3) .* posr_vals - 2 * x1_vals) ./ 24 .* Kernk1, 'all');
        
        % Vectorized update for Line(k-1)
        Line(k-1) = Line(k-1) - sum((4 * (-4 + 4 * posr_vals + posr_vals.^2 - posr_vals.^3) .* posr_vals + 16 * x1_vals) ./ 24 .* Kernk1, 'all');
        
        % Vectorized update for Line(k)
        Line(k) = Line(k) - sum((posr_vals.^2 - 5) .* posr_vals.^2 / 4 .* Kernk1, 'all');
        
        % Conditional vectorized updates for Line(k+1) and Line(k+2)
        if k < nr - 0.5
            Line(k+1) = Line(k+1) - sum((4 * (4 + 4 * posr_vals - posr_vals.^2 - posr_vals.^3) .* posr_vals - 16 * x1_vals) ./ 24 .* Kernk1, 'all');
            if k < nr - 1.5
                Line(k+2) = Line(k+2) - sum(((-2 - posr_vals + 2 * posr_vals.^2 + posr_vals.^3) .* posr_vals + 2 * x1_vals) ./ 24 .* Kernk1, 'all');
            end
        end
    
        i_valsk2 = ((2*refp+1):((rn(k)+nr)*refp))';
        
        % Compute Kernels and radn values for all combinations of ii and l
        Kernk2 = 2 * (1 ./ (i_valsk2 - 1/2) - 1 ./ (i_valsk2 + 1/2));        
        
        radn = abs(sqrt((rn(k) + i_valsk2 .* cos_l_batch / refp).^2 + (i_valsk2 .* sin_l_batch / refp).^2));
        idxs = floor(radn); 
        w1 = min(max(0, radn - idxs), 1); w2 = w1.^2; w3 = w1.^3;
        
        % Initialize Line vector update values
        LineUpdatek = zeros(1, nr);
        
        % Accumulate updates for each index in Line based on conditions
        cond1 = idxs < 0.5; 
        
        LineUpdatek = LineUpdatek + sliceAndSum(-(3 * w3 / 4 - 7 * w2 / 4 + 1) .* Kernk2, cond1, idxs+1, nr);
        LineUpdatek = LineUpdatek + sliceAndSum(-(-w3 + 2 * w2)  .* Kernk2, cond1, idxs+2, nr);
        LineUpdatek = LineUpdatek + sliceAndSum(-(w3 - w2) / 4 .* Kernk2, cond1, idxs+3, nr);
        
        
        cond2 = (idxs >= 0.5) & (idxs < nr); 
        
        LineUpdatek = LineUpdatek + sliceAndSum(-(-w3 / 6 + w2 / 2 - w1/3) .* Kernk2, cond2, idxs, nr);
        LineUpdatek = LineUpdatek + sliceAndSum(-(w3 / 2 - w2 - w1 / 2 + 1) .* Kernk2, cond2, idxs+1, nr);
        
        cond3 = (idxs < nr - 1) & (idxs >= 0.5); 
        cond4 = (idxs < nr - 2) & (idxs >= 0.5); 
        LineUpdatek = LineUpdatek + sliceAndSum(-(-w3 / 2 + w2 / 2 + w1) .* Kernk2, cond3, idxs+2, nr);
        LineUpdatek = LineUpdatek + sliceAndSum(-(w3 - w1) / 6 .* Kernk2, cond4, idxs+3, nr);
        
        Line = Line + LineUpdatek;
        
        % Final adjustments to Line
        Line = dtheta / (2 * pi * drp) * Line;
        DTNnew345(k,:) = DTNnew345(k, :) + Line;
    end

end
%disp(full_counter);
k = 2;
DTNnew345(k, :) = dtheta / (2 * pi * drp) * (DTNnew345(k, :) + LineUpdate2);
DTNnew345(k, k) = DTNnew345(k, k) + 2 / (4 * dr + drp);
k = 3;
DTNnew345(k, :) = dtheta / (2 * pi * drp) * (DTNnew345(k, :) + LineUpdate3);
DTNnew345(k, k) = DTNnew345(k, k) + 2 / (4 * dr + drp);

I = (eye(nr)==1); I(1:3, 1:3) = 0;
DTNnew345(I) = DTNnew345(I) + 2/(4*dr+drp);

save(['DTNnew345nr',num2str(nr),'D',num2str(D),'refp',num2str(refp),'.mat'],'DTNnew345')

end
