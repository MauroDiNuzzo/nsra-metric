
function print_diagnostics(name,meas,pred)
    mae = mean(abs(pred-meas));
    score = nsra(meas,pred);
    fprintf('%s\n',name);
    score_check = nsra(meas,pred,Method="full-pairwise");
    if ~isnan(score) && (score ~= score_check)
        fprintf('  NSRA (ref): %.3f\n',score_check);
        warning("Score equivalence check failed.");
    end
    fprintf('  NSRA: %.3f\n',score);
    fprintf('  MAE : %.3f\n\n',mae);
end

clc
fprintf("Running NSRA unit tests...\n\n");

% Perfect ordering
meas = [2; 1; 0; -1; -2];
pred = [10; 5; 0; -5; -9];
assert(nsra(meas,pred) == 1.0);
print_diagnostics('Perfect ordering',meas,pred);

% Reversed ordering
pred = -meas;
assert(nsra(meas,pred) == 0.0);
print_diagnostics('Reverse ordering',meas,pred);

% All NULL
meas = zeros(10,1);
pred = randn(10,1);
assert(isnan(nsra(meas,pred)));
print_diagnostics('All null',meas,pred);

% Predicted ties
meas = [1; 0; -1];
pred = [0; 0; 0];
assert(abs(nsra(meas,pred)-0.5) < 1e-12);
print_diagnostics('Predicted ties',meas,pred);

% False positives below true signal
meas = [2; -2; zeros(20,1)];
pred = [10; -10; linspace(-1,1,20)'];
assert(nsra(meas,pred) > 0.9);
print_diagnostics('False positives (below true signal)',meas,pred);

% False positives interleaved
pred = [0.5; -0.5; linspace(-10,10,20)'];
assert(nsra(meas,pred) < 0.7);
print_diagnostics('False positives (interleaved)',meas,pred);

% False positives ranked safely
meas = [2; -2; zeros(20,1)];
pred = [10; -10; linspace(-1,1,20)'];
assert(nsra(meas,pred) > 0.9);
print_diagnostics('False positives (ranked safely)',meas,pred);

% False negatives
meas = [2; -2; zeros(10,1)];
pred = zeros(size(meas));
assert(abs(nsra(meas,pred) - 0.5) < 1e-12);
print_diagnostics('False negatives (collapsed to NULL)',meas,pred);

% Sparse single-gene perturbation (correct)
meas = [3; zeros(50,1)];
pred = [5; randn(50,1)];
assert(nsra(meas,pred) > 0.9);
print_diagnostics('Sparse perturbation (correct)',meas,pred);

% Sparse single-gene perturbation (wrong)
pred(1) = -5;
assert(nsra(meas,pred) < 0.1);
print_diagnostics('Sparse perturbation (wrong direction)',meas,pred);

% Equivalence and time-to-complete
rng(1);
meas = randn(20000,1);
pred = randn(20000,1);
N = 1;
timerVal = tic();
for i = 1:N
    s1 = nsra(meas,pred,Epsilon=0.1,TieScore=0.5,Method="full-pairwise");  % reference (slow)
end
full_pairwise_ttc = toc(timerVal);
timerVal = tic();
for i = 1:N
    s2 = nsra(meas,pred,Epsilon=0.1,TieScore=0.5,Method="rank-reduced"); % default (fast)
end
rank_reduced_ttc = toc(timerVal);
assert(abs(s1 - s2) < 1e-12);
timerVal = tic();
for i = 1:N
    s3 = nsra_mex(meas,pred, 0.1, 0.5);
end
compiled_ttc = toc(timerVal);
fprintf("[MATLAB] Elapsed times: %.3f (full-pairwise) vs %.3f (rank-reduced) vs %.3f (compiled)\n\n", ...
    full_pairwise_ttc,rank_reduced_ttc,compiled_ttc);

fprintf("\nAll tests passed successfully!\n")  ;