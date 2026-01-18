addpath("../matlab/");
%%
run('../tests/matlab/tests.m'); %[output:177f396e]

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline"}
%---
%[output:177f396e]
%   data: {"dataType":"text","outputData":{"text":"Running NSRA unit tests...\n\nPerfect ordering\n  NSRA: 1.000\n  MAE : 4.600\n\nReverse ordering\n  NSRA: 0.000\n  MAE : 2.400\n\nAll null\n  NSRA: NaN\n  MAE : 0.785\n\nPredicted ties\n  NSRA: 0.500\n  MAE : 0.667\n\nFalse positives (below true signal)\n  NSRA: 1.000\n  MAE : 1.206\n\nFalse positives (interleaved)\n  NSRA: 0.512\n  MAE : 4.921\n\nFalse positives (ranked safely)\n  NSRA: 1.000\n  MAE : 1.206\n\nFalse negatives (collapsed to NULL)\n  NSRA: 0.500\n  MAE : 0.333\n\nSparse perturbation (correct)\n  NSRA: 1.000\n  MAE : 0.817\n\nSparse perturbation (wrong direction)\n  NSRA: 0.000\n  MAE : 0.934\n\n[MATLAB] Elapsed times: 13.317 (full-pairwise) vs 0.005 (rank-reduced) vs 1.972 (compiled)\n\n\nAll tests passed successfully!\n","truncated":false}}
%---
