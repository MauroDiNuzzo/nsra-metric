# test_nsra.py
import numpy as np
import time
import warnings
from nsra import nsra
from nsra_cpp import nsra_cpp

def test_all():
    """Run all NSRA tests."""
    test_perfect_ordering()
    test_reversed_ordering()
    test_all_null()
    test_predicted_ties()
    test_false_positives_ranked_below_true_U()
    test_false_positives_interleaved()
    test_false_positives_ranked_safely()
    test_false_negatives_collapse_to_null()
    test_sparse_single_gene_correct()
    test_sparse_single_gene_wrong_direction()
    test_epsilon_effect()       
    test_dense_equivalence_random()
    # add other test functions here as needed
    print("\nAll tests passed successfully!")    

def print_diagnostics(name, meas, pred, epsilon=0.0):
    mae = np.mean(np.abs(pred - meas))
    score = nsra(meas, pred, epsilon=epsilon)
    print(f"{name}")
    score_check = nsra(meas, pred, method="full-pairwise")
    if not np.isnan(score) and (score != score_check):
        print(f"  NSRA (ref): {score_check:.3f}")
        warnings.warn("Score equivalence check failed.")
    print(f"  NSRA: {score:.3f}")
    print(f"  MAE : {mae:.3f}")
    print()

def test_perfect_ordering():
    meas = np.array([ 2,  1,  0, -1, -2])
    pred = np.array([10,  5,  0, -5, -9])
    assert nsra(meas, pred) == 1.0
    print_diagnostics("Perfect ordering",
                      meas, pred)    

def test_reversed_ordering():
    meas = np.array([ 2,  1,  0, -1, -2])
    pred = -meas
    assert nsra(meas, pred) == 0.0
    print_diagnostics("Reverse ordering",
                      meas, pred)    

def test_all_null():
    meas = np.zeros(10)
    pred = np.random.randn(10)
    assert np.isnan(nsra(meas, pred))
    print_diagnostics("All null",
                      meas, pred)

def test_predicted_ties():
    meas = np.array([ 1,  0, -1])
    pred = np.array([ 0,  0,  0])
    score = nsra(meas, pred, tie_score=0.5)
    assert np.isclose(score, 0.5)
    print_diagnostics("Predicted ties",
                      meas, pred)    

def test_false_positives_ranked_below_true_U():
    """
    Many false positives exist but are ranked below true U genes.
    NSRA should remain high.
    """
    meas = np.array([2.0, -2.0] + [0.0]*20)  # 1 U, 1 D, many NULL
    pred = np.array([10.0, -10.0] + list(np.linspace(-1, 1, 20)))
    score = nsra(meas, pred)
    assert score > 0.9
    print_diagnostics("False positives (ranked below true U)",
                      meas, pred)    

def test_false_positives_interleaved():
    """
    False positives interleave with true U/D.
    NSRA should penalize this.
    """
    meas = np.array([2.0, -2.0] + [0.0]*20)
    pred = np.array([0.5, -0.5] + list(np.linspace(-10, 10, 20)))
    score = nsra(meas, pred)
    assert score < 0.7
    print_diagnostics("False positives (interleaved)",
                      meas, pred)    

def test_false_positives_ranked_safely():
    meas = np.array([2.0, -2.0] + [0.0]*20)
    pred = np.array([10.0, -10.0] + list(np.linspace(-1, 1, 20)))
    score = nsra(meas, pred)
    assert score > 0.9
    print_diagnostics("False positives (ranked safely)",
                      meas, pred)

def test_false_negatives_collapse_to_null():
    """
    True U/D genes predicted near zero.
    """
    meas = np.array([2.0, -2.0] + [0.0]*10)
    pred = np.zeros_like(meas)
    score = nsra(meas, pred)
    assert np.isclose(score, 0.5)  # only ties
    print_diagnostics("False negatives (collapsed to NULL)",
                      meas, pred)

def test_sparse_single_gene_correct():
    """
    One true U gene, correctly ranked.
    """
    meas = np.array([3.0] + [0.0]*50)
    pred = np.array([5.0] + list(np.random.randn(50)))
    score = nsra(meas, pred)
    assert score > 0.9
    print_diagnostics("Sparse (single gene correct)",
                      meas, pred)    

def test_sparse_single_gene_wrong_direction():
    """
    One true U gene predicted as strongly down.
    """
    meas = np.array([3.0] + [0.0]*50)
    pred = np.array([-5.0] + list(np.random.randn(50)))
    score = nsra(meas, pred)
    assert score < 0.1
    print_diagnostics("Sparse (single gene wrong direction)",
                      meas, pred)

def test_epsilon_effect():
    meas = np.array([0.01, -0.01, 1.0, -1.0])
    pred = np.array([0.2, -0.2, 0.5, -0.5])
    s_small_eps = nsra(meas, pred, epsilon=0.0)
    s_large_eps = nsra(meas, pred, epsilon=0.1)
    assert s_large_eps >= s_small_eps 

def test_dense_equivalence_random():
    # Dense equivalence (small random)
    rng = np.random.default_rng(1)
    meas = rng.standard_normal(20000)
    pred = rng.standard_normal(20000)
    N = 1
    epsilon = 0.1
    tie_score = 0.5
    # dense version
    t0 = time.perf_counter()
    for _ in range(N):
        s1 = nsra(meas, pred, epsilon, tie_score, method="full-pairwise")
    full_pairwise_ttc = time.perf_counter() - t0
    # fast version
    t0 = time.perf_counter()
    for _ in range(N):
        s2 = nsra(meas, pred, epsilon, tie_score)
    rank_reduced_ttc = time.perf_counter() - t0
    # equivalence check
    assert abs(s1 - s2) < 1e-12
    # fast compiled version
    t0 = time.perf_counter()
    for _ in range(N):
        s3 = nsra_cpp(meas, pred, epsilon, tie_score)   
    compiled_ttc = time.perf_counter() - t0 
    print(f"[PYTHON] Elapsed times: {full_pairwise_ttc:.3f} (full-pairwise) vs {rank_reduced_ttc:.3f} (rank-reduced) vs {compiled_ttc:.3f} (compiled)\n")

if __name__ == "__main__":
    print("Running NSRA unit tests...\n")
    test_all()      