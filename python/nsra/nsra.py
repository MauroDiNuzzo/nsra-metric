# nsra.py
import numpy as np

def nsra(meas_deltas, pred_deltas, epsilon=0.0, tie_score=0.5, method="rank-reduced"):
    """
    Scalable Null-Stratified Rank Accuracy (NSRA)
    
    Parameters
    ----------
    meas_deltas : (G,) array
        Measured signed deltas vs control
    pred_deltas : (G,) array
        Predicted signed deltas vs control
    epsilon : float
        Soft sign tolerance
    tie_score : float
        Score for ties / ambiguous comparisons (default 0.5)

    Returns
    -------
    nsra_score : float
    """    

    meas_deltas = np.asarray(meas_deltas)
    pred_deltas = np.asarray(pred_deltas)
    G = meas_deltas.size

    # Measured strata
    s = np.zeros(G, dtype=np.int8)
    s[meas_deltas >  epsilon] =  1   # U
    s[meas_deltas < -epsilon] = -1   # D

    if (method == "full-pairwise"):

        # Full-pairwise (Reference definition)

        # - Explicit evaluation of all admissible gene pairs
        # - Pairwise correctness defined directly
        # - Computationally quadratic: O(G^2)        

        # Pairwise measured sign matrix (upper triangle only) 
        # sign_gt[i,j] = sign(s[i] - s[j])
        # but 0 if both are NULL
        diff_gt = s[:, None] - s[None, :]
        sign_gt = np.sign(diff_gt)

        # Mask NULL-NULL pairs
        null_null = (s[:, None] == 0) & (s[None, :] == 0)
        sign_gt[null_null] = 0

        # Keep upper triangle only
        iu = np.triu_indices(G, k=1)
        sign_gt = sign_gt[iu]

        # Predicted ordering
        diff_pred = pred_deltas[:, None] - pred_deltas[None, :]
        sign_pred = np.sign(diff_pred)[iu]

        # Valid comparisons
        valid = sign_gt != 0
        if not np.any(valid):
            return np.nan

        sg = sign_gt[valid]
        sp = sign_pred[valid]

        # Scoring
        correct = (sp == sg)
        tie = (sp == 0)

        score = np.sum(correct) + tie_score * np.sum(tie)
        return score/sg.size

    elif (method == "rank-reduced"):

        # Rank-reduced

        # - Exploits order structure induced by signed deltas
        # - Reduces pairwise consistency to rank comparisons
        # - Computationally sub-quadratic: O(G log G) time, O(G) memory        
            
        # Counts
        nU = np.sum(s ==  1)
        nN = np.sum(s ==  0)
        nD = np.sum(s == -1)

        # Total number of valid comparisons
        total_pairs = nU * (nN + nD) + nN * nD
        if total_pairs == 0:
            return np.nan

        # Sort by predicted deltas
        order = np.argsort(pred_deltas, kind="mergesort")
        s_sorted = s[order]
        p_sorted = pred_deltas[order]

        # Cumulative counters
        seen_D = 0
        seen_N = 0

        correct = 0.0

        i = 0
        while i < G:

            # Handle ties in predicted deltas
            j = i
            while j < G and p_sorted[j] == p_sorted[i]:
                j += 1

            block = s_sorted[i:j]

            nU_b = np.sum(block ==  1)
            nN_b = np.sum(block ==  0)
            nD_b = np.sum(block == -1)

            # U vs (D, N) seen so far
            correct += nU_b * seen_D
            correct += nU_b * seen_N

            # NULL vs D seen so far
            correct += nN_b * seen_D

            # ties within block
            # U vs D, U vs N, N vs D
            correct += tie_score * (nU_b * nD_b + nU_b * nN_b + nN_b * nD_b)

            # Update counters
            seen_D += nD_b
            seen_N += nN_b

            i = j

        return correct/total_pairs

    else:

        raise ValueError("Invalid method " + method)