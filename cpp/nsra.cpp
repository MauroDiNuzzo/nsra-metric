// nsra_omp.cpp
// Reference C++ implementation of Null-stratified Rank Accuracy (NSRA)

#include <cmath>
#include <cstddef>
#include <limits>
#include <algorithm>
#include <vector>
#include <omp.h>   // OpenMP

inline double sign_with_epsilon(double x, double epsilon) {
    if (x > epsilon) return 1.0;
    if (x < -epsilon) return -1.0;
    return 0.0;
}

// Direct pointer version, avoids std::vector copies
double nsra_core_omp(const double* meas,
                     const double* pred,
                     size_t n,
                     double epsilon = 0.0,
                     double tie_score = 0.5) {

    // classify measured deltas
    std::vector<int> meas_sign(n);
    size_t count_U = 0, count_D = 0;
    for (size_t i = 0; i < n; ++i) {
        double s = sign_with_epsilon(meas[i], epsilon);
        meas_sign[i] = static_cast<int>(s);
        if (s > 0) count_U++;
        else if (s < 0) count_D++;
    }

    size_t n_active = count_U + count_D;
    if (n_active == 0) return std::numeric_limits<double>::quiet_NaN();

    // sort measured deltas descending
    std::vector<size_t> idx(n);
    for (size_t i = 0; i < n; ++i) idx[i] = i;
    std::sort(idx.begin(), idx.end(), [&meas](size_t i1, size_t i2) { return meas[i1] > meas[i2]; });

    // sort predicted deltas descending
    std::vector<size_t> pred_idx(n);
    for (size_t i = 0; i < n; ++i) pred_idx[i] = i;
    std::sort(pred_idx.begin(), pred_idx.end(), [&pred](size_t i1, size_t i2) { return pred[i1] > pred[i2]; });

    // build rank maps
    std::vector<size_t> pred_rank(n);
    for (size_t r = 0; r < n; ++r) pred_rank[pred_idx[r]] = r;

    // pairwise agreement
    double sum_score = 0.0;
    size_t count_pairs = 0;

    #pragma omp parallel for reduction(+:sum_score,count_pairs) schedule(dynamic)
    for (size_t i = 0; i < n; ++i) {
        if (meas_sign[i] == 0) continue;
        for (size_t j = i + 1; j < n; ++j) {
            if (meas_sign[j] == 0) continue;

            double gt_meas = meas_sign[i] - meas_sign[j];
            double gt_pred = 0.0;
            if (pred_rank[i] < pred_rank[j]) gt_pred = 1.0;
            else if (pred_rank[i] > pred_rank[j]) gt_pred = -1.0;
            else gt_pred = tie_score;

            if (gt_meas * gt_pred > 0) sum_score += 1.0;
            else if (gt_meas * gt_pred == 0) sum_score += tie_score;

            count_pairs++;
        }
    }

    if (count_pairs == 0) return std::numeric_limits<double>::quiet_NaN();
    return sum_score / static_cast<double>(count_pairs);
}
