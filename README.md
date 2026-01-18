
# Null-Stratified Rank Accuracy (NSRA): a rank-based metric for evaluating Perturb-seq predictions

![status](https://img.shields.io/badge/status-proposal--stage-yellow)


## Metric overview

Null-Stratified Rank Accuracy (NSRA) is a rank-based evaluation metric for perturbation-induced gene expression prediction, designed for sparse and heterogeneous Perturb-seq data.

NSRA measures pairwise order consistency between predicted and measured signed gene expression changes relative to a control state, under a partial order in which non-differentially expressed (non-responding) genes form an unconstrained equivalence class. Specifically, the metric enforces the biologically meaningful ordering
up-regulated > non-changing > down-regulated genes, while ignoring the internal ordering of non-changing genes.

By collapsing the noise-dominated null region into a single equivalence class, NSRA avoids spurious comparisons among non-responding genes and yields a threshold-free, magnitude-agnostic score that remains sensitive to sign errors and false positives. The metric is intended to complement magnitude-based scores (e.g., MAE) by assessing whether a model captures the structural organization of the transcriptional response induced by a perturbation.

### DISCLAIMER

> The development of evaluation metrics for Perturb-seq is an open, community-driven effort, and NSRA is proposed in this spirit. This work is not linked to any academic research project, consortium, or private for-profit initiative, and it has not undergone peer review. At its current stage, NSRA should be regarded as a proposed metric, rather than a validated or established standard.
>
> NSRA has not yet been systematically tested on large, diverse real-world Perturb-seq datasets, nor has it been comprehensively benchmarked against existing magnitude-based or rank-based metrics. While the formulation is motivated by biological and statistical considerations, it is possible that the metric exhibits unforeseen limitations, edge cases, or failure modes in specific experimental regimes. Accordingly, its suitability for any particular application is not warranted at this stage.
>
> The primary goal of releasing NSRA at this point is to make the idea and implementation available to the community for scrutiny, experimentation, and discussion. Feedback, empirical testing, critical assessment, and contributions (whether conceptual, methodological, or implementation-related) are very welcome and strongly encouraged. Community input will be essential to determine whether NSRA is useful in practice, how it should be refined, or whether alternative formulations are preferable.

## Background and motivation

Gene expression responses to perturbations are high-dimensional, sparse, and heterogeneous across experiments. In the context of perturbation-based gene expression prediction, no single scalar metric can fully capture model performance across the diverse regimes encountered in real data. Magnitude-based scores such as mean absolute error (MAE) or weighted variants are essential for assessing numerical accuracy, but they necessarily conflate multiple effects, including scale calibration, noise modeling, and the relative importance assigned (explicitly or implicitly) to different genes. As a result, magnitude-based metrics can be dominated by non-changing genes and may obscure biologically meaningful structural errors, such as incorrect directionality or spurious activation of irrelevant genes.

As a complement to these measures, it is useful to consider a purely rank-based metric that evaluates whether a model correctly captures the *relative biological ordering* induced by a perturbation, independent of absolute effect sizes. Such a metric focuses on the structure of the response, i.e., which genes are more up- or down-regulated relative to others, while avoiding thresholds, weights, or penalties that may encode subjective assumptions about biological relevance.

Rank-based evaluation is not intended to replace magnitude-sensitive metrics. Indeed, magnitude-based and rank-based metrics capture orthogonal aspects of performance: the former assesses quantitative fidelity, while the latter evaluates structural and directional correctness. Models that perform well on both are therefore accurate in magnitude and faithful to the underlying biological response, whereas discrepancies between the two can reveal specific failure modes. This complementary evaluation strategy supports fair benchmarking, robust leaderboard comparisons, and more transparent interpretation of model behavior.

Purely rank-based metrics, however, have important drawbacks in the context of Perturb-seq data. First, perturbations are often sparse: many interventions induce few or no reproducible transcriptional changes, and most genes lie in a long tail dominated by noise. In such regimes, global rankings over all genes become ill-defined, as the relative order of non-changing genes is effectively arbitrary and should not contribute to model evaluation. Standard rank-based scores, such as mean rank error or Kendall-type correlations computed over full rankings, are therefore overly sensitive to noise in this null region and can behave unpredictably when the true signal is weak. Second, naive ranking schemes struggle to penalize false positives: a model that correctly ranks a small set of truly responsive genes but assigns spurious structure to hundreds of non-changing genes may achieve a deceptively good rank-based score, despite poor biological specificity. Finally, by construction, rank-based metrics discard magnitude information and thus cannot distinguish between large and marginal effects; they must therefore be interpreted strictly as measures of structural consistency rather than practical utility for downstream decision-making.

To address these limitations while retaining a threshold-free, magnitude-agnostic evaluation, we introduce a specific formulation of rank-based accuracy. The key idea is to replace the full ordering over all genes with a biologically motivated partial order: truly up-regulated genes should rank above non-changing genes, which in turn should rank above truly down-regulated genes, while the internal ordering of non-changing genes is treated as unconstrained. This effectively collapses the null region into a single equivalence class, removing noise-dominated comparisons while preserving meaningful ordering constraints. Model performance is then assessed by counting pairwise ordering agreements only among informative gene pairs, yielding a rank-based score that remains sensitive to false positives and sign errors, yet robust to sparse perturbations and long-tailed noise. This formulation avoids thresholds, weights, or penalties in the evaluation itself, and instead focuses on satisfying biologically meaningful ordering constraints.

Conceptually, this idea is related to existing approaches in information retrieval (e.g., partially ordered relevance judgments), survival analysis with tied ranks, and constrained or masked variants of Kendall’s tau. However, explicitly treating non-differentially expressed genes as an unconstrained equivalence class within a signed, perturbation-centric ranking is not standard practice in Perturb-seq benchmarking and is rarely made explicit in gene expression prediction challenges. In this sense, NSRA can be viewed as a principled adaptation of existing rank-based concepts to the statistical and biological structure of perturbational transcriptomics.


## Formalization

Model predictions are evaluated against the measured perturbation response, defined as signed gene expression changes relative to a control state. Genes are stratified into upregulated, downregulated, and non-differentially expressed sets based on the measured data. Evaluation is performed using a null-stratified rank-based metric that assesses consistency with this measured partial order.

The objective is to evaluate not only predictive correctness, but also the interpretability and structural fidelity of the predicted transcriptional response.

### Setup and notation

Let:

* $`G = \{g_1,\dots,g_N\}`$ denote the set of genes,
* $`x^{(0)} \in \mathbb{R}^N`$ the reference (control) expression vector,
* $`x^{(\mathrm{meas})} \in \mathbb{R}^N`$ the measured perturbed expression,
* $`x^{(\mathrm{pred})} \in \mathbb{R}^N`$ a model prediction.

For each gene $`g \in G`$, we define measured and predicted signed expression deltas relative to control as:

```math
\Delta^{(\mathrm{meas})}_g = x^{(\mathrm{meas})}_g - x^{(0)}_g, \qquad \Delta^{(\mathrm{pred})}_g = x^{(\mathrm{pred})}_g - x^{(0)}_g,
```

when computed in linear space, or

```math
\Delta^{(\mathrm{meas})}_g = \log(x^{(\mathrm{meas})}_g + 1) - \log(x^{(0)}_g + 1), \qquad
\Delta^{(\mathrm{pred})}_g = \log(x^{(\mathrm{pred})}_g + 1) - \log(x^{(0)}_g + 1),
```

when computed in log-transformed space.

Because NSRA depends only on ordering, the choice of delta definition matters only insofar as it affects monotonicity across genes. In practice, up- and down-regulated gene sets are assumed to be defined using signed log1p expression deltas relative to control, ensuring scale-invariant detection of regulatory effects. However, because log transformation compresses effect sizes, this choice influences which genes enter or exit the null set (see below).

Throughout, all rankings are assumed to be descending in $`\Delta_g`$: most upregulated genes first, most downregulated genes last.

### Gene stratification and partial order

Using the full experimental dataset (including replicates), genes are partitioned once and for all into three disjoint sets:

* $`U`$: significantly upregulated genes,
* $`D`$: significantly downregulated genes,
* $`N = G \setminus (U \cup D)`$: non-differentially expressed (null) genes.

This partition is external to the metric: it is not model-dependent and is fixed across all prediction evaluations.

The biological ordering implied by the data is a partial order, not a total order. Between classes, the following relations hold:

```math
\forall u \in U, \forall n \in N, \forall d \in D:\quad
u \succ n \succ d.
```

Within classes, ordering is defined as follows:

* For $`u,u\prime \in U`$: ordered by $`\Delta^{(\mathrm{meas})}`$,
* For $`d,d\prime \in D`$: ordered by $`\Delta^{(\mathrm{meas})}`$,
* For $`n,n\prime \in N`$: no ordering constraint.

Thus, $`N`$ forms a single equivalence class with unconstrained internal ordering.

For convenience, we define a class label function:

```math
c(g) =
\begin{cases}
+1, & g \in U, \\
0, & g \in N, \\
-1, & g \in D.
\end{cases}
```

### Informative pair set

Evaluation is based on pairwise comparisons that are meaningful under the partial order. We define the set of informative ordered gene pairs as:

```math
\mathcal{P} =
(U \times U)
\cup
(D \times D)
\cup
(U \times N)
\cup
(N \times D)
\cup
(U \times D),
```

corresponding respectively to:

* up-up comparisons,
* down-down comparisons,
* up vs null,
* null vs down,
* up vs down.

Equivalently, this set can be written compactly as:

```math
\mathcal{P} = {(g,h) \in G \times G : g \neq h, \neg(c(g)=0 \wedge c(h)=0)}.
```

By construction, pairs in $`N \times N`$ are excluded.

Because genes in $`N`$ form an unconstrained equivalence class, there is no unique total ranking consistent with the biological structure. Assigning explicit ranks to null genes would therefore introduce arbitrary ordering artifacts. For this reason, NSRA is formulated in terms of pairwise signs rather than explicit ranks. Importantly, this does not lose information: explicit ranks encode global positions, while signs of pairwise differences encode the same ordering information locally. The pairwise formulation allows biologically meaningful comparisons to be counted while excluding null–null comparisons entirely.

### Ordering sign

For any pair $`(g,h) \in \mathcal{P}`$, we define the true ordering:

```math
\mathrm{sign}_{\mathrm{meas}}(g,h) =
\begin{cases}
\mathop{\text{sign}}\left(\Delta^{(\mathrm{meas})}_g - \Delta^{(\mathrm{meas})}_h\right),
& \text{if } c(g)=c(h)\in{+1,-1} \\
\mathop{\text{sign}}\left(c(g)-c(h)\right),
& \text{if } c(g)\neq c(h)
\end{cases}
```

This makes it explicit that class separation dominates, with deltas only refining ordering within classes.

By definition, $`\mathrm{sign}_{\mathrm{meas}}(g,h)`$ encodes the biologically correct ordering sign for the ordered pair $`(g,h)`$:

* $`+1`$ means that gene $`g`$ should be ranked above gene $`h`$ (i.e., $`g \succ h`$, or $`\Delta_g > \Delta_h`$)

* $`−1`$ means that gene $`g`$ should be ranked below gene $`h`$ (i.e., $`g \prec h`$, or $`\Delta_g < \Delta_h`$)

* $`0`$ means no constraint (or tie, see below)

Thus, the definition is antisymmetric:

```math
\mathrm{sign}_{\mathrm{meas}}(g,h) = -\mathrm{sign}_{\mathrm{meas}}(h,g)
```

It is noted that $`N`$ is an equivalence class internally, but it is ordered relative to $`U`$ and $`D$. Cross-class ordering uses only class membership, while within-class ordering uses signed deltas.

For clarity, here is the full cross-class behavior:

| $`g`$ | $`h`$ | $`c(g)`$ | $`c(h)`$ | $`\mathrm{sign}_{\mathrm{meas}}(g,h)`$ | Interpretation |
| --- | --- | ------ | ------ | -------------- | -------------- |
| $`U`$   | $`N`$   | $`+1`$     | $`0`$      | $`+1`$             | up > null      |
| $`N`$   | $`U`$   | $`0`$      | $`+1`$     | $`−1`$             | null < up      |
| $`U`$   | $`D`$   | $`+1`$     | $`−1`$     | $`+1`$             | up > down      |
| $`D`$   | $`U`$   | $`−1`$     | $`+1`$     | $`−1`$             | down < up      |
| $`N`$   | $`D`$   | $`0`$      | $`−1`$     | $`+1`$             | null > down    |
| $`D`$   | $`N`$   | $`−1`$     | $`0`$      | $`−1`$             | down < null    |
| $`U`$   | $`U`$   | $`+1`$     | $`+1`$     | $`\mathop{\text{sign}}(\Delta_g−\Delta_h)`$    | order within $`U`$ |
| $`D`$   | $`D`$   | $`−1`$     | $`−1`$     | $`\mathop{\text{sign}}(\Delta_g−\Delta_h)`$    | order within $`D`$ |

The predicted ordering uses only the model output (no class information enters here), which is simply:

```math
\mathrm{sign}_{\mathrm{pred}}(g,h) =
\mathop{\text{sign}}(\Delta^{(\mathrm{pred})}_g - \Delta^{(\mathrm{pred})}_h).
```

### Pairwise correctness indicator

Now, we can define correctness as:

```math
C(g,h) =
\begin{cases}
1, & \mathrm{sign}_{\mathrm{pred}}(g,h) = \mathrm{sign}_{\mathrm{meas}}(g,h) \\
0, & \text{otherwise}
\end{cases}
```

or, more compactly:

```math
C(g,h) = 
\mathbb{I}\left[
\mathrm{sign}_{\mathrm{pred}}(g,h)
=
\mathrm{sign}_{\mathrm{meas}}(g,h)
\right]
```

where $`\mathbb{I}[\cdot]`$ is the indicator function.

Ties in the pairwise comparison can arise from two sources, as follows. 

First, we can have measured within-class ties, i.e., genes in $`U`$ or $`D`$ may have identical or nearly identical $`\Delta^{(\mathrm{meas})}`$, especially after averaging or pseudo-bulk aggregation. These are true ties, and no ordering is biologically meaningful within the tied set. Formally, for a pair $`(g,h)\in U\times U`$ or $`D\times D`$ we have:

```math
\text{if } \Delta^{(\mathrm{meas})}_g = \Delta^{(\mathrm{meas})}_h \implies \mathrm{sign}_{\mathrm{meas}}(g,h) = 0
```

which means no constraint (any predicted ordering is acceptable).

Second, a model may predict identical $`\Delta^{(\mathrm{pred})}`$ for multiple genes (common in sparse models or rounding), which should be penalized only if it contradicts the biologically meaningful ordering. Formally, for a pair $`(g,h)`$ we have:

```math
\mathrm{sign}_{\mathrm{pred}}(g,h) =
\begin{cases}
+1 & \Delta^{(\mathrm{pred})}_g > \Delta^{(\mathrm{pred})}_h \\
-1 & \Delta^{(\mathrm{pred})}_g < \Delta^{(\mathrm{pred})}_h \\
0 & \Delta^{(\mathrm{pred})}_g = \Delta^{(\mathrm{pred})}_h
\end{cases}
```

If a model predicts a tie where a strict ordering exists in the measured data, it is partially correct (it does not reverse the order, but fails to distinguish the two genes). Thus, we can assign this case a partial credit, which is a standard approach in rank correlation metrics (e.g., Kendall's $`\tau_b`$).

In particular, with ties the correctness indicator $`C(g,h)`$ can be written as:

```math
C(g,h) =
\begin{cases}
1, & \mathrm{sign}_{\mathrm{meas}}(g,h) = 0 \quad \text{(measured tie)} \\
1, & \mathrm{sign}_{\mathrm{meas}}(g,h) \neq 0 \land \mathrm{sign}_{\mathrm{pred}}(g,h) = \mathrm{sign}_{\mathrm{meas}}(g,h) \\
0.5, & \mathrm{sign}_{\mathrm{meas}}(g,h) \neq 0 \land \mathrm{sign}_{\mathrm{pred}}(g,h) = 0 \\
0, & \text{otherwise}
\end{cases}
```

We also incorporate soft sign-tolerance to reduce sign instability near zero. Indeed, in real expression data, tiny differences may be numerically non-zero but biologically irrelevant. That could be taken into account by introducing a small tolerance $`\epsilon > 0`$:

```math
\Delta^{(\mathrm{meas})}_g - \Delta^{(\mathrm{meas})}_h \in [-\epsilon, \epsilon] \implies \mathrm{sign}_{\mathrm{meas}}(g,h) = 0
```

and similarly for predictions:

```math
\Delta^{(\mathrm{pred})}_g - \Delta^{(\mathrm{pred})}_h \in [-\epsilon, \epsilon] \implies \mathrm{sign}_{\mathrm{pred}}(g,h) = 0
```

The choice of $`\epsilon`$ might depends on measurement noise level, pseudo-bulk averaging, and/or normalization scale. In particular, $`\epsilon`$ should approximate the scale of irreducible noise (i.e., noise floor) in $`\Delta`$ calculated in the relevant (linear or log) space. Notably, log-space deltas admit a globally meaningful $`\epsilon`$, while linear-space deltas do not unless made relative.

### Score definition

At this point, we can finally define the metric as the average correctness over all informative pairs:

```math
\boxed{
\mathrm{NSRA} =
\frac{1}{|\mathcal{P}|}
\sum_{(g,h) \in \mathcal{P}} C_{\epsilon}(g,h) \quad \in [0,1]
}
```

where $`C_{\epsilon}(g,h)`$ is computed using tolerance-aware signs as above.

This is somewhat equivalent to a Kendall-style pairwise accuracy computed on a restricted pair set. Accordingly, NSRA could be expressed as $`\tau_{\text{NSRA}} = 2 \cdot \mathrm{NSRA} - 1 \in [-1,1]`$ if a signed correlation-like score is preferred.

NSRA is to be interpreted as a biological ordering consistency metric representing a structural fidelity score, with the following outcomes: 

* biologically faithful (near perfect) ordering: $`\approx 1.0`$
* largely correct structure with local errors: $`\approx 0.7-0.9`$
* random ordering (no usable ordering information): $`\approx 0.5`$
* systematically misleading predictions: $`<0.5`$

These interpretations are independent of scale, normalization, and noise variance, i.e., precisely the regime where magnitude-based metrics become ambiguous.

The metric has the followed properties:

* fully rank-based
* no weights
* no penalties
* no thresholds inside the metric (DEG partition is fixed externally)
* robustness to numerical noise via $`\epsilon`$

Characteristic behaviors include:

* partial credit for ties to avoiding over-penalization
* natural down-weighting of sparse perturbations (i.e., few informative pairs)
* strong penalization of hallucinated signal, random ordering, and sign flips
* symmetry between false positives and false negatives (penalty linear in the number of missed constraints)
* undefined score when all genes are inside null (no biological signal to evaluate)

Importantly, the metric is not a measure of effect-size accuracy and should not be used in isolation for decision-making. It is, however, well-suited for model comparison, benchmarking, and leaderboards where structural biological fidelity is the primary concern.

## Implementation

### Computational challenge

The naive evaluation of NSRA via explicit pairwise comparisons scales as
$`O(G^2)`$,
where $`G`$ is the number of genes. For $`G \approx 20{,}000`$, the number of unordered pairs is:

```math
\frac{G(G-1)}{2} \approx 2 \times 10^8.
```

Materializing even a single dense array of this size is prohibitive: a `float64` matrix requires $`\sim 1.6`$ GB of memory, while an `int8` representation still requires $`\sim 200`$ MB. Since multiple intermediate arrays would be needed, peak memory usage quickly exceeds what is practical on typical machines, and the total number of operations leads to runtimes on the order of tens of seconds to minutes.

As a result, a dense pairwise implementation of NSRA is neither memory-efficient nor scalable to genome-scale data.

### Aggregated formulation

Although NSRA is defined in terms of pairwise comparisons, it does not require explicit enumeration of all gene pairs. In particular, pairs within the null stratum $`(N \times N)`$ are excluded by construction, and in typical Perturb-seq experiments the number of truly responsive genes is small:

```math
|U| + |D| \ll G.
```

All informative comparisons involve at least one gene in $`U`$ or $`D`$. Exploiting this structure, NSRA can be rewritten in aggregated form as:

```math
\mathrm{NSRA} =
\frac{
\sum_{g \in U}\sum_{h \in N \cup D} \mathbb{I}\left[\Delta^{(\mathrm{pred})}_g > \Delta^{(\mathrm{pred})}_h\right]
+
\sum_{g \in N}\sum_{h \in D} \mathbb{I}\left[\Delta^{(\mathrm{pred})}_g > \Delta^{(\mathrm{pred})}_h\right]
}{|\mathcal{P}|},
```

where $`|\mathcal{P}|`$ is the number of informative ordered pairs.

This formulation is **exactly equivalent** to the dense pairwise sign test defined in the Formalization section, but avoids constructing pairwise matrices altogether. It relies only on the relative ordering of predicted deltas and the measured gene stratification.

### Algorithmic implementation

NSRA is computed using a rank-based aggregation procedure with the following steps:

1. **Measured stratification**
   Compute the measured class label
   $`s_{\mathrm{meas}}(g) \in {+1, 0, -1}`$
   using signed deltas and tolerance $`\epsilon`$, assigning genes to $`U`$, $`N`$, or $`D`$.

2. **Sorting by predictions**
   Sort all genes by predicted deltas $`\Delta^{(\mathrm{pred})}`$ in descending order.

3. **Single pass aggregation**
   Traverse genes in predicted order while maintaining counters for how many $`N`$ and $`D`$ genes have already been encountered.

4. **Accumulation of correct comparisons**

   * When encountering a gene in $`U`$, count all correctly ordered $`U-N`$ and $`U-D`$ pairs.
   * When encountering a gene in $`N`$, count correctly ordered $`N-D`$ pairs.
   * Pairs within the same class $`(U-U), (D-D), (N-N)`$ are ignored, consistent with the partial order definition.

5. **Handling predicted ties**
   Predicted ties are handled explicitly by assigning fractional credit (typically 0.5), consistent with the pairwise definition.

This procedure has time complexity $`O(G \log G)`$ due to sorting, and memory complexity $`O(G)`$, making it suitable for genome-scale Perturb-seq data.

### Practical implementations

We provide implementations in:

* **Python (NumPy)**: minimal dependencies, clarity-first implementation.
* **MATLAB**: benefits from JIT compilation and efficient vectorized operations.
* **C++ backend** with:

  * MATLAB (`mex`) wrapper
  * Python (`pybind11`) wrapper

While all implementations use the same aggregated algorithm, performance differs across environments. In practice:

* MATLAB is fastest for the pure interpreted implementation due to JIT optimization.
* The compiled Python (`pybind11`) backend is faster than the MATLAB MEX version, which incurs additional function-call and memory-access overhead.

The C++ backend is optional, and the pure Python implementation remains feasible for genome-scale data. The Python implementation requires only NumPy, no additional dependencies are needed.


## FAQs
#### Why not use a single metric?
Gene expression responses to perturbations are sparse, noisy, and heterogeneous. A single scalar metric inevitably collapses multiple, qualitatively different failure modes into one number. Magnitude-based metrics and rank-based metrics probe orthogonal aspects of model behavior; using both is useful to distinguish numerical accuracy from structural biological correctness.
#### Why isn't MAE sufficient?
MAE measures numerical agreement but is dominated by the large set of non-changing genes. In sparse perturbations, a model can achieve a low MAE while at the same time missing most true responders, predicting the wrong direction of change, or hallucinating structured responses in irrelevant genes. MAE is therefore necessary but not sufficient to assess biological fidelity.
#### Why not use correlation (Pearson or Spearman)?
Correlation metrics measure linear or monotonic association but are insensitive to absolute errors and can be unstable when true signal is weak. In particular, Pearson correlation conflates scale calibration and directionality, and Spearman correlation implicitly ranks all genes, including those dominated by noise, which reintroduces the same problems as naive rank-based scores. Correlation is informative diagnostically but ill-suited as a primary leaderboard metric.
#### Why use rank-based evaluation at all?
Rank-based evaluation focuses on whether a model captures the relative biological structure of a perturbation, i.e., which genes are more up- or down-regulated than others, independent of effect size calibration. This is particularly valuable when effect sizes vary widely across perturbations, absolute magnitudes are difficult to compare across conditions, or the downstream use case depends on prioritization rather than precise quantification.
#### Why not rank all genes directly?
In Perturb-seq data, most genes do not respond to a given perturbation. Their relative ordering is largely driven by noise and is biologically meaningless. Full rankings therefore introduce arbitrary constraints and can overwhelm the contribution of truly responsive genes. Any rank-based metric that treats all gene pairs equally will be unstable in sparse regimes.
#### What problem does the null-awareness solve?
The NSRA metric explicitly acknowledges that non-differentially expressed genes form a null region with no meaningful internal ordering. By collapsing these genes into a single equivalence class, noise-dominated comparisons are excluded, sparse perturbations are handled robustly, and evaluation focuses on biologically meaningful ordering constraints (up > null > down). This preserves the advantages of rank-based evaluation without forcing arbitrary decisions in the null region.
#### How are false positives penalized without thresholds?
False positives are penalized through pairwise ordering violations. If a gene that is truly non-changing is predicted as strongly up- or down-regulated, it will incorrectly outrank true responders or be incorrectly outranked by them, contributing negatively to the score. No explicit thresholds or weights are required in the metric itself.
#### Why use pairwise comparisons instead of explicit ranks?
Pairwise comparisons encode ordering constraints directly and naturally support partial orders. Explicit ranks require assigning a total order, which is inappropriate when large subsets of items are intentionally unordered. Pairwise formulations are standard in ordinal and preference-based evaluation and allow the null region to be treated as unconstrained.
#### Why are signed deltas used instead of absolute expression levels?
Perturb-seq evaluation is inherently comparative: the biological signal of interest is the change induced by a perturbation relative to a control state. Signed deltas encode directionality explicitly, allow up- and down-regulation to be treated symmetrically, and are consistent across bulk and pseudo-bulk representations.
#### Should deltas be absolute or relative?
Both are admissible, but they emphasize different aspects. Absolute deltas favor genes with large baseline expression and reflect raw transcriptional impact. Relative deltas normalize for baseline abundance and emphasize proportional change. The metric definition is agnostic to this, and the choice should be driven by biological intent and consistency with the ground-truth processing.
#### Is this metric meant to replace other scores?
No. The NSRA metric is explicitly complementary to magnitude-based metrics. It is designed to reveal structural and directional errors that are otherwise hidden, not to summarize overall model quality on its own.
#### Is this approach entirely new?
The underlying ideas (pairwise ranking, partial orders, equivalence classes) are well established in statistics and information retrieval. However, their explicit combination into a signed, perturbation-centric rank metric with a null region is not standard practice in Perturb-seq benchmarking and represents a principled adaptation to the structure of the data.


## How to cite

If you use NSRA in academic work, we recommend citing this repository and, where appropriate, describing the metric explicitly in the methods section.

#### Suggested citation (software):
```
DiNuzzo M. Null-Stratified Rank Accuracy (NSRA): a rank-based metric for evaluating Perturb-seq predictions. 
GitHub repository: https://github.com/MauroDiNuzzo/nsra-metric
```

A machine-readable citation is provided via GitHub's citation feature (CITATION.cff). For LaTeX users, a BibTeX entry is also provided below.

```bibtex
@software{nsra_metric,
  author = {DiNuzzo, Mauro},
  title = {Null-Stratified Rank Accuracy (NSRA)},
  year = {2026},
  url = {https://github.com/MauroDiNuzzo/nsra-metric},
  note = {Rank-based metric for evaluating Perturb-seq predictions}
}
```

#### Suggested methods description:

*Model performance was evaluated using Null-Stratified Rank Accuracy (NSRA), a rank-based metric that assesses pairwise order consistency of predicted perturbation-induced gene expression changes under a partial order in which non-differentially expressed genes form an unconstrained equivalence class.*
