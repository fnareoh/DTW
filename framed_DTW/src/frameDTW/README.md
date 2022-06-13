

# Computes the Dynamic Time Warping distances for every alignment of a pattern P and and a text T in O(#runP x #runT x k) time



In the matrix a bloc is defined by a run in P of lenght height and a run of T of length width. 

The integer k is a threshold for the DTW distance.

The border of a block is computed in O(max(k, height + width)) time.(Garance: shouldn't it be min ?)

Cf Paper: Kuszmaul https://arxiv.org/pdf/1904.09690.pdf

Cf future paper from Gourdel et al.



Computations ideas are sketched here [Note 2 févr. 2022](https://notability.com/n/M9iqdy50C4dnBZoPSJDaT)



## Usage:

`python python DTW_blocks.py P T <max_value>`

Compute the borders of the blocks of the dynamic programming matrix filled with the DTW distances. If a threshold k is fixed, the algorithm does not compute distances larger than k to save time.
Output is only for validation: for each block the 4 extreme values are shown as well as the positions where the value of the distance increases.



## Validation

`python validation.py nb_tests min_size_P_T homopol_size_bound`

Enables to validate this approach as compared to a fully computed dynamic matrix.
Compares running times.

Performs  

* nb_tests (eg 10)

* with P and T sizes at least equal to `min_size_P_T` (eg 1000)

* each letter is randomly repeated between 1 and `homopol_size_bound` time (eg 10)

**Rq:** with low homopol_size_bound (such as 1) the classical approach is faster, et vice versa. This is true for homopol_size_bound <= 5, this is due to the cost of instantiating an object block of small size.

    




