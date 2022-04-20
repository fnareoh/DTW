

# Computes a Dynamic Time Warp distance between two sequences Q and T in O(#runQ x #runT x k) time



In the matrix a bloc is defined by a a^height x b^width letters to be compared

k bounds the computations. This is a maximal value.

Computations of a specific block is in O(max(k, height + width)).(Garance: shouldn't it be min ?)

Cf Paper: Kuszmaul https://arxiv.org/pdf/1904.09690.pdf

Cf future paper from Gourdel et al.



Computations ideas are sketched here [Note 2 févr. 2022](https://notability.com/n/M9iqdy50C4dnBZoPSJDaT)



## Usage:

`python python DTW_blocks.py Q T <max_value>`

Compute the framed DTW of sequences Q versus T. If a max value is fixed, the maximal value in the matrix does not get higher than this value. This saves times
Output is only for validation : for each block the 4 extreme values are shown as well as *cuts* positions (that is to say where value increase along the frames).



## Validation

`python validation.py nb_tests min_size_Q_T homopol_size_bound`

Enables to validate this approach as compared to a fully computed dynamic matrix.
Compare running times.

Performs  

* nb_tests (eg 10)

* with Q and T sizes at least equal to `min_size_Q_T` (eg 1000)

* each letter is randomly repeated between 1 and `homopol_size_bound` time (eg 10)

**Rq:** with low homopol_size_bound (such as 1) the classical approach is faster, et vice versa. This is true for homopol_size_bound <= 5, this is due to the cost of instantiating an object block of small size.

    




