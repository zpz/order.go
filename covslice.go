package stats

import (
    "math"
//    "github.com/skelterjohn/go.matrix"
//    "github.com/gonum/floats"
//    "fmt"
)

// We store the lower-triangular elements
// of a covariance matrix in column-major order,
// including elements on the main diagonal.
// For example, for a 3x3 cov matrix, we store these elements:
//
//     0
//     1  3
//     2  4  5
//
// We call these elements a "cov slice".
// The length of the cov slice for a k by k cov matrix is
// n = (k * k + k) / 2 = k * (k+1) / 2.
//
// Given a cov slice of length n, we can determine
// the dimensionality of the cov matrix as
// k = floor(sqrt(2 * n)).


// covslice_ndim returns the dimensionality
// of the cov matrix given the length of the cov slice.
func covslice_ndim(n int) int {
    return int(math.Floor(math.Sqrt(float64(2 * n))))
}



// covslice_expand expands a cov slice to a cov matrix (as a slice).
// The resultant slice contains elements of the cov matrix
// in either row major or column major, b/c the matrix is symmetric,
// and its row-major and col-major listing of elements are the same.
func covslice_expand(src []float64) []float64 {
    ndim := covslice_ndim(len(src))
    dst := make([]float64, ndim * ndim)
    covslice_expand_fill(src, dst)
    return dst
}




func covslice_expand_fill(src, dst []float64) {
    ndim := covslice_ndim(len(src))
    k := 0
    for col := 0; col < ndim; col++ {
        dst[col * ndim + col] = src[k]
            // Element (col, col).
        k++
        for row := col + 1; row < ndim; row++ {
            dst[row * ndim + col] = src[k]
            dst[col * ndim + row] = src[k]
                // Elements (row, col) and (col, row).
            k++
        }
    }
}

/*

// covslice_ij_idx returns the index of the element
// in a cov slice that represents the element <row, col>
// in the full cov matrix.
func covslice_ij_idx(ndim, row, col int) int {
    // Total number of elements up to, but not including,
    // column j in the cov slice:
    //        n + (n-1) +...+ (n-j+1)
    //      = (n + n-j+1)/2 * (n - (n-j+1) + 1)
    //      = j * (n + n - j + 1) / 2
    //
    // Check:
    //   j = 0: --> 0
    //   j = 1: --> n
    //   j = 2: --> n + n - 1
    //   j = 3: --> 3 * (n - 1) = n + (n-1) + (n-2)

    if row < col {
        row, col = col, row
    }

    return col * (ndim + ndim - col + 1) / 2 + row - col
}





// covslice_subset_xx_idx returns the indices of elements
// in a cov slice that would create the cov slice for
// the subset of dimensions specified in dims.
// The elements in dims are not assumed to be ordered.
func covslice_subset_xx_idx(ndim int, dims []int) []int {
    n := len(dims)
    idx := make([]int, (n*n + n)/2)
    k := 0
    for col := 0; col < n; col++ {
        for row := col; row < n; row++ {
            idx[k] = covslice_ij_idx(ndim, dims[row], dims[col])
            k++
        }
    }
    return idx
}

*/
