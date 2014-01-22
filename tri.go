package stats

import (
	"github.com/zpz/matrix.go/dense"
)

// Lower-triangular elements
// of a square matrix in column-major order,
// including elements on the main diagonal.
// For example, for a 3x3 matrix, we store these elements:
//
//     0
//     1  3
//     2  4  5
//
// The length of the LowerTri for a k by k matrix is
// n = (k * k + k) / 2 = k * (k+1) / 2.
type LowerTri struct {
	ndim int
	data []float64
}

func NewLowerTri(k int, data []float64) *LowerTri {
	if data != nil {
		if len(data) != k*(k+1)/2 {
			panic("wrong length")
		}
		return &LowerTri{
			ndim: k,
			data: data}
	}
	return &LowerTri{
		ndim: k,
		data: make([]float64, k*(k+1)/2)}
}

// Dim returns the dimensionality
// of the matrix corresponding to a LowerTri.
func (x *LowerTri) Dim() int {
	return x.ndim
}

func (x *LowerTri) At(row, col int) float64 {
	return x.data[x.ij2idx(row, col)]
}

func (x *LowerTri) Set(row, col int, v float64) {
	x.data[x.ij2idx(row, col)] = v
}

// IJ2Idx returns the index of the element
// in x that represents the element (row, col)
// in the full matrix.
func (x *LowerTri) ij2idx(row, col int) int {
	// Total number of elements up to, but not including,
	// column j in the slice:
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

	return col*(x.ndim+x.ndim-col+1)/2 + row - col
}

func (x *LowerTri) idx2ij(idx int) (int, int) {
	panic("to be implemented")
	row, col := 0, 0
	return row, col
}

func (x *LowerTri) CopySlice(out []float64) []float64 {
	out = use_slice(out, len(x.data))
	copy(out, x.data)
	return out
}

// covslice_expand expands a cov slice to a cov matrix (as a slice).
// The resultant slice contains elements of the cov matrix
// in either row major or column major, b/c the matrix is symmetric,
// and its row-major and col-major listing of elements are the same.
func (x *LowerTri) Expand(out *dense.Dense) *dense.Dense {
	out = use_matrix(out, x.ndim, x.ndim)

	for row := 0; row < x.ndim; row++ {
		for col := 0; col < x.ndim; col++ {
			out.Set(row, col, x.At(row, col))
		}
	}
	// FIXME: some speedup may be possible if we
	// use the fact that x.data is col major and out is row major,
	// and do some slice copying.

	return out
}

// symm_fill fills a symmetric matrix with data, which contains
// its on- and below-diagonal elements in col major.
func symm_fill(data []float64, mat *dense.Dense) *dense.Dense {
	// len(data) = (n*n + n)/2 = n * (n+1) / 2
	n := int(math.Sqrt(2 * len(data)))
	mat = use_matrix(mat, n, n)
	for k, i := n, 0; k > 0; k-- {
		copy(mat.RowView(i)[i:], data[:k])
		for j := 1; j < k; j++ {
			mat.Set(n-k+j, i, data[j])
		}
		i++
		data = data[k:]
	}
	return mat
}
