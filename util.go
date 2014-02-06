package stats

import (
	"github.com/zpz/matrix.go/dense"
)

func assert(test bool, msg string) {
	if !test {
		panic(msg)
	}
}

func use_float_slice(in []float64, n int) []float64 {
	if in == nil {
		return make([]float64, n)
	}
	return in[:n]
}

func use_int_slice(in []int, n int) []int {
	if in == nil {
		return make([]int, n)
	}
	return in[:n]
}

func use_dense(in *dense.Dense, rows, cols int) *dense.Dense {
	if in == nil {
		return dense.NewDense(rows, cols)
	}
	r, c := in.Dims()
	assert(r == rows && c == cols, "Provided matrix has wrong shape")
	return in
}

func float_clone(in []float64) []float64 {
	out := make([]float64, len(in))
	copy(out, in)
	return out
}

func int_clone(in []int) []int {
	out := make([]int, len(in))
	copy(out, in)
	return out
}

func smaller(a, b int) int {
	if a < b {
		return a
	}
	return b
}

func bigger(a, b int) int {
	if a > b {
		return a
	}
	return b
}

// exclude returns integers between 0 (inclusive) and
// n (exclusive), excluding those in index.
// Does not assume elements in include are sorted
// (otherwise the code can be more efficient).
func exclude(n int, index []int) []int {
	idx := make([]int, n)
	for _, i := range index {
		idx[i] = 1
	}
	k := 0
	for j := 0; j < n; j++ {
		if idx[j] == 0 {
			idx[k] = j
			k++
		}
	}
	return idx[:k]
}

// pick_floats returns a subslice with the elements
// at the specified indices.
func pick_floats(x []float64, index []int, y []float64) []float64 {
	y = use_float_slice(y, len(index))
	for i, j := range index {
		y[i] = x[j]
	}
	return y
}

// xtAy computes the scalar value that is
//  x^t * A * y
// where 'x' and 'y' are treated as column vectors,
// and '*' is matrix multiplication.
func xtAy(x []float64, A *dense.Dense, y []float64) float64 {
	v := 0.0
	k := len(x)
	for i := 0; i < k; i++ {
		v += FloatDot(A.RowView(i), y) * x[i]
	}
	return v
}
