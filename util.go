package stats

import (
	"github.com/zpz/matrix.go/dense"
)

func assert(test bool, msg string) {
	if !test {
		panic(msg)
	}
}

func use_slice(in []float64, n int) []float64 {
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

func use_matrix(in *dense.Dense, rows, cols int) *dense.Dense {
	if in == nil {
		return dense.NewDense(rows, cols)
	}
	r, c := in.Dims()
	assert(r == rows && c == cols, "Provided matrix has wrong shape")
	return in
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
