package stats

import (
	"github.com/zpz/matrix.go/dense"
)

// Mahalanobis computes squared Mahalanobis distance
// between each row of x and the vector y, with covariance matrix sigma.
func Mahalanobis(
	x *dense.Dense,
	y []float64,
	sigma *dense.Dense,
	out []float64) []float64 {

	n, p := x.Dims()
	assert(p == len(y),
		"Dimensionalities of 'x' and 'y' mismatch")
	r, c := sigma.Dims()
	assert(r == p && c == p,
		"Dimensionalities of input data and sigma mismatch")

	out = use_slice(out, n)

	// Diff matrix between x and y.
	xt := dense.NewDense(p, n)
	for row := 0; row < p; row++ {
		for col := 0; col < n; col++ {
			xt.Set(row, col, x.Get(col, row)-y[row])
		}
	}

	yy := dense.Solve(sigma, xt, nil)

	// TODO: an 'element-wise multiplication' function
	// would be better; or if dense.Dense has a 'RowView',
	// we can use Multiply row by row.
	for row := 0; row < p; row++ {
		for col := 0; col < n; col++ {
			yy.Set(row, col, yy.Get(row, col)*xt.Get(row, col))
		}
	}

	// Col sums of yy.
	for col := 0; col < n; col++ {
		out[col] = 0
		for row := 0; row < p; row++ {
			out[col] += yy.Get(row, col)
		}
	}

	return out
}
