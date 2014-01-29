package stats

import (
	"github.com/zpz/matrix.go/dense"
)

// Mahalanobis computes _squared_ Mahalanobis distance
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

	out = use_float_slice(out, n)

	// Diff matrix between x and y.
	// Each col is an observation; each row is a dimension.
	// xt has dimensions p, n
	xt := dense.T(x, nil)
	for row := 0; row < p; row++ {
		FloatShift(xt.RowView(row), -y[row], xt.RowView(row))
	}

	// yy has dimensions p, n
	if chol, ok := dense.Chol(sigma); ok {
		yy := chol.Solve(dense.Clone(xt))

		// Element-wise multiplication.
		yy.Elemult(xt)

		// Col sums of yy.
		yy.GetRow(0, out)
		for row := 1; row < p; row++ {
			FloatAdd(out, yy.RowView(row), out)
		}

		return out
	}

	panic("Cholesky failed on cov matrix")
}
