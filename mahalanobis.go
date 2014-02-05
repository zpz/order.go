package stats

import (
	"github.com/zpz/matrix.go/dense"
)

func xtAx(x []float64, A *dense.Dense) float64 {
	v := 0.0
	k := len(x)
	for i := 0; i < k; i++ {
		v += FloatDot(A.RowView(i), x) * x[i]
	}
	return v
}

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

	if chol, ok := dense.Chol(sigma); ok {
		sigma = chol.Inv(nil)
	} else {
		panic("Cholesky failed on cov matrix")
	}

	out = use_float_slice(out, n)

	diff := make([]float64, p)
	for i := 0; i < n; i++ {
		FloatSubtract(x.RowView(i), y, diff)
		out[i] = xtAx(diff, sigma)
	}

	return out
}
