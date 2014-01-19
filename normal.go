package stats

import (
	"github.com/zpz/matrix.go/dense"
	"math"
	"math/rand"
)

// Normal is a struct for a normal distribution.
type Normal struct {
	mean float64
	sd   float64
}

// NewNormal creates a Normal object and returns a pointer to it.
func NewNormal(mean, sd float64) *Normal {
	assert(sd > 0, "sd must be positive")
	return &Normal{mean: mean, sd: sd}
}

// Density computes the probability density of the given
// values in the normal distribution norm.
// out either has the correct length or is nil,
// in which case a slice is allocated and used.
func (norm *Normal) Density(
	x []float64,
	out []float64) []float64 {

	coef := 1.0 / (Sqrt2Pi * norm.sd)

	out = use_slice(out, len(x))

	for i, v := range x {
		z := (v - norm.mean) / norm.sd
		out[i] = coef * math.Exp(-0.5*z*z)
	}

	return out
}

// Random generates n random numbers from the normal distribution norm.
func (norm *Normal) Random(
	n int,
	_ RandomNumberGenerator,
	out []float64) []float64 {

	out = use_slice(out, n)
	for i := range out {
		out[i] = rand.NormFloat64()*norm.sd + norm.mean
	}
	return out
}

// Mvnormal is a data structure for multivariate normal distribution.
type Mvnormal struct {
	mean        []float64
	cov         *dense.Dense
	cov_det     float64
	log_cov_det float64
}

// NewMvnormal creates a new Mvnormal object and returns a pointer to
// it. Note that the input mean slice and cov matrix are not deeply
// copied. They are simply pointed to by the created Mvnormal object.
func NewMvnormal(mean []float64, cov *dense.Dense) *Mvnormal {
	r, c := cov.Dims()
	assert(r == c && c == len(mean),
		"Dimensionalities of 'mean' and 'cov' mismatch")
	cov_det := cov.Det()
	return &Mvnormal{
		mean:        mean,
		cov:         cov,
		cov_det:     cov_det,
		log_cov_det: math.Log(cov_det)}
}

// Dim returns the dimensionality of the Mvnormal distribution mvn.
func (mvn *Mvnormal) Dim() int {
	return len(mvn.mean)
}

// Density computes the prob density values for each row of x (as a
// case, or observation) in the Mvnormal distribution mvn.
// out either has the correct length (i.e. the number of rows in x)
// or is nil, in which case a slice is allocated and used.
func (mvn *Mvnormal) Density(
	x *dense.Dense,
	// Each row is a case.
	out []float64) []float64 {

	n, p := x.Dims()
	assert(len(mvn.mean) == p,
		"Dimensionalities of input 'x' and distribution mismatch")

	out = use_slice(out, n)

	// Mahalanobis distance
	Mahalanobis(x, mvn.mean, mvn.cov, out)

	coef := 1.0 / math.Sqrt(
		math.Pow(2*math.Pi, float64(p))*mvn.cov_det)

	for i, v := range out {
		out[i] = coef * math.Exp(-0.5*v)
	}

	return out
}

// Random generates random samples from the multivariate normal
// distribution mvn.
// The result is a matrix with each row being a case.
func (mvn *Mvnormal) Random(
	n int,
	_ RandomNumberGenerator,
	out *dense.Dense) *dense.Dense {

	p := len(mvn.mean)

	U, ok := dense.CholeskyR(mvn.cov)
	assert(ok, "Cholesky failed")

	z := dense.DenseView(NewNormal(0., 1.).Random(p*n, &RNG{}, nil), n, p)

	out = dense.Mult(z, U, out)

	for row := 0; row < n; row++ {
		Add(out.RowView(row), mvn.mean, out.RowView(row))
	}

	// Alternative algorithm: see function 'mvrnom' in the R package
	// MASS.

	return out
}
