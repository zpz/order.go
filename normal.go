package stats

import (
    "math"
    "math/rand"
    "github.com/gonum/matrix/mat64"
)



type Normal struct {
    mean float64
    sd float64
}



func NewNormal(mean, sd float64) *Normal {
    assert(sd > 0, "sd must be positive")
    return &Normal{mean: mean, sd: sd}
}





// DNorm computes the probability density of the given
// values in a normal distribution.
func (norm *Normal) Density(
    x []float64,
    out []float64) []float64 {

    coef := 1.0 / (Sqrt2Pi * norm.sd)

    out = use_slice(out, len(x))

    for i, v := range x {
        z := (v - norm.mean) / norm.sd
        out[i] = coef * math.Exp(-0.5 * z * z)
    }

    return out
}





// Random generates n random numbers from the normal distribution.
func (norm *Normal) Random(
    n int,
    _ RandomNumberGenerator,
    out []float64) []float64 {

    return Generate(n,
        func(int) float64 {
            return rand.NormFloat64() * norm.sd + norm.mean },
        out)
}




type Mvnormal struct {
    mean []float64
    cov *mat64.Dense
    cov_det float64
    log_cov_det float64
}




func NewMvnormal(mean []float64, cov *mat64.Dense) *Mvnormal {
    r, c := cov.Dims()
    assert(r == c && c == len(mean),
        "Dimensionalities of 'mean' and 'cov' mismatch")
    cov_det := mat64.Det(cov)
    return &Mvnormal{
        mean: mean,
        cov: cov,
        cov_det: cov_det,
        log_cov_det: math.Log(cov_det)}
}




func (mvn *Mvnormal) Dim() int {
    return len(mvn.mean)
}



// Density of a multivariate normal distribution.
// Input x contains one or more cases;
// output is densities for each case in x.
func (mvn *Mvnormal) Density(
    x *mat64.Dense,
        // Each row is a case.
    out []float64) []float64 {

    n, p := x.Dims()
    assert(len(mvn.mean) == p,
        "Dimensionalities of input 'x' and distribution mismatch")

    out = use_slice(out, n)

    // Mahalanobis distance
    Mahalanobis(x, mvn.mean, mvn.cov, out)

    coef := 1.0 / math.Sqrt(
        math.Pow(2*math.Pi, float64(p)) * mvn.cov_det)

    for i, v := range out {
        out[i] = coef * math.Exp(-0.5 * v)
    }

    return out
}




// Random generates random samples from a multivariate normal distribution.
// The result is a matrix with each row being a case.
func (mvn *Mvnormal) Random(
    n int,
    _ RandomNumberGenerator,
    out *mat64.Dense) *mat64.Dense {

    out = use_matrix(out, n, len(mvn.mean))

    U, ok := mat64.CholeskyR(mvn.cov)
    assert(ok, "Cholesky failed")

    p := len(mvn.mean)

    norm := NewNormal(0.0, 1.0)
    z := mat64.NewDense(n, p, norm.Random(p*n, &RNG{}, nil))

    out.Mul(z, U)

    for row := 0; row < n; row++ {
        for col := 0; col < p; col++ {
            out.Set(row, col, out.At(row, col) + mvn.mean[col])
        }
    }
        // TODO: vectorize using out.RowView when it's available,
        // or a 'sweep' function.

    // Alternative algorithm: see function 'mvrnom' in the R package
    // MASS.

    return out
}


