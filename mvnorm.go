package stats

import (
    "math"
    "github.com/skelterjohn/go.matrix"
    "github.com/gonum/floats"
)



type Mvnorm struct {
    cov []float64
        // When a Mvnorm is used repeatedly
        // for distributions with the same dimensionality
        // (and varying mean and cov), the size of
        // 'cov' does not change, hence does not need re-allocation.
}



func MakeMvnorm() *Mvnorm {
    return new(Mvnorm)
}



// Density of a multivariate normal distribution.
// Input x contains one or more cases;
// output is densities for each case in x.
func (mvn *Mvnorm) Density(
    mean []float64,
        // Mean vector; its length reflects the dimensionality.
    cov []float64,
        // Lower triangular elements of the cov matrix,
        // in column major.
        // Length is (d*d + d)/2, where d is len(mean).
    x []float64,
        // One case after another,
        // each case is a vector of length equal to len(mean).
    return_log bool,
        // If true, return log density;
        // otherwise, return density.
    ) []float64 {

    z := make([]float64, len(x) / len(mean))
    mvn.DensityFill(mean, cov, x, return_log, z)
    return z
}




func (mvn *Mvnorm) DensityFill(
    mean []float64,
    cov []float64,
    x []float64,
    return_log bool,
    z []float64) {

    p := len(mean)
        // dimensionality

    Sigma := mvn.cov_matrix(cov)

    // Mahalanobis distance
    mahalanobis_fill(x, mean, Sigma, z)

    logdet := Sigma.Det()
    assert(logdet > 0, "Mvnorm.Density: cov matrix is not p.d.")
    logdet = math.Log(logdet)

    floats.AddConst(float64(p) * Ln2Pi + logdet, z)
    floats.Scale(-0.5, z)

    if !return_log {
        floats.Apply(math.Exp, z)
    }
}





// Random generates random samples from a multivariate normal distribution.
// The result is a slice of the required number of random cases
// stacked one after another.
func (mvn *Mvnorm) Random(
    mean []float64,
        // Mean vector; its length reflects the dimensionality.
    cov []float64,
        // Lower triangular elements of the cov matrix,
        // in column major.
        // Length is (d*d + d)/2, where d is len(mean).
    n int,
        // Number of random cases to be generated.
    ) []float64 {

    z := make([]float64, len(mean) * n)
    mvn.RandomFill(mean, cov, n, z)
    return z
}



func (mvn *Mvnorm) RandomFill(
    mean []float64,
    cov []float64,
    n int,
    z []float64) {

    if len(mean) == 1 {
        RNormFill(n, mean[0], cov[0], z)
        return
    }


    Sigma := mvn.cov_matrix(cov)

    C, err := Sigma.Cholesky()
    assert(err == nil, "Mvnorm.Random: Sigma.Cholesky() failed")

    p := len(mean)
    assert(Sigma.Rows() == p, "Mvnorm.Random: mu and Sigma sizes mismatch")

    zz := matrix.MakeDenseMatrix(RNorm(p * n, 0., 1.), p, n)
    zz, err = C.TimesDense(zz)
    assert(err == nil, "error in TimesDense")

    for col := 0; col < n; col++ {
        zz.BufferCol(col, z[col*p :])
        floats.Add(z[col*p : (col+1)*p], mean)
    }

    // Alternative algorithm: see function 'mvrnom' in the R package
    // MASS.
}




// cov_matrix takes a lower-tri cov slice and returns a cov matrix.
//
// Input is a slice containing lower triangular elements of the cov
// matrix, including elements on the main diagonal,
// stored in column major.
// Output is the full cov matrix.
func (mvn *Mvnorm) cov_matrix(
    cov []float64,
    ) *matrix.DenseMatrix {

    // Let d be the dimensionality,
    // then len(cov) is (d * d + d) / 2.

    p := covslice_ndim(len(cov))

    if cap(mvn.cov) < p * p {
        mvn.cov = make([]float64, p * p)
    }

    mvn.cov = mvn.cov[: p * p]

    covslice_expand_fill(cov, mvn.cov)
    Sigma := matrix.MakeDenseMatrix(mvn.cov[: p * p], p, p)
    return Sigma
}



