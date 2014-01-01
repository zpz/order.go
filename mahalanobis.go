package stats 

import (
    "github.com/skelterjohn/go.matrix"
)


// Mahalanobis computes squared Mahalanobis distance
// between each row of X and the vector y, with covariance matrix Sigma.
func mahalanobis(
    x []float64,
        // One observation after another.
    y []float64,
    Sigma *matrix.DenseMatrix) []float64 {

    z := make([]float64, len(x) / len(y))
    mahalanobis_fill(x, y, Sigma, z)
    return z
}



func mahalanobis_fill(
    x []float64,
        // One observation after another.
    y []float64,
    Sigma *matrix.DenseMatrix,
    z []float64) {

    X := matrix.MakeDenseMatrix(x, len(x) / len(y), len(y))
        // This does not copy x.

    Xt := X.Transpose()
        // This does allocate new space.

    for row := 0; row < Xt.Rows(); row++ {
        for col := 0; col < Xt.Cols(); col++ {
            Xt.Set(row, col, Xt.Get(row, col) - y[row])
        }
    }

    yy, err := Sigma.Solve(Xt)
    assert(err == nil, "mahalanobis: solver failed")
    yy.ScaleMatrixDense(Xt)

    // Col sums of yy.
    for j := 0; j < yy.Cols(); j++ {
        z[j] = 0.
        for i := 0; i < yy.Rows(); i++ {
            z[j] += yy.Get(i, j)
        }
    }
}


