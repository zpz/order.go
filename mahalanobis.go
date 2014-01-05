package stats 

import (
    "github.com/gonum/matrix/mat64"
)


// Mahalanobis computes squared Mahalanobis distance
// between each row of x and the vector y, with covariance matrix sigma.
func Mahalanobis(
    x *mat64.Dense,
    y []float64,
    sigma *mat64.Dense,
    out []float64) []float64 {

    n, p := x.Dims()
    assert(p == len(y),
        "Dimensionalities of 'x' and 'y' mismatch")
    r, c := sigma.Dims()
    assert(r == p && c == p,
        "Dimensionalities of input data and sigma mismatch")

    out = use_slice(out, n)

    // Diff matrix between x and y.
    xt, _ := mat64.NewDense(p, n, nil)
    for row := 0; row < p; row++ {
        for col := 0; col < n; col++ {
            xt.Set(row, col, x.At(col, row) - y[row])
        }
    }

    yy := mat64.Solve(sigma, xt)

    // TODO: an 'element-wise multiplication' function
    // would be better; or if mat64.Dense has a 'RowView',
    // we can use Multiply row by row.
    for row := 0; row < p; row++ {
        for col := 0; col < n; col++ {
            yy.Set(row, col, yy.At(row, col) * xt.At(row, col))
        }
    }


    // Col sums of yy.
    for col := 0; col < n; col++ {
        out[col] = 0
        for row := 0; row < p; row++ {
            out[col] += yy.At(row, col)
        }
    }

    return out
}


