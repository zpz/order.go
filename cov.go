package stats

import (
    "fmt"
    "github.com/gonum/floats"
)


// p is dimensionality;
// n is number of points.
func mean(x []float64, p, n int) []float64 {
    if p * n != len(x) {
        panic(fmt.Errorf("length of input data is wrong"))
    }

    v := make([]float64, p)

    k := 0
    for i := 0; i < n; i++ {
        for j := 0; j < p; j++ {
            v[j] = v[j] * float64(i) / float64(i + 1.0) +
                    x[k] / float64(i + 1.0)
            k++
        }
    }

    return v
}



func centerize(x []float64, p, n int) {
    x_mean := mean(x, p, n)
    for i, j := 0, 0; i < n; i++ {
        floats.Sub(x[j : (j+p)], x_mean)
        j += p
    }
}



// cov_xy returns elements of the cov matrix in column major.
func cov_xy(x []float64, p_x int, y []float64, p_y int, n int) []float64 {
    x_copy := append([]float64{}, x...)
    centerize(x_copy, p_x, n)
    y_copy := append([]float64{}, y...)
    centerize(y_copy, p_y, n)

    v := make([]float64, p_x * p_y)

    i_v := 0
    for i_y := 0; i_y < p_y; i_y++ {
        for i_x := 0; i_x < p_x; i_x++ {
            z := 0.0
            for k := 0; k < n; k++ {
                z += x_copy[k * p_x + i_x] * y_copy[k * p_y + i_y]
            }
            z /= float64(n - 1.0)
            v[i_v] = z
            i_v++
        }
    }

    return v
}



// cov_xx returns lower-triangular elements of the cov matrix,
// including the diagonal elements, in column major.
func cov_xx(x []float64, p, n int) []float64 {
    x_copy := append([]float64{}, x...)
    centerize(x_copy, p, n)

    v := make([]float64, (p * p + p)/2)
    i_v := 0
    for j := 0; j < n; j++ {
        for i := j; i < n; i++ {
            // Covariance between dim i and j:
            z := 0.0
            for k := 0; k < n; k++ {
                kk := k * p
                z += x_copy[kk + i] * x_copy[kk + j]
            }
            z /= float64(n - 1.0)
            v[i_v] = z
            i_v++
        }
    }

    return v
}


