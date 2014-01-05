package stats

import (
    "github.com/gonum/matrix/mat64"
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




func use_matrix(in *mat64.Dense, rows, cols int) *mat64.Dense {
    if in == nil {
        out, ok := mat64.NewDense(rows, cols, nil)
        assert(ok == nil, "mat64.NewDense failed")
        return out
    }
    r, c := in.Dims()
    assert(r == rows && c == cols, "Provided matrix has wrong shape")
    return in
}


