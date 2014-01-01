package stats

import (
    "math"
    "math/rand"
    "github.com/gonum/floats"
    "fmt"
)



// DNorm computes the probability density of given
// values in a normal distribution.
func DNorm(
    x []float64,
    mean, variance float64, return_log bool) []float64 {

    z := make([]float64, len(x))
    DNormFill(x, mean, variance, return_log, z)
    return z
}




// DNormFill is the same as DNorm except that result is
// returned in the provided slice.
func DNormFill(
    x []float64,
    mean, variance float64, return_log bool,
    z []float64) {

    a := -0.5 * (Ln2Pi + math.Log(variance))
    for i := range x {
        z[i] = a - 0.5 * (x[i] - mean) * (x[i] - mean) / variance
    }

    if return_log {
        return
    }

    floats.Apply(math.Exp, z[0 : len(x) - 1])
}





// RNorm generates n random numbers from the normal distribution.
func RNorm(n int, mean, variance float64) []float64 {
    z := make([]float64, n)
    RNormFill(n, mean, variance, z)
    return z
}




// RNormFill is the same as RNorm except that the result
// is returned in the provided slice.
func RNormFill(n int, mean, variance float64, z []float64) {
    assert(variance > 0., fmt.Sprintf("RNormFill: variance is %v", variance))
    sd := math.Sqrt(variance)
    for i := 0; i < n; i++ {
        z[i] = rand.NormFloat64() * sd + mean
    }
}


