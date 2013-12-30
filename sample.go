package stats

import (
    "math"
    "math/rand"
    "github.com/gonum/floats"
    "code.google.com/p/probab/dst"
)



func Sample(n, k int) []int {
    idx := make([]int, k)
    SampleFill(n, k, idx)
    return idx
}



func SampleFill(n, k int, idx []int) {
    for i := 0; i < k; i++ {
        idx[i] = rand.Intn(n)
    }
}




// WeightedSample returns a random sample of n integers in the set
// {0, 1,..., len(w) - 1}, with replacement, according to weights
// w assigned to the members of the set. The weights in w must sum to 1;
// this is assumed and not checked.
func WeightedSample(w []float64, n int) []int {
    z := make([]int, n)
    WeightedSampleFill(w, n, z)
    return z
}




// WeightedSampleFill returns a random sample in the provide slice idx;
// otherwise it's identical to WeightedSample.
func WeightedSampleFill(w []float64, n int, idx []int) {
    for i := 0; i < n; i++ {
        idx[i] = int(dst.ChoiceNext(w))
    }
        // FIXME: there might be room for improvement,
        // in a way that generates all n numbers in a batch.
}




// LogweightedSample returns a random sample of n integers in the set
// {0, 1,..., len(w) - 1}, with replacement, according to log-weights
// logw assigned to the members of the set.
// The exponentials of the log-weights in logw need not sum to 1,
// as normalization is performed internally.
func LogweightedSample(logw []float64, n int) []int {
    z := make([]int, n)
    LogweightedSampleFill(logw, n, z)
    return z
}




// LogweightedSampleFill returns a random sample in the provided slide idx;
// otherwise it's identical to LogSample.
func LogweightedSampleFill(logw []float64, n int, idx []int) {
    w := make([]float64, len(logw))

    max, _ := floats.Max(logw)
    var sum float64
    for i, lw := range logw {
        w[i] = math.Exp(lw - max)
        sum += w[i]
    }
    floats.Scale(1/sum, w)

    WeightedSampleFill(w, n, idx)
}


