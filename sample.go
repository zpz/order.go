package stats

import (
    "math"
    "math/rand"
    "code.google.com/p/probab/dst"
)



func Sample(n, k int, idx []int) []int {
    if idx == nil {
        idx = make([]int, k)
    }
    for i := 0; i < k; i++ {
        idx[i] = rand.Intn(n)
    }
    return idx[0:k]
}




// WeightedSample returns a random sample of n integers in the set
// {0, 1,..., len(w) - 1}, with replacement, according to weights
// w assigned to the members of the set. The weights in w must sum to 1;
// this is assumed and not checked.
func WeightedSample(w []float64, n int, idx []int) []int {
    if idx == nil {
        idx = make([]int, n)
    }
    for i := 0; i < n; i++ {
        idx[i] = int(dst.ChoiceNext(w))
    }
        // FIXME: there might be room for improvement,
        // in a way that generates all n numbers in a batch.
    return idx[0:n]
}




// LogweightedSample returns a random sample of n integers in the set
// {0, 1,..., len(w) - 1}, with replacement, according to log-weights
// logw assigned to the members of the set.
// The exponentials of the log-weights in logw need not sum to 1,
// as normalization is performed internally.
func LogweightedSample(logw []float64, n int, idx []int) []int {
    w := make([]float64, len(logw))
    max, _ := Max(logw)
    var sum float64
    for i, lw := range logw {
        w[i] = math.Exp(lw - max)
        sum += w[i]
    }
    Scale(w, 1/sum, w)

    return WeightedSample(w, n, idx)
}


