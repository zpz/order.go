package stats

import (
    "math"
)



func Sum_(x Numeric) float64 {
    sum := 0.0
    for i, n := 0, x.Len(); i < n; i++ {
        sum += x.Get(i)
    }
    return sum
}



func Sum(x []float64) float64 {
    return Sum_(Float64Slice(x))
}




func Prod_(x Numeric) float64 {
    prod := 1.0
    for i, n := 0, x.Len(); i < n; i++ {
        prod *= x.Get(i)
    }
    return prod
}



func Prod(x []float64) float64 {
    return Prod_(Float64Slice(x))
}




func Min_(x Numeric) (float64, int) {
    min := x.Get(0)
    min_idx := 0
    for i, n := 1, x.Len(); i < n; i++ {
        v := x.Get(i)
        if v < min {
            min, min_idx = v, i
        }
    }
    return min, min_idx
}



func Min(x []float64) (float64, int) {
    return Min_(Float64Slice(x))
}






func Max_(x Numeric) (float64, int) {
    max := x.Get(0)
    max_idx := 0
    for i, n := 1, x.Len(); i < n; i++ {
        v := x.Get(i)
        if v > max {
            max, max_idx = v, i
        }
    }
    return max, max_idx
}



func Max(x []float64) (float64, int) {
    return Max_(Float64Slice(x))
}




func MinMax_(x Numeric) (float64, float64, int, int) {
    min := x.Get(0)
    max := min
    min_idx := 0
    max_idx := 0
    for i, n := 1, x.Len(); i < n; i++ {
        v := x.Get(i)
        if v < min {
            min, min_idx = v, i
        }
        if v > max {
            max, max_idx = v, i
        }
    }
    return min, max, min_idx, max_idx
}




func MinMax(x []float64) (float64, float64, int, int) {
    return MinMax_(Float64Slice(x))
}




func Range_(x Numeric) (float64, float64) {
    min, max, _, _ := MinMax_(x)
    return min, max
}



func Range(x []float64) (float64, float64) {
    return Range_(Float64Slice(x))
}





func Mean_(x Numeric) float64 {
    mean := x.Get(0)
    for i, n := 1, x.Len(); i < n; i++ {
        mean += (x.Get(i) - mean) / float64(i + 1)
    }
    return mean
}



func Mean(x []float64) float64 {
    return Mean_(Float64Slice(x))
}




func Average_(x Numeric, f func(float64) float64) float64 {
    res := 0.0
    for i, n := 0, x.Len(); i < n; i++ {
        v := f(x.Get(i))
        res += (v - res) / float64(i + 1.0)
    }
    return res
}



func Average(x []float64, f func(float64) float64) float64 {
    return Average_(Float64Slice(x), f)
}




// MSD_ computes mean squared deviation about a specified value.
func MSD_(x Numeric, center float64) float64 {
    return Average_(x, func(x float64) float64 {
        return (x - center) * (x - center) })
}





func MSD(x []float64, center float64) float64 {
    return MSD_(Float64Slice(x), center)
}




func Var_(x Numeric) (float64, float64) {
    mean, mean_xx := 0.0, 0.0
    n := x.Len()
    for i := 0; i < n; i++ {
        v := x.Get(i)
        mean += (v - mean) / float64(i + 1.0)
        mean_xx += (v * v - mean_xx) / float64(i + 1.0)
    }
    variance := (mean_xx - mean * mean) * float64(n) / float64(n-1)
    return variance, mean
}




func Var(x []float64) (float64, float64) {
    return Var_(Float64Slice(x))
}




func Sd_(x Numeric) (float64, float64) {
    v, m := Var_(x)
    return math.Sqrt(v), m
}



func Sd(x []float64) (float64, float64) {
    return Sd_(Float64Slice(x))
}




// out can be x, a prepared output slice, or nil,
// in which case an output slice is created.
func Center(x []float64, out []float64) []float64 {
    return Shift(x, -Mean(x), out)
}




// out can be x, a prepared output slice, or nil,
// in which case an output slice is created.
func Standardize(x []float64, out []float64) []float64 {
    if out == nil {
        out = make([]float64, len(x))
    }
    sd, mean := Sd(x)
    for i, val := range x {
        out[i] = (val - mean) / sd
    }
    return out
}




// out can be a prepared output slice, or nil,
// in which case an output slice is created.
func Generate(n int, f func(int) float64, out []float64) []float64 {
    if out == nil {
        out = make([]float64, n)
    }
    for i := 0; i < n; i++ {
        out[i] = f(i)
    }
    return out
}




// out can be a prepared output slice, or nil,
// in which case an output slice is created.
func Seq(from, to, step float64, out []float64) []float64 {
    n := int((to - from) / step) + 1
    return Generate(
        n,
        func (i int) float64 { return from + float64(i) * step },
        out)
}




func All(x []float64, f func(float64) bool) bool {
    for _, val := range x {
        if !f(val) {
            return false
        }
    }
    return true
}




func Any(x []float64, f func(float64) bool) bool {
    for _, val := range x {
        if f(val) {
            return true
        }
    }
    return false
}




func None(x []float64, f func(float64) bool) bool {
    return !Any(x, f)
}




// Count returns the number of elements satisfying specific criterion.
func Count(x []float64, f func(float64) bool) int {
    n := 0
    for _, val := range x {
        if f(val) {
            n++
        }
    }
    return n
}




// Which returns the indices of elements satisfying specific criterion.
func Which(x []float64, f func(float64) bool, out []int) []int {
    if out == nil {
        out = make([]int, 0, len(x))
            // FIXME: this could be wasteful in memory
            // if only a small number of elements pass the filter.
    } else {
        out = out[0:0]
    }
    for i, val := range x {
        if f(val) {
            out = append(out, i)
        }
    }
    return out
}




// PickByIndex picks elements by indices.
// out may be a prepared output slice or nil,
// in which case an output slice is created.
func PickByIndex(x []float64, idx []int, out []float64) []float64 {
    if out == nil {
        out = make([]float64, len(idx))
    }
    for i, j := range idx {
        out[i] = x[j]
    }
    return out
}




func Filter(x []float64, f func(float64) bool, out []float64) []float64 {
    if out == nil {
        out = make([]float64, 0, len(x))
            // FIXME: this could be wasteful in memory
            // if only a small number of elements pass the filter.
    } else {
        out = out[0:0]
    }
    for _, val := range x {
        if f(val) {
            out = append(out, val)
        }
    }
    return out
}


