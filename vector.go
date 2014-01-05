package stats


// out can be x, or a prepared output slice, or nil,
// in which case an output slice is created.
func Shift(x []float64, amnt float64, out []float64) []float64 {
    out = use_slice(out, len(x))
    for i, val := range x {
        out[i] = val + amnt
    }
    return out
}




// out can be x, or a prepared output slice, or nil,
// in which case an output slice is created.
func Scale(x []float64, factor float64, out []float64) []float64 {
    out = use_slice(out, len(x))
    for i, val := range x {
        out[i] = val * factor
    }
    return out
}




// out can be x, or y, or a prepared output slice, or nil,
// in which case an output slice is created.
func Add(x, y []float64, out []float64) []float64 {
    out = use_slice(out, len(x))
    for i, val := range x {
        out[i] = val + y[i]
    }
    return out
}




// out can be x, or y, or a prepared output slice, or nil,
// in which case an output slice is created.
func AddScaled(x, y []float64, yfactor float64, out []float64) []float64 {
    out = use_slice(out, len(x))
    for i, val := range x {
        out[i] = val + y[i] * yfactor
    }
    return out
}




// out can be x, or y, or a prepared output slice, or nil,
// in which case an output slice is created.
func Subtract(x, y []float64, out []float64) []float64 {
    out = use_slice(out, len(x))
    for i, val := range x {
        out[i] = val - y[i]
    }
    return out
}




// out can be x, or y, or a prepared output slice, or nil,
// in which case an output slice is created.
func Multiply(x, y []float64, out []float64) []float64 {
    out = use_slice(out, len(x))
    for i, val := range x {
        out[i] = val * y[i]
    }
    return out
}




// out can be x, or y, or a prepared output slice, or nil,
// in which case an output slice is created.
func Divide(x, y []float64, out []float64) []float64 {
    out = use_slice(out, len(x))
    for i, val := range x {
        out[i] = val / y[i]
    }
    return out
}




// out can be x, or a prepared output slice, or nil,
// in which case an output slice is created.
func Transform(x []float64, f func(float64) float64, out []float64) []float64 {
    out = use_slice(out, len(x))
    for i, val := range x {
        out[i] = f(val)
    }
    return out
}




func Dot(x, y []float64) float64 {
    res := 0.0
    for i, v := range x {
        res += v * y[i]
    }
    return res
}







// out can be x, a prepared output slice, or nil,
// in which case an output slice is created.
func Center(x []float64, out []float64) []float64 {
    return Shift(x, -Mean(x), out)
}




// out can be x, a prepared output slice, or nil,
// in which case an output slice is created.
func Standardize(x []float64, out []float64) []float64 {
    out = use_slice(out, len(x))
    sd, mean := Sd(x)
    for i, val := range x {
        out[i] = (val - mean) / sd
    }
    return out
}




// out can be a prepared output slice, or nil,
// in which case an output slice is created.
func Generate(n int, f func(int) float64, out []float64) []float64 {
    out = use_slice(out, n)
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
    out = use_slice(out, len(idx))
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


