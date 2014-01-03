package stats


// out can be x, or a prepared output slice, or nil,
// in which case an output slice is created.
func Shift(x []float64, amnt float64, out []float64) []float64 {
    if out == nil {
        out = make([]float64, len(x))
    }
    for i, val := range x {
        out[i] = val + amnt
    }
    return out
}




// out can be x, or a prepared output slice, or nil,
// in which case an output slice is created.
func Scale(x []float64, factor float64, out []float64) []float64 {
    if out == nil {
        out = make([]float64, len(x))
    }
    for i, val := range x {
        out[i] = val * factor
    }
    return out
}




// out can be x, or y, or a prepared output slice, or nil,
// in which case an output slice is created.
func Add(x, y []float64, out []float64) []float64 {
    if out == nil {
        out = make([]float64, len(x))
    }
    for i, val := range x {
        out[i] = val + y[i]
    }
    return out
}




// out can be x, or y, or a prepared output slice, or nil,
// in which case an output slice is created.
func AddScaled(x, y []float64, yfactor float64, out []float64) []float64 {
    if out == nil {
        out = make([]float64, len(x))
    }
    for i, val := range x {
        out[i] = val + y[i] * yfactor
    }
    return out
}




// out can be x, or y, or a prepared output slice, or nil,
// in which case an output slice is created.
func Subtract(x, y []float64, out []float64) []float64 {
    if out == nil {
        out = make([]float64, len(x))
    }
    for i, val := range x {
        out[i] = val - y[i]
    }
    return out
}




// out can be x, or y, or a prepared output slice, or nil,
// in which case an output slice is created.
func Multiply(x, y []float64, out []float64) []float64 {
    if out == nil {
        out = make([]float64, len(x))
    }
    for i, val := range x {
        out[i] = val * y[i]
    }
    return out
}




// out can be x, or y, or a prepared output slice, or nil,
// in which case an output slice is created.
func Divide(x, y []float64, out []float64) []float64 {
    if out == nil {
        out = make([]float64, len(x))
    }
    for i, val := range x {
        out[i] = val / y[i]
    }
    return out
}




// out can be x, or a prepared output slice, or nil,
// in which case an output slice is created.
func Apply(x []float64, f func(float64) float64, out []float64) []float64 {
    if out == nil {
        out = make([]float64, len(x))
    }
    for i, val := range x {
        out[i] = f(val)
    }
    return out
}


