package order

import "sort"

// OrderInterface is a collection type that satisfies order.
// Typically this is a slice of structs that contain
// an orderable field, such as int/float64/string, and an int field.
// The int field is a helper that should not be used by the caller.
// Comparison of elements in the collection should rely solely on the
// orderable field and should not involve the int field.
// The int field is used by the order function as a scratch field.
type OrderInterface interface {
	sort.Interface
	SetIndex(i, j int)
	// Write j to the int field of the ith element.
	GetIndex(i int) int
	// Get the int field of the ith element.
}

// Order returns an index array which would place the elements
// of the input collection, i.e. data, into ascending order.
// For example, Order([]int{3, 5, 2, 6}) returns []int{2, 0, 1, 3}.
func Order(data OrderInterface) []int {
	idx := make([]int, data.Len())
	OrderFill(data, idx)
	return idx
}

// OrderFill is the same as Order, except that the return fills in the
// provided int slice, hence avoiding some memory allocation.
func OrderFill(data OrderInterface, idx []int) {
	for i := 0; i < data.Len(); i++ {
		data.SetIndex(i, i)
	}
	sort.Sort(data)
	for i := 0; i < data.Len(); i++ {
		idx[i] = data.GetIndex(i)
	}
}



func StableOrder(data OrderInterface) []int {
    idx := make([]int, data.Len())
    StableOrderFill(data, idx)
    return idx
}




func StableOrderFill(data OrderInterface, idx []int) {
    for i := 0; i < data.Len(); i++ {
        data.SetIndex(i, i)
    }
    sort.Stable(data)
    for i := 0; i < data.Len(); i++ {
        idx[i] = data.GetIndex(i)
    }
}



type float64Index struct {
	x   float64
	idx int
}


// isNaN is a copy of math.IsNaN to avoid dependence on the math
// package.
func isNaN(x float64) bool { return x != x }

type float64IndexSlice []float64Index

func (x float64IndexSlice) Len() int { return len(x) }

func (x float64IndexSlice) Less(i, j int) bool {
	return x[i].x < x[j].x || isNaN(x[i].x) && !isNaN(x[j].x)
}

func (x float64IndexSlice) Swap(i, j int) { x[i], x[j] = x[j], x[i] }

func (x float64IndexSlice) SetIndex(i, j int) { x[i].idx = j }

func (x float64IndexSlice) GetIndex(i int) int { return x[i].idx }

func Float64sOrder(x []float64) []int {
	idx := make([]int, len(x))
	Float64sOrderFill(x, idx)
	return idx
}

func Float64sOrderFill(x []float64, idx []int) {
	xx := make([]float64Index, len(x))
	for i := range x {
		xx[i].x = x[i]
	}
	OrderFill(float64IndexSlice(xx), idx)
}



func Float64sStableOrder(x []float64) []int {
    idx := make([]int, len(x))
    Float64sStableOrderFill(x, idx)
    return idx
}


func Float64sStableOrderFill(x []float64, idx []int) {
    xx := make([]float64Index, len(x))
    for i := range x {
        xx[i].x = x[i]
    }
    StableOrderFill(float64IndexSlice(xx), idx)
}


type intIndex struct {
	x, idx int
}

type intIndexSlice []intIndex

func (x intIndexSlice) Len() int { return len(x) }

func (x intIndexSlice) Less(i, j int) bool { return x[i].x < x[j].x }

func (x intIndexSlice) Swap(i, j int) { x[i], x[j] = x[j], x[i] }

func (x intIndexSlice) SetIndex(i, j int) { x[i].idx = j }

func (x intIndexSlice) GetIndex(i int) int { return x[i].idx }

func IntsOrder(x []int) []int {
	idx := make([]int, len(x))
	IntsOrderFill(x, idx)
	return idx
}

func IntsOrderFill(x []int, idx []int) {
	xx := make([]intIndex, len(x))
	for i := range x {
		xx[i].x = x[i]
	}
	OrderFill(intIndexSlice(xx), idx)
}



func IntsStableOrder(x []int) []int {
    idx := make([]int, len(x))
    IntsStableOrderFill(x, idx)
    return idx
}


func IntsStableOrderFill(x []int, idx []int) {
    xx := make([]intIndex, len(x))
    for i := range x {
        xx[i].x = x[i]
    }
    StableOrderFill(intIndexSlice(xx), idx)
}


type stringIndex struct {
	x   string
	idx int
}

type stringIndexSlice []stringIndex

func (x stringIndexSlice) Len() int { return len(x) }

func (x stringIndexSlice) Less(i, j int) bool { return x[i].x < x[j].x }

func (x stringIndexSlice) Swap(i, j int) { x[i], x[j] = x[j], x[i] }

func (x stringIndexSlice) SetIndex(i, j int) { x[i].idx = j }

func (x stringIndexSlice) GetIndex(i int) int { return x[i].idx }

func StringsOrder(x []string) []int {
	idx := make([]int, len(x))
	StringsOrderFill(x, idx)
	return idx
}

func StringsOrderFill(x []string, idx []int) {
	xx := make([]stringIndex, len(x))
	for i := range x {
		xx[i].x = x[i]
	}
	OrderFill(stringIndexSlice(xx), idx)
}



func StringsStableOrder(x []string) []int {
    idx := make([]int, len(x))
    StringsStableOrderFill(x, idx)
    return idx
}


func StringsStableOrderFill(x []string, idx []int) {
    xx := make([]stringIndex, len(x))
    for i := range x {
        xx[i].x = x[i]
    }
    StableOrderFill(stringIndexSlice(xx), idx)
}


