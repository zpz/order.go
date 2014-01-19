package stats

import (
	"math"
	"sort"
)

// OrderInterface is a collection type that satisfies order calculations.
// Typically this is a slice of structs that contain
// an ordered field, such as int/float64/string, and an int field.
// The int field is an internal helper.
// Comparison of elements in the collection should rely solely on the
// ordered field and should not involve the int field.
// The int field is used by the order function as a scratch field.
type OrderInterface interface {
	sort.Interface

	// Write j to the int field of the ith element.
	SetIndex(i, j int)

	// Get the int field of the ith element.
	GetIndex(i int) int
}

func order(
	data OrderInterface,
	idx []int,
	stable bool) []int {

	n := data.Len()

	for i := 0; i < n; i++ {
		data.SetIndex(i, i)
	}

	if stable {
		sort.Stable(data)
	} else {
		sort.Sort(data)
	}

	idx = use_int_slice(idx, n)
	for i := 0; i < n; i++ {
		idx[i] = data.GetIndex(i)
	}
	return idx
}

// Order_ returns an index array which would place the elements
// of the input collection, i.e. data, into ascending order.
// For example, conceptually,
// Order([]int{3, 5, 2, 6}) returns []int{2, 0, 1, 3}.
//
// If idx upon entry is nil, a new index slice is allocated, filled, and
// returned.
// If idx upon entry is not nil, it will be filled and returned.
// If idx is not nil but does not have enough space, an error will result.
//
// Upon return, the elements in data have been re-ordered.
func Order_(data OrderInterface, idx []int) []int {
	return order(data, idx, false)
}

func StableOrder_(data OrderInterface, idx []int) []int {
	return order(data, idx, true)
}

type float64_index struct {
	x   float64
	idx int
}

// Type float64_index_slide implements OrderInterface.
type float64_index_slice []float64_index

func (x float64_index_slice) Len() int { return len(x) }
func (x float64_index_slice) Less(i, j int) bool {
	return x[i].x < x[j].x || math.IsNaN(x[j].x)
}
func (x float64_index_slice) Swap(i, j int)      { x[i], x[j] = x[j], x[i] }
func (x float64_index_slice) SetIndex(i, j int)  { x[i].idx = j }
func (x float64_index_slice) GetIndex(i int) int { return x[i].idx }

func make_float64_index_slice(x []float64) float64_index_slice {
	xx := make([]float64_index, len(x))
	for i := range x {
		xx[i].x = x[i]
	}
	return float64_index_slice(xx)
}

// The input data x is not changed in this operation.
func Order(x []float64, idx []int) []int {
	return Order_(make_float64_index_slice(x), idx)
}

// The input data x is not changed in this operation.
func StableOrder(x []float64, idx []int) []int {
	return StableOrder_(make_float64_index_slice(x), idx)
}

type int_index struct {
	x, idx int
}

// int_index_slice implements OrderInterface.
type int_index_slice []int_index

func (x int_index_slice) Len() int           { return len(x) }
func (x int_index_slice) Less(i, j int) bool { return x[i].x < x[j].x }
func (x int_index_slice) Swap(i, j int)      { x[i], x[j] = x[j], x[i] }
func (x int_index_slice) SetIndex(i, j int)  { x[i].idx = j }
func (x int_index_slice) GetIndex(i int) int { return x[i].idx }

func make_int_index_slice(x []int) int_index_slice {
	xx := make([]int_index, len(x))
	for i := range x {
		xx[i].x = x[i]
	}
	return int_index_slice(xx)
}

func IntOrder(x []int, idx []int) []int {
	return Order_(make_int_index_slice(x), idx)
}

func IntStableOrder(x []int, idx []int) []int {
	return StableOrder_(make_int_index_slice(x), idx)
}

type string_index struct {
	x   string
	idx int
}

// string_index_slice implements OrderInterface.
type string_index_slice []string_index

func (x string_index_slice) Len() int           { return len(x) }
func (x string_index_slice) Less(i, j int) bool { return x[i].x < x[j].x }
func (x string_index_slice) Swap(i, j int)      { x[i], x[j] = x[j], x[i] }
func (x string_index_slice) SetIndex(i, j int)  { x[i].idx = j }
func (x string_index_slice) GetIndex(i int) int { return x[i].idx }

func make_string_index_slice(x []string) string_index_slice {
	xx := make([]string_index, len(x))
	for i := range x {
		xx[i].x = x[i]
	}
	return string_index_slice(xx)
}

func StringOrder(x []string, idx []int) []int {
	return Order_(make_string_index_slice(x), idx)
}

func StringStableOrder(x []string, idx []int) []int {
	return StableOrder_(make_string_index_slice(x), idx)
}
