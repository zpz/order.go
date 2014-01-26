package stats

import (
	"sort"
)

// PartialSort sorts data so that at the n-th position
// is the element that would be there had data be fully sorted,
// and all elements on or after the n-th position are no smaller
// than any element before location n.
// However, on either side of position n, the elements need not be
// sorted.
func PartialSort(data sort.Interface, n int) {
	nth_element(data, n)
}

// nth_element rearranges elements in data such that at index n
// is the element that would have been in that position
// if data were sorted. All elements at index >= n are no smaller than
// any element at index < n.
func nth_element(data sort.Interface, n int) {
	nn := data.Len()
	if n < 0 || n >= nn {
		return
	}
	max_depth := 0
	for i := nn; i > 0; i >>= 1 {
		max_depth++
	}
	max_depth *= 2
	intro_select(data, 0, n, nn, max_depth)
}

func intro_select(data sort.Interface, a, n, b, max_depth int) {
	for b-a > 7 {
		if max_depth == 0 {
			heap_select(data, a, n+1, b)

			// Place the nth largest element in its final position.
			data.Swap(a, n)
			return
		}
		max_depth--
		mlo, mhi := doPivot(data, a, b)

		// The following logic needs some thinking.
		// Not totally sure.
		if mhi <= n {
			a = mhi
		} else if n < mlo {
			b = mlo
		} else {
			return
		}
	}
	if b-a > 1 {
		insertionSort(data, a, b)
	}
}

func heap_select(data sort.Interface, a, n, b int) {
	first := a
	lo := 0
	hi := n - a

	// Build heap with greatest element
	// on [a, n) at position a.
	for i := (hi - 1) / 2; i >= lo; i-- {
		siftDown(data, i, hi, first)
	}

	for i := n; i < b; i++ {
		if data.Less(i, a) {
			data.Swap(i, a)
			siftDown(data, lo, hi, first)
		}
	}
}

// (Taken from Go 1.2 std package sort.)
// Insertion sort
func insertionSort(data sort.Interface, a, b int) {
	for i := a + 1; i < b; i++ {
		for j := i; j > a && data.Less(j, j-1); j-- {
			data.Swap(j, j-1)
		}
	}
}

// (Taken from Go 1.2 std package sort.)
// medianOfThree moves the median of the three values data[a], data[b], data[c] into data[a].
func medianOfThree(data sort.Interface, a, b, c int) {
	m0 := b
	m1 := a
	m2 := c
	// bubble sort on 3 elements
	if data.Less(m1, m0) {
		data.Swap(m1, m0)
	}
	if data.Less(m2, m1) {
		data.Swap(m2, m1)
	}
	if data.Less(m1, m0) {
		data.Swap(m1, m0)
	}
	// now data[m0] <= data[m1] <= data[m2]
}

// (Taken from Go 1.2 std package sort.)
func swapRange(data sort.Interface, a, b, n int) {
	for i := 0; i < n; i++ {
		data.Swap(a+i, b+i)
	}
}

// (Taken from Go 1.2 std package sort.)
func doPivot(data sort.Interface, lo, hi int) (midlo, midhi int) {
	m := lo + (hi-lo)/2 // Written like this to avoid integer overflow.
	if hi-lo > 40 {
		// Tukey's ``Ninther,'' median of three medians of three.
		s := (hi - lo) / 8
		medianOfThree(data, lo, lo+s, lo+2*s)
		medianOfThree(data, m, m-s, m+s)
		medianOfThree(data, hi-1, hi-1-s, hi-1-2*s)
	}
	medianOfThree(data, lo, m, hi-1)

	// Invariants are:
	//	data[lo] = pivot (set up by ChoosePivot)
	//	data[lo <= i < a] = pivot
	//	data[a <= i < b] < pivot
	//	data[b <= i < c] is unexamined
	//	data[c <= i < d] > pivot
	//	data[d <= i < hi] = pivot
	//
	// Once b meets c, can swap the "= pivot" sections
	// into the middle of the slice.
	pivot := lo
	a, b, c, d := lo+1, lo+1, hi, hi
	for {
		for b < c {
			if data.Less(b, pivot) { // data[b] < pivot
				b++
			} else if !data.Less(pivot, b) { // data[b] = pivot
				data.Swap(a, b)
				a++
				b++
			} else {
				break
			}
		}
		for b < c {
			if data.Less(pivot, c-1) { // data[c-1] > pivot
				c--
			} else if !data.Less(c-1, pivot) { // data[c-1] = pivot
				data.Swap(c-1, d-1)
				c--
				d--
			} else {
				break
			}
		}
		if b >= c {
			break
		}
		// data[b] > pivot; data[c-1] < pivot
		data.Swap(b, c-1)
		b++
		c--
	}

	n := smaller(b-a, a-lo)
	swapRange(data, lo, b-n, n)

	n = smaller(hi-d, d-c)
	swapRange(data, c, hi-n, n)

	return lo + b - a, hi - (d - c)
}

// (Taken from Go 1.2 std package sort.)
// siftDown implements the heap property on data[lo, hi).
// first is an offset into the array where the root of the heap lies.
func siftDown(data sort.Interface, lo, hi, first int) {
	root := lo
	for {
		child := 2*root + 1
		if child >= hi {
			break
		}
		if child+1 < hi && data.Less(first+child, first+child+1) {
			child++
		}
		if !data.Less(first+root, first+child) {
			return
		}
		data.Swap(first+root, first+child)
		root = child
	}
}
