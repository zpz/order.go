package stats

import "sort"


type Numeric interface {
    Get(int) float64
    Set(int, float64)
    sort.Interface
}





// Float64Slice implements Numeric.
type Float64Slice []float64

func (f Float64Slice) Get(i int) float64 { return f[i] }
func (f Float64Slice) Set(i int, v float64) { f[i] = v }
func (f Float64Slice) Len() int { return len(f) }
func (f Float64Slice) Less(i, j int) bool { return f[i] < f[j] }
func (f Float64Slice) Swap(i, j int) { f[i], f[j] = f[j], f[i] }



// IntSlice implements Numeric.
type IntSlice []int

func (f IntSlice) Get(i int) float64 { return float64(f[i]) }
func (f IntSlice) Set(i int, v float64) { f[i] = int(v) }
func (f IntSlice) Len() int { return len(f) }
func (f IntSlice) Less(i, j int) bool { return f[i] < f[j] }
func (f IntSlice) Swap(i, j int) { f[i], f[j] = f[j], f[i] }

