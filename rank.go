package order



// Rank returns the sample ranks of the elements in the collection data.
// For example, Rank([]int{3, 5, 2, 6}) returns []int{1, 2, 0, 3}.
func Rank(data OrderInterface) []int {
    idx := make([]int, data.Len())
    RankFill(data, idx)
    return idx
}




// RankFill is the same as Rank, except that the return fills in the
// provided int slice, hence avoiding some memory allocation.
func RankFill(data OrderInterface, idx []int) {
    z := Order(data)
    order_to_rank(z, idx)
}





func StableRank(data OrderInterface) []int {
    idx := make([]int, data.Len())
    StableRankFill(data, idx)
    return idx
}




func StableRankFill(data OrderInterface, idx []int) {
    z := StableOrder(data)
    order_to_rank(z, idx)
}




func Float64sRank(x []float64) []int {
    idx := make([]int, len(x))
    Float64sRankFill(x, idx)
    return idx
}


func Float64sRankFill(x []float64, idx []int) {
    z := Float64sOrder(x)
    order_to_rank(z, idx)
}




func Float64sStableRank(x []float64) []int {
    idx := make([]int, len(x))
    Float64sStableRankFill(x, idx)
    return idx
}


func Float64sStableRankFill(x []float64, idx []int) {
    z := Float64sStableOrder(x)
    order_to_rank(z, idx)
}



func IntsRank(x []int) []int {
    idx := make([]int, len(x))
    IntsRankFill(x, idx)
    return idx
}


func IntsRankFill(x []int, idx []int) {
    z := IntsOrder(x)
    order_to_rank(z, idx)
}



func IntsStableRank(x []int) []int {
    idx := make([]int, len(x))
    IntsStableRankFill(x, idx)
    return idx
}


func IntsStableRankFill(x []int, idx []int) {
    z := IntsStableOrder(x)
    order_to_rank(z, idx)
}




func StringsRank(x []string) []int {
    idx := make([]int, len(x))
    StringsRankFill(x, idx)
    return idx
}


func StringsRankFill(x []string, idx []int) {
    z := StringsOrder(x)
    order_to_rank(z, idx)
}




func StringsStableRank(x []string) []int {
    idx := make([]int, len(x))
    StringsStableRankFill(x, idx)
    return idx
}


func StringsStableRankFill(x []string, idx []int) {
    z := StringsStableOrder(x)
    order_to_rank(z, idx)
}



func order_to_rank(in, out []int) {
    for i := range in {
        out[in[i]] = i
    }
}
