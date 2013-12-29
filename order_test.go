package order

import (
	"fmt"
	"math"
	"testing"
)



func TestInts(t *testing.T) {
    var ints = []int{5, 3, 4, 6}
    var ints_stable = []int{5, 3, 3, 4, 6, 6, 2}

    fmt.Println("Ints: ", ints)
    fmt.Println("Order: ", IntsOrder(ints))
    fmt.Println("Correct order: 1, 2, 0, 3")
    fmt.Println("Ints: ", ints_stable)
    fmt.Println("Order: ", IntsStableOrder(ints_stable))
}


func TestFloat64s(t *testing.T) {
    var float64s = []float64{74.3, 26.2, 58., math.NaN(),  40.}
    var float64s_stable = []float64{74.3, 26.2, 58., 58., math.NaN(), math.NaN(), 58., 40.}

    fmt.Println("Float64s: ", float64s)
    fmt.Println("Order: ", Float64sOrder(float64s))
    fmt.Println("Correct order: 3, 1, 4, 2, 0")
    fmt.Println("Ints: ", float64s_stable)
    fmt.Println("Order: ", Float64sStableOrder(float64s_stable))
}


func TestStrings(t *testing.T) {
    var strings = []string{"ad", "ed", "af", "89"}
    var strings_stable = []string{"ad", "ed", "af", "af", "89", "af"}

    fmt.Println("Strings: ", strings)
    fmt.Println("Order: ", StringsOrder(strings))
    fmt.Println("Correct order:  3 0 2 1")
    fmt.Println("Strings: ", strings_stable)
    fmt.Println("Order: ", StringsStableOrder(strings_stable))
}
