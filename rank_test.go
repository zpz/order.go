package stats

import (
	"fmt"
	"math"
	"testing"
)

func TestIntsRank(t *testing.T) {
	var ints = []int{5, 3, 4, 6}
	var ints_stable = []int{5, 3, 3, 4, 6, 6, 2}

	fmt.Println("Ints: ", ints)
	fmt.Println("Rank: ", IntRank(ints, nil))
	fmt.Println("Correct rank: 2, 0, 1, 3")
	fmt.Println("Ints: ", ints_stable)
	fmt.Println("Rank: ", IntStableRank(ints_stable, nil))
}

func TestFloat64sRank(t *testing.T) {
	var float64s = []float64{74.3, 26.2, 58., math.NaN(), 40.}
	var float64s_stable = []float64{74.3, 26.2, 58., 58., math.NaN(), math.NaN(), 58., 40.}

	fmt.Println("Float64s: ", float64s)
	fmt.Println("Rank: ", FloatRank(float64s, nil))
	fmt.Println("Correct Rank: 4, 1, 3, 0, 2")
	fmt.Println("Ints: ", float64s_stable)
	fmt.Println("Rank: ", FloatStableRank(float64s_stable, nil))
}

func TestStringsRank(t *testing.T) {
	var strings = []string{"ad", "ed", "af", "89"}
	var strings_stable = []string{"ad", "ed", "af", "af", "89", "af"}

	fmt.Println("Strings: ", strings)
	fmt.Println("Rank ", StringRank(strings, nil))
	fmt.Println("Correct Rank:  1 3 2 0")
	fmt.Println("Strings: ", strings_stable)
	fmt.Println("Rank: ", StringStableRank(strings_stable, nil))
}
