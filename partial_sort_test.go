package stats

import (
	"fmt"
	"testing"
)

func TestMedian(t *testing.T) {
	data := []float64{2, 3, 1.6, 7, 5.3, 4.8}
	var nilslice []float64

	fmt.Println("Median of ", data, " is ",
		Median(append(nilslice, data...)))

	fmt.Println("Median of ", data[1:], " is ",
		Median(append(nilslice, data[1:]...)))
}
