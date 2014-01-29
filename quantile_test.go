package stats

import (
	"fmt"
	"testing"
)

func TestMedian(t *testing.T) {
	data := []float64{2, 3, 1.6, 7, 5.3, 4.8}

	fmt.Println("Median of ", data, " is ",
		FloatMedian(float_clone(data)))

	fmt.Println("Median of ", data[1:], " is ",
		FloatMedian(float_clone(data[1:])))
}
