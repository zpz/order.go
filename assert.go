package stats

import "fmt"


func assert(test bool, msg string) {
    if !test {
        panic(fmt.Errorf(msg))
    }
}


