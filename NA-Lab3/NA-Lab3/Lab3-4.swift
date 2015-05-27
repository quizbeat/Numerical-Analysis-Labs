//
//  Lab3-4.swift
//  NA-Lab3
//
//  Created by Nikita Makarov on 21/05/15.
//  Copyright (c) 2015 quiz.corp. All rights reserved.
//

import Foundation

func derivative(X: [Double], Y: [Double], rate: Int) -> ((Double) -> Double) {
    assert(rate == 1 || rate == 2, "rate must be equal 1 or 2")
    return { (x) in
        // find range
        var j = -1
        for i in 1...X.count {
            if (X[i - 1] <= x && x <= X[i]) {
                j = i - 1
                break
            }
        }
        if (rate == 1) {
            return (Y[j + 1] - Y[j]) / (X[j + 1] - X[j])
        }
        else {
            assert(j + 2 < X.count, "out of range")
            var res = 0.0
            res = (Y[j + 2] - Y[j + 1]) / (X[j + 2] - X[j + 1])
            res -= (Y[j + 1] - Y[j]) / (X[j + 1] - X[j])
            res /= X[j + 2] - X[j]
            return 2 * res
        }
    }
}