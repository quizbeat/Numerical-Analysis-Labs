//
//  Lab3-5.swift
//  NA-Lab3
//
//  Created by Nikita Makarov on 23/05/15.
//  Copyright (c) 2015 quiz.corp. All rights reserved.
//

import Foundation

func integralRectangleMethod(f: (Double) -> Double, h: Double, x0: Double, x1: Double) -> Double {
    var N = (x1 - x0) / h
    var result = 0.0
    var x = x0 + h
    for i in 1...Int(N) {
        result += f(((x - h) + x) / 2)
        x += h
    }
    return h * result
}

func integralTrapezoidMethod(f: (Double) -> Double, h: Double, x0: Double, x1: Double) -> Double {
    var N = Int((x1 - x0) / h)
    var result = 0.0
    var x = x0 + h
    for i in 1...N {
        result += (f(x - h) + f(x))
        x += h
    }
    return 0.5 * h * result
}

func integralSimpsonMethod(f: (Double) -> Double, h: Double, x0: Double, x1: Double) -> Double {
    var result = 0.0
    var x = x0 + 2 * h
    while x <= x1 {
        result += f(x - 2 * h) + 4 * f(x - h) + f(x)
        x += 2 * h
    }
    return h / 3.0 * result
}

func integralRungeRombergMethod(res1: Double, res2:Double, h1: Double, h2: Double, p: Int) -> Double {
    let k = h2 / h1
    return res1 + (res1 - res2) / (pow(k, Double(p)) - 1)
}