//
//  Lab2.swift
//  NA-Labs
//
//  Created by Nikita Makarov on 15/04/15.
//  Copyright (c) 2015 com.scott. All rights reserved.
//

import Foundation

func d(f: (Double) -> Double, x: Double) -> Double {
    let dx = 1e-3
    return (f(x + dx) - f(x)) / dx
}

func d2(f: (Double) -> Double, x: Double) -> Double {
    let dx = 1e-3
    return (f(x - dx) - 2 * f(x) + f(x + dx)) / (dx * dx)
}

// f(x) = e^(2x) + 3x - 4 = 0
// equal form: x = (4 - e^(2x)) / 3
//         or: x = log(4 - 3x) / 2

func iterativeMethod(f: (Double) -> Double, a: Double, b: Double, eps: Double) -> Double {
    var xPrev = 0.0
    var x = (a + b) / 2
    var q = d(f, x)
    q /= 1 - q
    println("q = \(q)")
    do {
        xPrev = x
        x = f(xPrev)
        println(x)
    } while q * abs(x - xPrev) > eps
    return x
}

func NewtonMethod(f: (Double) -> Double, a: Double, b: Double, eps: Double) -> Double {
    var x = 0.0
    if f(a) * d2(f, a) > 0 {
        x = a
    }
    else if f(b) * d2(f, b) > 0 {
        x = b
    }
    else {
        assert(true, "can't find first approximation")
    }
    
    var k = 0
    var delta = inf
    do {
        k++
        delta = f(x) / d(f, x)
        x -= delta
    } while abs(delta) >= eps
    
    println("\(k) iterations")
    println(x)
    
    return x
}