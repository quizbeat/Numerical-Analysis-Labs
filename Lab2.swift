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

func iterativeMethod(f: (Double) -> Double, a: Double, b: Double, eps: Double) -> Double {
    var xPrev = 0.0
    var x = (a + b) / 2
    var q = 0.0
    // find q
    for (var i = a; i <= b; i += eps) {
        let k = fabs(d(f, i))
        if (k > q) {
            q = k
        }
    }
    println("q = \(q)")
    
    var k = 0
    do {
        k++
        if (k > 100) {
            println("oops...")
            break
        }
        xPrev = x
        x = f(xPrev)
        println(x)
    } while (1 / (1 - q) * fabs(x - xPrev)) > eps // 35 iterations
    //} while (pow((1-q)/q,-1) * fabs(x - xPrev)) > eps // 8 iterations
    
    println("\(k) iterations")
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
        if (k > 100) {
            println("oops...")
            break
        }
        delta = f(x) / d(f, x)
        x -= delta
        println(x)
    } while fabs(delta) >= eps
    
    println("\(k) iterations")
    return x
}