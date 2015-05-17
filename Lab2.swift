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
        println("error = \(((1 / (1 - q)) * fabs(x - xPrev)))")
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

func df1x1(x1: Double, x2: Double) -> Double {
    return 1
}

func df1x2(x1: Double, x2: Double) -> Double {
    return sin(x2)
}

func df2x1(x1: Double, x2: Double) -> Double {
    return -cos(x1)
}

func df2x2(x1: Double, x2: Double) -> Double {
    return 1
}

func NewtonMethodForSystem(f1: (Double, Double) -> Double, f2: (Double, Double) -> Double, x01: Double, x02: Double, eps: Double) -> (Double, Double) {
    var A1 = zeros((2, 2))
    var A2 = zeros((2, 2))
    var J = zeros((2, 2))
    
    var x1 = x01
    var x1Prev = x1
    var x2 = x02
    var x2Prev = x1
    
    func computeJ(x1: Double, x2: Double) {
        J[0, 0] = df1x1(x1, x2)
        J[0, 1] = df1x2(x1, x2)
        J[1, 0] = df2x1(x1, x2)
        J[1, 1] = df2x2(x1, x2)
    }
    func computeA(x1: Double, x2: Double) {
        A1[0, 0] = f1(x1, x2)
        A1[0, 1] = df1x2(x1, x2)
        A1[1, 0] = f2(x1, x2)
        A1[1, 1] = df2x2(x1, x2)
        
        A2[0, 0] = df1x1(x1, x2)
        A2[0, 1] = f1(x1, x2)
        A2[1, 0] = df2x1(x1, x2)
        A2[1, 1] = f2(x1, x2)
    }
    
    var k = 0
    do {
        k++
        if (k > 100) {
            println("infinity cycle")
            break
        }
        x1Prev = x1
        x2Prev = x2
        computeJ(x1Prev, x2Prev)
        computeA(x1Prev, x2Prev)
        x1 = x1Prev - det(A1) / det(J)
        x2 = x2Prev - det(A2) / det(J)
    } while (max(fabs(x1 - x1Prev), fabs(x2 - x2Prev)) > eps)
    
    println("\(k) iterations")
    return (x1, x2)
}