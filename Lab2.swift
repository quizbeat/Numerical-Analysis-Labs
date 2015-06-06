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
    
    var k = 0
    do {
        k++
        if (k > 100) {
            println("oops...")
            break
        }
        xPrev = x
        x = f(xPrev)
    //} while (1.0/(1.0-q) * fabs(x - xPrev)) > eps // 35 iterations
    } while (q/(1.0-q)) * fabs(x - xPrev) > eps // 8 iterations
    
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
    
    var f = zeros(2)

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
        f[0] = -f1(x1, x2)
        f[1] = -f2(x1, x2)
        let x = solve(J, f)
        x1 = x1Prev + x[0]
        x2 = x2Prev + x[1]
    } while (max(fabs(x1 - x1Prev), fabs(x2 - x2Prev)) > eps)
    
    println("\(k) iterations")
    return (x1, x2)
}

func dphi1x1(x1: Double, x2: Double) -> Double {
    return 0
}

func dphi1x2(x1: Double, x2: Double) -> Double {
    return -sin(x2)
}

func dphi2x1(x1: Double, x2: Double) -> Double {
    return cos(x1)
}

func dphi2x2(x1: Double, x2: Double) -> Double {
    return 0
}


func iterativeMethodForSystem(f1: (Double, Double) -> Double, f2: (Double, Double) -> Double,
                              x01: Double, a1: Double, b1: Double,
                              x02: Double, a2: Double, b2: Double,
                              eps: Double) -> (Double, Double) {
    var phi = zeros((2, 2))
    func computePhi(x1: Double, x2: Double) {
        phi[0, 0] = dphi1x1(x1, x2)
        phi[0, 1] = dphi1x2(x1, x2)
        phi[1, 0] = dphi2x1(x1, x2)
        phi[1, 1] = dphi2x2(x1, x2)
    }
    
    let d1 = fabs(a1 - b1) / 100
    let d2 = fabs(a2 - b2) / 100
    var q = -inf
                                
    for (var x1 = a1; x1 <= b1; x1 += d1) {
        for (var x2 = a2; x2 <= b2; x2 += d2) {
            computePhi(x1, x2)
            let detPhi = det(phi)
            if (detPhi > q) {
                q = detPhi
            }
        }
    }
    assert(q < 1, "q must be less than 1")
                            
    var x1 = x01
    var x1Prev = x1
    var x2 = x02
    var x2Prev = x1
                                
    let qk = q / (1.0 - q)
    var k = 0
                                
    do {
        k++
        x1Prev = x1
        x2Prev = x2
        x1 = f1(x1Prev, x2Prev)
        x2 = f2(x1Prev, x2Prev)
    } while (qk * max(fabs(x1 - x1Prev), fabs(x2 - x2Prev)) > eps)
                               
    println("\(k) iterations")
    return (x1, x2)
}