//
//  Lab4-1.swift
//  NA-Lab4
//
//  Created by Nikita Makarov on 29/05/15.
//  Copyright (c) 2015 quiz.corp. All rights reserved.
//

import Foundation

typealias func3arg = (Double, Double, Double) -> Double

func diffEquationEulerMethod(f: func3arg, g: func3arg, y0: Double, z0: Double, a: Double, b: Double, h: Double) -> ([Double], [Double]) {
    let n = Int(fabs(b - a) / h)
    
    var X: [Double] = []
    var Y: [Double] = []
    var Z: [Double] = []
    
    var x = a
    var y = y0
    var z = z0
    
    X.append(x)
    Y.append(y)
    Z.append(z)
    
    for i in 1...n  {
        X.append(x + h)
        Y.append(Y[i - 1] + h * f(x, y, z))
        Z.append(Z[i - 1] + h * g(x, y, z))
        x += h
        y = Y[i]
        z = Z[i]
    }
    return (X, Y)
}

// Runge-Kutta method

func K(f: func3arg, g: func3arg, k: Int, h: Double, x: Double, y: Double, z: Double) -> Double {
    if (k == 1) {
        return h * f(x, y, z)
    }
    else if (k == 2) {
        return h * f(x + 0.5 * h, y + 0.5 * K(f, g, 1, h, x, y, z), z + 0.5 * L(f, g, 1, h, x, y, z))
    }
    else if (k == 3) {
        return h * f(x + 0.5 * h, y + 0.5 * K(f, g, 2, h, x, y, z), z + 0.5 * L(f, g, 2, h, x, y, z))
    }
    else {
        return h * f(x + h, y + K(f, g, 3, h, x, y, z), z + L(f, g, 3, h, x, y, z))
    }
}

func L(f: func3arg, g: func3arg, k: Int, h: Double, x: Double, y: Double, z: Double) -> Double {
    if (k == 1) {
        return h * g(x, y, z)
    }
    else if (k == 2) {
        return h * g(x + 0.5 * h, y + 0.5 * K(f, g, 1, h, x, y, z), z + 0.5 * L(f, g, 1, h, x, y, z))
    }
    else if (k == 3) {
        return h * g(x + 0.5 * h, y + 0.5 * K(f, g, 2, h, x, y, z), z + 0.5 * L(f, g, 2, h, x, y, z))
    }
    else {
        return h * g(x + h, y + K(f, g, 3, h, x, y, z), z + L(f, g, 3, h, x, y, z))
    }
}

func dy(f: func3arg, g: func3arg, h: Double, x: Double, y: Double, z: Double) -> Double {
    return (1/6.0) * (K(f, g, 1, h, x, y, z) + 2 * K(f, g, 2, h, x, y, z) + 2 * K(f, g, 3, h, x, y, z) + K(f, g, 4, h, x, y, z))
}

func dz(f: func3arg, g: func3arg, h: Double, x: Double, y: Double, z: Double) -> Double {
    return (1/6.0) * (L(f, g, 1, h, x, y, z) + 2 * L(f, g, 2, h, x, y, z) + 2 * L(f, g, 3, h, x, y, z) + L(f, g, 4, h, x, y, z))
}

func diffEquationRungeKuttaMethod(f: func3arg, g: func3arg, y0: Double, z0: Double, a: Double, b: Double, h: Double) -> ([Double], [Double], [Double]) {
    
    let n = Int(fabs(b - a) / h)
    
    var X: [Double] = []
    var Y: [Double] = []
    var Z: [Double] = []
    
    var x = a
    var y = y0
    var z = z0
    
    X.append(x)
    Y.append(y)
    Z.append(z)
    
    for i in 1...n  {
        X.append(x + h)
        Y.append(Y[i - 1] + dy(f, g, h, x, y, z))
        Z.append(Z[i - 1] + dz(f, g, h, x, y, z))
        x += h
        y = Y[i]
        z = Z[i]
    }
    return (X, Y, Z)
}

func diffEquationAdamsMethod(f: func3arg, g: func3arg, y0: Double, z0: Double, a: Double, b: Double, h: Double) -> ([Double], [Double]) {
    let b0 = a + 3 * h
    var (X, Y, Z) = diffEquationRungeKuttaMethod(f, g, y0, z0, a, b0, h)
    let n = Int(fabs(b - a) / h)
    
    for i in 4...n {
        X.append(X.last! + h)
        
        var yk = 55*f(X[i-1], Y[i-1], Z[i-1]) - 59*f(X[i-2], Y[i-2], Z[i-2])
        yk += 37*f(X[i-3], Y[i-3], Z[i-3]) - 9*f(X[i-4], Y[i-4], Z[i-4])
        Y.append(Y[i-1] + (h/24.0) * yk)
        
        var zk = 55*g(X[i-1], Y[i-1], Z[i-1]) - 59*g(X[i-2], Y[i-2], Z[i-2])
        zk += 37*g(X[i-3], Y[i-3], Z[i-3]) - 9*g(X[i-4], Y[i-4], Z[i-4])
        Z.append(Z[i-1] + (h/24.0) * zk)
    }
    
    return (X, Y)
}

func RungeRombergError(res1: Double, res2:Double, h1: Double, h2: Double, p: Int) -> Double {
    let k = h2 / h1
    return (res1 - res2) / (pow(k, Double(p)) - 1)
}