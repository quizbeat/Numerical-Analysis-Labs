//
//  File.swift
//  NA-Labs
//
//  Created by Nikita Makarov on 19/05/15.
//  Copyright (c) 2015 com.scott. All rights reserved.
//

import Foundation

func interpolationCubicSpline(X: [Double], Y: [Double]) -> ((Double) -> Double) {
    /*
    struct splineFactors {
        var a: [Double] = []
        var b: [Double] = []
        var c: [Double] = []
        var d: [Double] = []
    }
    */
    func h(i: Int) -> Double {
        return X[i] - X[i - 1]
    }
    
    let size = X.count
    let n = size - 1
    
    //var spline = splineFactors() // [__, S1, S2, S3], size = 4, S0 inactive
    var a = Y // a_i = f(x_i)

    // compute c_i
    var C: matrix = zeros((size, size))
    C[0, 0] = 0
    C[0, 1] = h(1)
    C[n, n - 1] = h(n)
    C[n, n] = 0
    for i in 1...(n - 1) {
        C[i, i - 1] = h(i)
        C[i, i] = 2 * (h(i) + h(i + 1))
        C[i, i + 1] = h(i + 1)
    }
    
    var r = zeros(size)
    for i in 1...(n - 1) {
        r[i] = (((Y[i+1] - Y[i]) / h(i+1)) - ((Y[i] - Y[i-1]) / h(i)))
        r[i] *= 6
    }
    //var c = (TDMA(C, r)).grid
    var c = solve(C, r).grid

    // compute d_i
    var d = [Double](count: size, repeatedValue: 0.0)
    for i in 1...n {
        d[i] = (c[i] - c[i - 1]) / h(i)
    }
    
    // compute b_i
    var b = [Double](count: size, repeatedValue: 0.0)
    for i in 1...n {
        b[i] = (Y[i] - Y[i - 1]) / h(i)
        b[i] += h(i) * (2.0 * c[i] + c[i - 1]) / 6.0
    }

    return { (x) in
        assert(X[0] <= x && x <= X[n], "x value out of range.")
        var res = 0.0
        // find range
        var j = -1
        for i in 0...n {
            if x >= X[i] {
                j = i
                break
            }
        }
        
        for i in 0...n {
            println("a{\(i)} = \(a[i])")
            println("b{\(i)} = \(b[i])")
            println("c{\(i)} = \(c[i])")
            println("d{\(i)} = \(d[i])\n")
        }
        
        j++
        res += a[j]
        res += b[j] * (x - X[j])
        res += c[j] / 2 * pow(x - X[j], 2.0)
        res += d[j] / 6 * pow(x - X[j], 3.0)
        return res
    }
}

func interpolationOLSPolynom(X: [Double], Y: [Double], rate: Int) -> ((Double) -> Double) {
    var x: matrix = zeros((rate + 1, rate + 1))
    var a: ndarray = zeros(rate + 1)
    var y: ndarray = zeros(rate + 1)
    
    var ndX = asarray(X)
    var ndY = asarray(Y)
    let n = rate
    
    for i in 0...n {
        for j in 0...n {
            x[i, j] = sum(ndX^Double(i + j))
        }
        y[i] = sum(ndY * (ndX^Double(i)))
    }
    a = solve(x, y)
    return { (x) in
        var res = 0.0
        for i in 0...n {
            res += a[i] * pow(x, Double(i))
        }
        return res
    }
}
