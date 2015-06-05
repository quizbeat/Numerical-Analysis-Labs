//
//  Lab4.swift
//  NA-Labs
//
//  Created by Nikita Makarov on 03/06/15.
//  Copyright (c) 2015 com.scott. All rights reserved.
//

import Foundation

func finiteDifferenceMethod(p: (Double) -> Double, q: (Double) -> Double, f: (Double) -> Double, var x0: Double, x1: Double, h: Double, y0: Double, y1: Double) -> ([Double], [Double]) {
    
    let n = Int(fabs(x1 - x0) / h) + 1
    
    var X = [Double]()
    for i in 0..<n {
        X.append(x0)
        x0 += h
    }
    
    var A: matrix = zeros((n, n))
    var b: ndarray = zeros(n)
    
    A[0, 0] = -1.0/h
    A[0, 1] = 1.0/h
    b[0] = y0
    
    for i in 1...(n-2) {
        A[i, i-1] = 1.0 - h*p(X[i])/2.0
        A[i, i] = -2.0 + q(X[i])*h*h
        A[i, i+1] = 1.0 + h*p(X[i])/2.0
        b[i] = f(X[i])*h*h
    }
    
    A[n-1, n-2] = -1.0/h
    A[n-1, n-1] = 1.0/h - 2
    b[n-1] = y1
    
    var Y = TDMA(A, b)
    
    return (X, Y.grid)
}