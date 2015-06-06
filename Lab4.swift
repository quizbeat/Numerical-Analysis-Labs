//
//  Lab4.swift
//  NA-Labs
//
//  Created by Nikita Makarov on 03/06/15.
//  Copyright (c) 2015 com.scott. All rights reserved.
//

import Foundation

func finiteDifferenceMethod(p: (Double) -> Double, q: (Double) -> Double, t: (Double) -> Double, g: (Double) -> Double, var x0: Double, x1: Double, alpha1: Double, alpha2: Double, alpha: Double, beta1: Double, beta2: Double, beta: Double, h: Double) -> ([Double], [Double]) {
    
    let n = Int(fabs(x1 - x0) / h) + 1
    
    var X = [Double]()
    for i in 0..<n {
        X.append(x0)
        x0 += h
    }
    
    var A: matrix = zeros((n, n))
    var b: ndarray = zeros(n)
    
    A[0, 0] = alpha1 - alpha2/h
    A[0, 1] = alpha2/h
    b[0] = alpha
    
    for i in 1...(n-2) {
        let x = X[i]
        A[i, i-1] = p(x)/(h*h) - q(x)/(2*h)
        A[i, i] = -2.0*p(x)/(h*h) + t(x)
        A[i, i+1] = p(x)/(h*h) + q(x)/(2*h)
        b[i] = g(x)
    }
    
    A[n-1, n-2] = beta1 - beta2/h
    A[n-1, n-1] = beta2/h
    b[n-1] = beta
    
    var Y = solve(A, b)
    
    return (X, Y.grid)
}