//
//  Lab3.swift
//  NA-Lab3
//
//  Created by Nikita Makarov on 19/05/15.
//  Copyright (c) 2015 quiz.corp. All rights reserved.
//

import Foundation

func interpolationLagrangePolynom(X: [Double], Y: [Double]) -> ((Double) -> Double) {
    let n = X.count - 1
    func w(x: Double) -> Double {
        var res = 1.0
        for i in 0...n {
            res *= (x - X[i])
        }
        return res
    }
    func dw(x: Double, k: Int) -> Double {
        var res = 1.0
        for i in 0...n {
            if (i == k) {
                continue
            }
            res *= (X[k] - X[i])
        }
        return res
    }
    func computeL(w: (Double) -> Double, dw: (Double, Int) -> Double) -> ((Double) -> Double) {
        return { (x) in
            var res = 0.0
            for i in 0...n {
                res += Y[i] * (w(x) / ((x - X[i]) * dw(X[i], i)))
            }
            return res
        }
    }
    return computeL(w, dw)
}

func dividedDifferences(X: [Double], Y: [Double], i: Int, j: Int) -> Double {
    if (i == j - 1) {
        return (Y[i] - Y[j]) / (X[i] - X[j])
    }
    else {
        return (dividedDifferences(X, Y, i, j - 1) - dividedDifferences(X, Y, i + 1, j)) / (X[i] - X[j])
    }
}

func interpolationNewtonPolynom(X: [Double], Y: [Double]) -> ((Double) -> Double) {
    let n = X.count - 1
    func computeP() -> ((Double) -> Double) {
        return { (x) in
            var res = Y[0]
            for i in 1...n {
                var k = 1.0
                for var j = 0; j < i; j++ {
                    k *= (x - X[j])
                }
                res += k * dividedDifferences(X, Y, 0, i)
            }
            return res
        }
    }
    return computeP()
}