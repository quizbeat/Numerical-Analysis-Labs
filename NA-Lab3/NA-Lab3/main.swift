//
//  main.swift
//  NA-Lab3
//
//  Created by Nikita Makarov on 19/05/15.
//  Copyright (c) 2015 quiz.corp. All rights reserved.
//

import Foundation

func apply(f: (Double) -> Double, x: [Double]) -> [Double] {
    var y: [Double] = []
    let n = x.count - 1
    for i in 0...n {
        y.append(f(x[i]))
    }
    return y
}

func y(x: Double) -> Double {
    return M_PI_2 - atan(x)
}

var Xa: [Double] = [-3, -1, 1, 3]
var Xb: [Double] = [-3, 0, 1, 3]

var Ya: [Double] = apply(y, Xa)
var Yb: [Double] = apply(y, Xb)

var x = -0.5

// Lagrange
func Lagrange() {
    println("Computing interpolation with Lagrange polynom.")
    println("f(x) = arcctg(x)")
    println("X = [-3; -1; 1; 3]")
    println("Result:")
    var L = interpolationLagrangePolynom(Xa, Ya)
    let yLa = y(x)
    let lLa = L(x)
    println("y(-0.5) = \(yLa)")
    println("L(-0.5) = \(lLa)")
    println("∆(L(-0.5)) = \(fabs(yLa - lLa))\n")
    
    println("X = [-3; 0; 1; 3]")
    println("Result:")
    L = interpolationLagrangePolynom(Xb, Yb)
    let yLb = y(x)
    let lLb = L(x)
    println("y(-0.5) = \(yLb)")
    println("L(-0.5) = \(lLb)")
    println("∆(L(-0.5)) = \(fabs(yLb - lLb))\n\n\n")
}

func Newton() {
    println("Computing interpolation with Newton polynom.")
    println("f(x) = arcctg(x)")
    println("X = [-3; -1; 1; 3]")
    println("Result:")
    var P = interpolationNewtonPolynom(Xa, Ya)
    let yPa = y(x)
    let lPa = P(x)
    println("y(-0.5) = \(yPa)")
    println("P(-0.5) = \(lPa)")
    println("∆(P(-0.5)) = \(fabs(yPa - lPa))\n")
    
    println("X = [-3; 0; 1; 3]")
    println("Result:")
    P = interpolationNewtonPolynom(Xb, Yb)
    let yPb = y(x)
    let lPb = P(x)
    println("y(-0.5) = \(yPb)")
    println("P(-0.5) = \(lPb)")
    println("∆(P(-0.5)) = \(fabs(yPb - lPb))\n")
}

Lagrange()
Newton()
