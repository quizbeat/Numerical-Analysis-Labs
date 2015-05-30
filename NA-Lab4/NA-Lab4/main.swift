//
//  main.swift
//  NA-Lab4
//
//  Created by Nikita Makarov on 29/05/15.
//  Copyright (c) 2015 quiz.corp. All rights reserved.
//

import Foundation

// MARK: Lab 4.1
//
// Second order equation:
// y'' + y'tg(x) - ycos^2(x) = 0
//
// Equal first order system:
// y' = z
// z' = ycos^2(x) - ztg(x)
//

func lab_4_1() {
    func f(x: Double, y: Double, z: Double) -> Double {
        return z
    }
    func g(x: Double, y: Double, z: Double) -> Double {
        return y * cos(x) * cos(x) - z * tan(x)
    }
    func correct(x: Double) -> Double {
        return exp(sin(x)) + exp(-sin(x))
    }
    
    let a = 0.0
    let b = 1.0
    let h = 0.1
    let y0 = 2.0
    let z0 = 0.0
    
    let (X, Y) = diffEquationRungeKuttaMethod(f, g, y0, z0, a, b, h)
    
    var error = 0.0;
    for j in 0...10 {
        error += fabs(Y[j] - correct(X[j]))
    }
    println(error)
}

// MARK: Lab 4.2

func lab_4_2() {

}


