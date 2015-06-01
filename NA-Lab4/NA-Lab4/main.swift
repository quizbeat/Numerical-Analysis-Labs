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

func computeError(X: [Double], Y: [Double], correct: (Double) -> Double) -> Double {
    var error = 0.0
    for i in 0...(X.count - 1) {
        error += fabs(Y[i] - correct(X[i]))
    }
    return error
}

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
    
    var X = [Double]()
    var Y = [Double]()
    var Z = [Double]()
    
    var X2 = [Double]()
    var Y2 = [Double]()
    var Z2 = [Double]()
    
    let h2 = h / 2.0

    println("Solving differential equation: y'' + y'tg(x) - ycos^2(x) = 0")
    println("Conditions: x∈[\(a); \(b)], y(0)=\(y0), y'(0)=\(z0)\n")
    
    // Euler
    println("==============================")
    println("Computing with Euler method...")
    (X, Y) = diffEquationEulerMethod(f, g, y0, z0, a, b, h)
    
    println("Result for Euler method:")
    println("Xi\tYi")
    for i in 0...(X.count-1) {
        println("\(X[i])\t\(Y[i])")
    }
    println(" ")
    
    let errorEuler = computeError(X, Y, correct)
    println("Difference with correct solution: \(errorEuler)\n")
    
    println("Error with Runge-Romberg method:")
    (X2, Y2) = diffEquationEulerMethod(f, g, y0, z0, a, b, h2)
    var j = 0
    for i in 0...(X.count-1) {
        let currentError = RungeRombergError(Y[i], Y2[j], h, h2, 1) // какое нужно P ??
        println("X_\(i): \(currentError)")
        j += 2
    }
    println("==============================\n")
    
    X.removeAll(keepCapacity: false)
    Y.removeAll(keepCapacity: false)
    
    X2.removeAll(keepCapacity: false)
    Y2.removeAll(keepCapacity: false)

    // Runge-Kutta
    println("==============================")
    println("Computing with Runge-Kutta method...")
    (X, Y, Z) = diffEquationRungeKuttaMethod(f, g, y0, z0, a, b, h)
    
    println("Result for Runge-Kutta method:")
    println("Xi\tYi")
    for i in 0...(X.count-1) {
        println("\(X[i])\t\(Y[i])")
    }
    println(" ")
    
    let errorRungeKutta = computeError(X, Y, correct)
    println("Difference with correct solution: \(errorRungeKutta)\n")
    
    println("Error with Runge-Romberg method:")
    (X2, Y2, Z2) = diffEquationRungeKuttaMethod(f, g, y0, z0, a, b, h2)
    j = 0
    for i in 0...(X.count-1) {
        let currentError = RungeRombergError(Y[i], Y2[j], h, h2, 1) // какое нужно P ??
        println("X_\(i): \(currentError)")
        j += 2
    }
    println("==============================\n")
    
    X.removeAll(keepCapacity: false)
    Y.removeAll(keepCapacity: false)
    
    
    // Adams
    println("==============================")
    println("Computing with Runge-Kutta method...")
    (X, Y) = diffEquationAdamsMethod(f, g, y0, z0, a, b, h)
    
    println("Result for Adams method:")
    println("Xi\tYi")
    for i in 0...(X.count-1) {
        println("\(X[i])\t\(Y[i])")
    }
    println(" ")
    
    let errorAdams = computeError(X, Y, correct)
    println("Difference with correct solution: \(errorAdams)\n")
    
    println("Error with Runge-Romberg method:")
    (X2, Y2) = diffEquationAdamsMethod(f, g, y0, z0, a, b, h2)
    j = 0
    for i in 0...(X.count-1) {
        let currentError = RungeRombergError(Y[i], Y2[j], h, h2, 1) // какое нужно P ??
        println("X_\(i): \(currentError)")
        j += 2
    }
    println("==============================")
    
    X.removeAll(keepCapacity: false)
    Y.removeAll(keepCapacity: false)
    Z.removeAll(keepCapacity: false)
    
    X2.removeAll(keepCapacity: false)
    Y2.removeAll(keepCapacity: false)
    Z2.removeAll(keepCapacity: false)
}

//lab_4_1()

// MARK: Lab 4.2

func lab_4_2() {

}

