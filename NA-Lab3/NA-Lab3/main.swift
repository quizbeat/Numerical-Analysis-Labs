//
//  main.swift
//  NA-Lab3
//
//  Created by Nikita Makarov on 19/05/15.
//  Copyright (c) 2015 quiz.corp. All rights reserved.
//

import Foundation

// MARK: Lab 3.1

func y(x: Double) -> Double {
    return M_PI_2 - atan(x)
}

var Xa: [Double] = [-3, -1, 1, 3]
var Xb: [Double] = [-3, 0, 1, 3]

var Ya: [Double] = apply(y, Xa)
var Yb: [Double] = apply(y, Xb)

var x = -0.5

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

func lab_3_1() {
    Lagrange()
    Newton()
}



// MARK: Lab 3.4

//let dx = [0.0, 0.5, 1.0, 1.5, 2.0]
//let dy = [1.0, 1.3776, 1.5403, 1.5707, 1.5839]

//let dx = [0.0, 1.0, 2.0, 3.0, 4.0]
//let dy = [0.0, 2.0, 3.4142, 4.7321, 6.0]

func lab_3_4() {
    let Xi = [0.0, 0.1, 0.2, 0.3, 0.4]
    let Yi = [1.0, 1.1052, 1.2214, 1.3499, 1.4918]
    let targetX = 0.2
    
    let firstDerivative = derivative(Xi, Yi, 1)(targetX)
    let secondDerivative = derivative(Xi, Yi, 2)(targetX)
    
    println("Computing 1st and 2nd derivatives:")
    println("X: \(Xi)")
    println("Y: \(Yi)")
    println("Target X = \(targetX)\n")
    println("1st derivative equal \(firstDerivative)")
    println("2nd derivative equal \(secondDerivative)\n\n")
}



// MARK: Lab 3.5

func lab_3_5() {
    func f(x: Double) -> Double {
        return x / pow(3 * x + 4, 2.0)
    }
    
    let x0 = -1.0
    let x1 = 1.0
    let h1 = 0.5
    let h2 = 0.25
    
    println("Computing definite integral of the function:")
    println("y = x / (3x + 4)^2   from \(x0) to \(x1)\n")
    
    println("=====[Rectangles method]=====")
    println("computing with step h1 = \(h1)...")
    let rectH1 = integralRectangleMethod(f, h1, x0, x1)
    println("∫ydx = \(rectH1)\n")
    println("computing with step h2 = \(h2)...")
    let rectH2 = integralRectangleMethod(f, h2, x0, x1)
    println("∫ydx = \(rectH2)")
    println("=============================\n")
    
    println("=====[Trapezoids method]=====")
    println("computing with step h1 = \(h1)...")
    let trapezoidH1 = integralTrapezoidMethod(f, h1, x0, x1)
    println("∫ydx = \(trapezoidH1)\n")
    println("computing with step h2 = \(h2)...")
    let trapezoidH2 = integralTrapezoidMethod(f, h2, x0, x1)
    println("∫ydx = \(trapezoidH2)")
    println("=============================\n")
    
    println("======[Simpson method]=======")
    println("computing with step h1 = \(h1)...")
    let simpsonH1 = integralSimpsonMethod(f, h1, x0, x1)
    println("∫ydx = \(simpsonH1)\n")
    println("computing with step h2 = \(h2)...")
    let simpsonH2 = integralSimpsonMethod(f, h2, x0, x1)
    println("∫ydx = \(simpsonH2)")
    println("=============================\n")
    
    println("Checking error value by Runge-Romberg method...\n")
    let rungeRect = RungeRombergError(rectH1, rectH2, h1, h2, 2)
    println("Result for rectangles method:")
    println("Error = \(rungeRect)\n")
    
    let rungeTrapezoid = RungeRombergError(trapezoidH1, trapezoidH2, h1, h2, 2)
    println("Result for trapezoids method:")
    println("Error = \(rungeTrapezoid)\n")
    
    let rungeSimpson = RungeRombergError(simpsonH1, simpsonH2, h1, h2, 2)
    println("Result for Simpson method:")
    println("Error = \(rungeSimpson)\n")
}

// ==================================================//
//                ENTER COMMANDS HERE:               //
// ==================================================//

//lab_3_1()
//lab_3_4()
//lab_3_5()
