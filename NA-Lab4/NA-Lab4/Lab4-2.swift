//
//  Lab4-2.swift
//  NA-Lab4
//
//  Created by Nikita Makarov on 29/05/15.
//  Copyright (c) 2015 quiz.corp. All rights reserved.
//

import Foundation

func diffEquationShootingMethod(f: func3arg, g: func3arg,
x0: Double, x1: Double,
y0: Double, y1: Double,
eta1: Double, eta2: Double,
h: Double, eps: Double) -> ([Double], [Double]) {
    var (X, Y, Z) = diffEquationRungeKuttaMethod(f, g, y0, eta1, x0, x1, h)
    var yt_1 = Y.last!
    
    (X, Y, Z) = diffEquationRungeKuttaMethod(f, g, y0, eta2, x0, x1, h)
    var yt_2 = Y.last!
    
    var eta_n1 = eta1
    var eta_n2 = eta2
    
    var error = fabs(yt_2 - y1)
    
    while (error > eps) {
        var eta_n = eta_n2 - (eta_n2 - eta_n1) / (yt_2 - yt_1) * (yt_2 - y1)
        (X, Y, Z) = diffEquationRungeKuttaMethod(f, g, y0, eta_n, x0, x1, h)
        yt_1 = yt_2
        yt_2 = Y.last!
        eta_n1 = eta_n2
        eta_n2 = eta_n
        error = fabs(yt_2 - y1)
    }
    
    return (X, Y)
}