//
//  HelperFunctions.swift
//  NA-Lab3
//
//  Created by Nikita Makarov on 27/05/15.
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

func der(f: (Double) -> Double) -> ((Double) -> Double) {
    let dx = 1e-6
    return { (x) in
        return (f(x + dx) - f(x)) / dx
    }
}

func der2(f: (Double) -> Double) -> ((Double) -> Double) {
    let dx = 1e-6
    return { (x) in
        return (f(x - dx) - 2 * f(x) + f(x + dx)) / (dx * dx)
    }
}
