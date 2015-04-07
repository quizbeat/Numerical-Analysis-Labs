//
//  main.swift
//  swix
//
//  Created by Nikita Makarov on 7/9/14.
//  Copyright (c) 2014 com.scott. All rights reserved.
//


import Foundation

let LAB1_MATRICES = "/Users/air/Documents/Numerical-Analysis/NA-Labs/lab1/"

let A: matrix = read_csv("A.csv", prefix: LAB1_MATRICES + "lup/")
let b: ndarray = read_csv("b.csv", prefix: LAB1_MATRICES + "lup/")

LUPDecomposition(A, b)