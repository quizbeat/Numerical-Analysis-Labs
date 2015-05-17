//
//  main.swift
//  swix
//
//  Created by Nikita Makarov on 7/9/14.
//  Copyright (c) 2014 com.scott. All rights reserved.
//

import Foundation

let projectFolder = "/Users/air/Documents/Numerical-Analysis/NA-Labs/"

// MARK: LAB 1

// LUP decomposition, tasks to do:
// 1. solve Ax = b
// 2. calculate inverse A
// 3. calculate det(A)
func LUPDecomposition() {
    let folder = projectFolder + "lab1/lup/"
    let A: matrix = read_csv("A.csv", prefix: folder)
    let b: ndarray = read_csv("b.csv", prefix: folder)
    
    let (L, U, P) = LUPDecomposition(A, b)
    let x = solveWithLUP(L, U, P, b)
    let invA = inverseMatrixLUP(L, U, P)
    let detA = ∏(U["diag"].grid)
    
    write_csv(L, filename: "L.csv", prefix: folder)
    write_csv(U, filename: "U.csv", prefix: folder)
    write_csv(P, filename: "P.csv", prefix: folder)
    write_csv(x, filename: "x.csv", prefix: folder)
    write_csv(invA, filename: "invA.csv", prefix: folder)
    println("det(A) = \(detA)")
}

// Tridiagonal matrix algorithm, tasks to do:
// 1. solve Ax = b
func TDMA() {
    let folder = projectFolder + "lab1/tdma/"
    let A: matrix = read_csv("A.csv", prefix: folder)
    let b: ndarray = read_csv("b.csv", prefix: folder)
    
    let x = TDMA(A, b)
    
    write_csv(x, filename: "x.csv", prefix: folder)
}

// Simple iterative method, tasks to do:
// 1. solve Ax = b
// 2. calculate amount of iterations
func iterativeMethod(eps: Double = 0.0001) {
    let folder = projectFolder + "lab1/iter_seidel/"
    let A: matrix = read_csv("A.csv", prefix: folder)
    let b: ndarray = read_csv("b.csv", prefix: folder)
    
    let x = iterativeMethod(A, b, eps)
    
    write_csv(x, filename: "x.csv", prefix: folder)
}

// Gauss-Seidel method, tasks to do:
// 1. solve Ax = b
// 2. calculate amount of iterations
func SeidelMethod(eps: Double = 0.0001) {
    let folder = projectFolder + "lab1/iter_seidel/"
    let A: matrix = read_csv("A.csv", prefix: folder)
    let b: ndarray = read_csv("b.csv", prefix: folder)
    
    let x = SeidelMethod(A, b, eps)
    
    write_csv(x, filename: "x.csv", prefix: folder)
}

// Jacobi rotations method, tasks to do:
// 1. find eigenvalues
// 2. find eigenvectors
func JacobiRotation(eps: Double = 0.0001) {
    let folder = projectFolder + "lab1/jacobi/"
    let A: matrix = read_csv("A.csv", prefix: folder)
    
    let (eigValues, eigVectors) = JacobiRotations(A, eps)
    
    //println(eigValues)
    //println(eigVectors)
    
    let n = A.rows - 1
    
    for i in 0...n {
        let l = eigValues[i]
        let v = eigVectors[0...n, i]
        println(A.dot(v)) // Ax
        println(l * v)   // lambda x
        println(" ")
    }
    
    write_csv(eigValues, filename: "eigValues.csv", prefix: folder)
    write_csv(eigVectors, filename: "eigVectors.csv", prefix: folder)
}

// QR decomposition algorithm, tasks to do:
// 1. find eigenvalues
func qr(eps: Double = 0.01) {
    let folder = projectFolder + "lab1/qr/"
    var A: matrix = read_csv("A.csv", prefix: folder)
    
    var (Q, R) = QRDecomposition(A)
    var APrev = A.copy()
    var e = inf
    let n = A.rows - 1
    
    func underDiagElems(A: matrix) -> ndarray {
        var underDiag = [Double]()
        for j in 0..<n {
            for i in (j + 1)...n {
                underDiag.append(A[i, j])
            }
        }
        return asarray(underDiag)
    }
    var k = 0
    
    do {
        k++
        if k == 100 {
            break
        }
        A = R.dot(Q)
        var underDiag = underDiagElems(A)
        e = √(∑(underDiag^2).grid)
        (Q, R) = QRDecomposition(A)
    } while e > eps
    
    if k == 100 {
        for j in 0...(n - 1) {
            var complexLambda = false
            for i in (j + 1)...n {
                if abs(A[i, j]) > eps {
                    complexLambda = true
                    let a = 1
                    let b = -(A[j, j] + A[j + 1, j + 1])
                    let c = (A[j, j] * A[j + 1, j + 1]) - (A[j + 1, j] * A[j, j + 1])
                    let D = b * b - 4 * a * c
                    //println(D)
                    if D > 0 {
                        var lambda = 0.0
                        lambda = (-b + sqrt(D)) / 2 * a
                        println("lambda = \(lambda)")
                        lambda = (-b - sqrt(D) / 2 * a)
                        println("lambda = \(lambda)")
                    }
                    else {
                        var lambdaRe = 0.0
                        var lambdaIm = 0.0
                        lambdaRe = -b / 2 * a
                        lambdaIm = sqrt(abs(D)) / 2 * a
                        println("lambda = \(lambdaRe) - \(lambdaIm) i")
                        println("lambda = \(lambdaRe) + \(lambdaIm) i")
                    }
                    break
                }
            }
            if (!complexLambda) {
                println("lambda = \(A[j, j])")
            }
        }
    }
    
    //write_csv(Q, filename: "Q.csv", prefix: folder)
    //write_csv(R, filename: "R.csv", prefix: folder)
}

//LUPDecomposition()
//TDMA()
//iterativeMethod()
//SeidelMethod()
//JacobiRotation()
//qr()



// MARK: LAB 2

func square(x: Double) -> Double {
    return x * x
}

func cube(x: Double) -> Double {
    return x * x * x
}

// test functions for 2.1
func f(x: Double) -> Double {
    return exp(x) - cube(x) + 3*square(x) - 2*x - 3
}

// equal form of f11
func phi(x: Double) -> Double {
    return log(cube(x) - 3*square(x) + 2*x + 3)
}

func SashaPhi(x: Double) -> Double {
    return sqrt((sin(x) + 0.5) / 2.0)
}

// Newton method, tasks to do:
// 1. find root of function
//let x = NewtonMethod(f, 0, 1.5, 1e-3)

// Simple iterative method
//let x = iterativeMethod(SashaPhi, 0, 2, 1e-3)

//println(x)

// test functions for part 2.2

func f1(x1: Double, x2: Double) -> Double {
    return x1 - cos(x2) - 2
}

func f2(x1: Double, x2: Double) -> Double {
    return x2 - sin(x1) - 2
}

let (x1, x2) = NewtonMethodForSystem(f1, f2, 1.25, 2.75, 1e-3)
println("x1 = \(x1)\nx2 = \(x2)")
