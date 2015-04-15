//
//  main.swift
//  swix
//
//  Created by Nikita Makarov on 7/9/14.
//  Copyright (c) 2014 com.scott. All rights reserved.
//

import Foundation

let projectFolder = "/Users/air/Documents/Numerical-Analysis/NA-Labs/"

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
func JacobiRotation(eps: Double = 0.3) {
    let folder = projectFolder + "lab1/jacobi/"
    let A: matrix = read_csv("A.csv", prefix: folder)
    
    let (eigValues, eigVectors) = JacobiRotations(A, eps)
    
    write_csv(eigValues, filename: "eigValues.csv", prefix: folder)
    write_csv(eigVectors, filename: "eigVectors.csv", prefix: folder)
}

// QR decomposition algorithm, tasks to do:
// 1. find eigenvalues
func qr(eps: Double = 0.01) {
    let folder = projectFolder + "lab1/qr/"
    var A: matrix = read_csv("A.csv", prefix: folder)
    
    var (Q, R) = qr(A)
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
        A = R.dot(Q)
        var underDiag = underDiagElems(A)
        e = √(∑(underDiag^2).grid)
        (Q, R) = qr(A)
        println(k)
    } while e > eps
    
    
    
    write_csv(Q, filename: "Q.csv", prefix: folder)
    write_csv(R, filename: "R.csv", prefix: folder)
}

//LUPDecomposition()
//TDMA()
//iterativeMethod()
//SeidelMethod()
//JacobiRotation()
qr()

