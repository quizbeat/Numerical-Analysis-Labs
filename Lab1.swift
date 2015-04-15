//
//  Lab1.swift
//  NM-Labs
//
//  Created by Nikita Makarov on 06/03/15.
//  Copyright (c) 2015 {company name} All rights reserved.
//

import Foundation

let prefix = "/Users/air/Documents/Numerical-Analysis/NA-Labs/lab1/"

func TDMA(A: matrix, b: ndarray) -> ndarray {
    var P = zeros(A.rows)
    var Q = zeros(A.rows)
    let n = A.rows - 1
    
    P[0] = (-1) * A[0, 1] / A[0, 0]
    Q[0] = b[0] / A[0, 0]
    
    for i in 1...(n - 1) {
        P[i] = (-1) * (A[i, i + 1]) / (A[i, i] + A[i, i - 1] * P[i - 1])
        Q[i] = (b[i] - A[i, i - 1] * Q[i - 1])
        Q[i] /= (A[i, i] + A[i, i - 1] * P[i - 1])
    }
    
    P[n] = 0
    Q[n] = b[n] - A[n, n - 1] * Q[n - 1]
    Q[n] /= A[n, n] + A[n, n - 1] * P[n - 1]
    
    var x = zeros(A.rows)
    x[n] = Q[n]
    for (var i = n; i > 0; i--) {
        x[i - 1] = P[i - 1] * x[i] + Q[i - 1]
    }
    
    return x
}

func LUDecomposition(A: matrix) -> (L: matrix, U: matrix) {
    let n = A.rows - 1
    
    var L = ones((A.rows, A.rows))
    var U = zeros((A.rows, A.rows))
    
    for i in 0...n {
        for j in 0...n {
            var sum = 0.0
            for k in 0..<i {
                sum += L[i, k] * U[k, j]
            }
            U[i, j] = A[i, j] - sum
            sum = 0
            for k in 0..<i {
                sum += L[j, k] * U[k, i]
            }
            L[j, i] = (A[j, i] - sum) / U[i, i]
        }
    }
    return (L, U)
}

func solveWithLUP(L: matrix, U: matrix, P: matrix, b: ndarray) -> ndarray {
    // solving system Ax = b -> LUx = Pb
    // first step: Ly = Pb
    // second step: Ux = y, where x - solution
    let n = L.rows - 1
    var x = zeros(L.rows)
    var y = zeros(L.rows)
    let Pb = P.dot(b)
    
    // 1st step
    y[0] = Pb[0]
    for i in 1...n {
        var sum = 0.0
        for j in 0...(i - 1) {
            sum += L[i, j] * y[j]
        }
        y[i] = Pb[i] - sum
    }
    // 2nd step
    x[n] = y[n] / U[n, n]
    for (var i = n - 1; i >= 0; i--) {
        var sum = 0.0
        for j in (i + 1)...n {
            sum += U[i, j] * x[j]
        }
        x[i] = 1 / U[i, i] * (y[i] - sum)
    }
    return x
}

func inverseMatrixLUP(L: matrix, U: matrix, P: matrix) -> matrix {
    let n = L.rows - 1
    var invA = zeros_like(L)
    for j in 0...n {
        let ei = ind(L.rows, j)
        invA[0...n, j] = solveWithLUP(L, U, P, ei)
    }
    return invA
}

func LUPDecomposition(A: matrix, b: ndarray) -> (L: matrix, U: matrix, P: matrix) {
    let n = A.rows - 1
    
    var C = A.copy()
    var P = eye(A.rows)
    
    // calculating matrix C = L + U - E
    for i in 0...n {
        // find pivot
        var pivot = 0.0
        var pivotIndex = -1
        for row in i...n {
            let c = fabs(C[row, i])
            if (c > pivot) {
                pivot = c
                pivotIndex = row
            }
        }
        if pivot == 0 {
            println("zero matrix")
            println(C)
            return (zeros_like(A), zeros_like(A), zeros_like(A))
        }
        // swap
        swapRows(&C, pivotIndex, i)
        swapRows(&P, pivotIndex, i)
        // calculs
        for (var j = i + 1; j <= n; j++) {
            C[j, i] /= C[i, i]
            for (var k = i + 1; k <= n; k++) {
                C[j, k] -= C[j, i] * C[i, k]
            }
        }
    }

    // make L
    var L = eye(A.rows)
    for i in 1...n {
        for j in 0..<i {
            L[i, j] = C[i, j]
        }
    }
    
    // make U
    var U = zeros((A.rows, A.rows))
    for i in 0...n {
        for j in i...n {
            U[i, j] = C[i, j]
        }
    }
    
    return (L, U, P)
}

func iterativeMethod(var A: matrix, b: ndarray, eps: Double) -> ndarray {
    // also known as Jacobi method
    var alpha = zeros_like(A)
    var beta = zeros_like(b)
    
    var C = A.copy()
    let n = A.rows - 1
    
    // init alpha and beta
    for i in 0...n {
        for j in 0...n {
            // fix zero diag element
            if (i == j && A[i, j] == 0) {
                for k in (i + 1)...n {
                    if A[k, j] != 0 {
                        swapRows(&A, i, k)
                    }
                }
            }
            beta[i] = b[i] / A[i, i]
            if (i != j) {
                alpha[i, j] = (-1) * A[i, j] / A[i, i]
            }
            else {
                alpha[i, j] = 0
            }
        }
    }
    
    func norm_c(A: matrix) -> Double {
        var A_norm = 0.0
        for i in 0..<A.rows {
            let sum = norm(A[i, 0..<A.columns], ord: 1)
            if sum > A_norm {
                A_norm = sum
            }
        }
        return A_norm
    }
    
    let alphaNorm = norm(alpha, ord: inf)
    let k = alphaNorm / (1 - alphaNorm)
    
    var z = 0
    
    func sufficiencyConditionNotHold(x: ndarray, xPrev: ndarray) -> Double {
        return norm(x - xPrev, ord: inf)
    }
    
    func sufficiencyConditionHold(x: ndarray, xPrev: ndarray) -> Double {
        return k * sufficiencyConditionNotHold(x, xPrev)
    }
    
    var error: (ndarray, ndarray) -> Double
    
    if diagonallyDominant(A) {
        error = sufficiencyConditionHold
    }
    else {
        error = sufficiencyConditionNotHold
    }
    
    var x = beta.copy() // first approximation
    var xPrev = zeros_like(beta)
    
    var iteration = 0
    var epsK = 0.0
    
    do {
        iteration++
        xPrev = x.copy()
        x = beta + alpha.dot(xPrev)
        epsK = error(x, xPrev)
    } while (epsK > eps)
    
    println("\(iteration) iterations")
    write_csv(x, filename: "solution.csv", prefix: prefix)
    
    return x
}

func SeidelMethod(var A: matrix, b: ndarray, eps: Double) -> ndarray {
    // alpha = B + C
    // B is lower triangular matrix, B[diag] = 0
    // C is upper triangular matrix, C[diag] ≢ 0
    var alpha = zeros_like(A)
    var beta = zeros_like(b)
    
    var A_original = A.copy()
    let n = A.rows - 1
    
    // init alpha and beta
    for i in 0...n {
        for j in 0...n {
            // fix zero diag element
            if (i == j && A[i, j] == 0) {
                for k in (i + 1)...n {
                    if A[k, j] != 0 {
                        swapRows(&A, i, k)
                    }
                }
            }
            beta[i] = b[i] / A[i, i]
            if (i != j) {
                alpha[i, j] = (-1) * A[i, j] / A[i, i]
            }
            else {
                alpha[i, j] = 0
            }
        }
    }
    
    var C = zeros_like(alpha)
    for i in 0...n {
        for j in i...n {
            C[i, j] = alpha[i, j]
        }
    }
    
    let cNorm = norm(C, ord: inf)
    let alphaNorm = norm(alpha, ord: inf)
    let k = cNorm / (1 - alphaNorm)
    
    func sufficiencyConditionNotHold(x: ndarray, xPrev: ndarray) -> Double {
        return norm(x - xPrev, ord: inf)
    }
    
    func sufficiencyConditionHold(x: ndarray, xPrev: ndarray) -> Double {
        return k * sufficiencyConditionNotHold(x, xPrev)
    }
    
    var error: (ndarray, ndarray) -> Double
    
    if diagonallyDominant(A) {
        error = sufficiencyConditionHold
    }
    else {
        error = sufficiencyConditionNotHold
    }
    
    var x = zeros_like(b)
    var xPrev = beta.copy() // first approximation
    
    var iteration = 0
    var epsK = 0.0
    
    do {
        iteration++
        for i in 0...n {
            x[i] += beta[i]
            for j in 0..<i {
                x[i] += alpha[i, j] * x[j]
            }
            for j in i...n {
                x[i] += alpha[i, j] * xPrev[j]
            }
        }
        epsK = error(x, xPrev)
        xPrev = x.copy()
        x = zeros_like(beta)
    } while (epsK > eps)
    
    x = xPrev.copy()
    println("\(iteration) iterations")
    
    return x
}

func JacobiRotations(var A: matrix, eps: Double) -> (eigValues: ndarray, eigVectors: matrix) {
    let n = A.rows - 1
    let originalA = A.copy()
    var APrev = A.copy()
    var eigVectors = eye(A.rows)
    
    // approximation error
    func t(B: matrix) -> Double {
        let n = B.rows - 1
        var sum = 0.0
        for i in 1...n {
            for j in 0..<i {
                sum += B[i, j] * B[i, j]
            }
        }
        return sqrt(sum)
    }
    
    var iteration = -1
    do {
        iteration++
        // finding max not diagonal element
        var (mi, mj) = (0, n)
        for i in 0...n {
            for j in 0...n {
                if i == j {
                    continue
                }
                if (fabs(A[i, j]) > fabs(A[mi, mj])) {
                    (mi, mj) = (i, j)
                }
            }
        }
        
        let phiK = 0.5 * atan((2 * A[mi, mj]) / (A[mi, mi] - A[mj, mj]))
        var U = eye(A.rows)
        U[mi, mj] = -sin(phiK)
        U[mj, mi] = sin(phiK)
        U[mi, mi] = cos(phiK)
        U[mj, mj] = cos(phiK)
        
        eigVectors = eigVectors.dot(U)
        
        APrev = A.copy()
        A = (U.T.dot(APrev)).dot(U)
    } while (t(A) > eps)
    
    println("\(iteration) iterations")
    
    return (A["diag"], eigVectors)
}

func qr(var A: matrix) -> (Q: matrix, R: matrix) {
    let n = A.rows - 1
    var APrev = zeros_like(A)
    var Q = eye(A.rows)
    
    for j in 0...(n - 1) {
        APrev = A.copy()
        var v = APrev[0...n, j]
        for i in 0..<j {
            v[i] = 0
        }
        v[j] += sign(v[j]) * norm(v, ord: 2)
        let mv = reshape(v, (v.n, 1)) // vector v in matrix view
        let mvT = reshape(v, (1, v.n)) // vector vT in matrix view
        let num = mv.dot(mvT) // numerator matrix, v * vT
        let den = (mvT.dot(mv)).flat.grid[0] // denominator number, vT * v
        var H = eye(A.rows) - 2 * (num / den) // Householder matrix
        A = H.dot(APrev)
        Q = Q.dot(H)
    }
    return (Q, A)
}
