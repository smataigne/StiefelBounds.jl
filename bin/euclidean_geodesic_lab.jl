using LinearAlgebra,SkewLinearAlgebra

A = π * sqrt(3)/3 * [0 1 -1;-1 0 1;1 -1 0]
B = π* sqrt(3)/3 *[1 1 1]

S = [2*A -B';B 0]

t = 0.4
display(exp(t * S))
display(exp(-t * A))
display(norm(exp(t * S)*Matrix(1I,4,3)*exp(-t * A)-Matrix(1I,4,3))^2)
display(12*sin((t*π/2))^2)
display((Matrix(1I,3,4)*exp(t * S)*Matrix(1I,4,3)*exp(-t * A)))
display(3 * cos(t*π))
#display(schur(exp(S*t)))
#display(schur(exp(-A*t)))