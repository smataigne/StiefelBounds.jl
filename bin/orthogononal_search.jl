using LinearAlgebra, SkewLinearAlgebra, Plots
import .Manifold as m

n = 3

N = 10000

Frobenius = zeros(N)
geodesic = zeros(N)

for i ∈ 1:N
    A = skewhermitian!((i / N) * 2 * randn(n, n))
    Q = exp(A)
    Frobenius[i] = norm(exp(A) - I)
    geodesic[i] = norm(m.skewlog(Q))
end
t = LinRange(0, 2*sqrt(n), 500)
t2 = LinRange(0, 2, 100)
P = plot()
scatter!(Frobenius ./ (2*sqrt(n)), geodesic./(π * sqrt(n)), markerstrokewidth= 0, color=:skyblue2)
plot!(t./ (2*sqrt(n)),  2*sqrt(n) * asin.(t ./ (2*sqrt(n)))./(π * sqrt(n)), label = "Lower bound", linewidth = 2, color=:black)
plot!(t2./ (2*sqrt(n)),  2 * asin.(t2 ./ 2)./(π * sqrt(n)), label = "Lower bound", linewidth = 2, color=:black)
display(P)