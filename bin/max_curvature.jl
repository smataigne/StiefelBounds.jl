using LinearAlgebra, SkewLinearAlgebra, Plots

n = 4
p = 2
N = 50
B₁ = [1 0; 0 0]
B₂ = [0 0; 1 0]

B = zeros(n- p, p)
B .= B₁ + 2 * B₂
S_ref = [zeros(p, p) -B₁'; B₁ zeros(n - p, n - p)]
S = [zeros(p, p) -B'; B zeros(n - p, n - p)]

t_ref = LinRange(0, π , N)
t = LinRange(0, π / norm(B), N)
Frobenius = zeros(N)
geodesic = zeros(N)
Frobenius_ref = zeros(N)
geodesic_ref = zeros(N)
for i ∈ 1:N
    Frobenius_ref[i] = norm(exp(t_ref[i] * S_ref)[:, 1:p] - Matrix(1.0I,n, p)) 
    geodesic_ref[i] = t_ref[i] * norm(B₁)
    Frobenius[i] =norm(exp(t[i] * S)[:, 1:p] - Matrix(1.0I,n, p)) 
    geodesic[i] = t[i] * norm(B)
end
P = plot()
plot!(Frobenius_ref, geodesic_ref, label = "ref")
scatter!(Frobenius, geodesic, label = "in plane")
plot!(Frobenius, 2 * asin.(Frobenius ./2))
display(P)
