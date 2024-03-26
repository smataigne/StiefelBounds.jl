using LinearAlgebra, SkewLinearAlgebra, Plots, LaTeXStrings
import .Manifold as m

#Choose N, the number of samples, and (n,p) the size of the Stiefel manifold
N = 10000
n = 4
p = 2

Frobenius = zeros(N)             #Memory for Frobenius distances
distancesEuclidean = zeros(N)    #Memory for Euclidean geodesic distances
distancesCanonical = zeros(N)    #Memory for canonical geodesic distances
S = zeros(2p,2p)
for i ∈ 1:N
    U = m.randorthonormal(Float64, n, n)
    U₁ = m.Orthonormal(U.Q[:,1:p], false)
    U₂ =  m.orthogonalize!((U.Q * exp(skewhermitian!(10*(i / N) *randn(n,n))))[:,1:p])
    
    Frobenius[i] = norm(Matrix(U₁ - U₂))
    #Compute Euclidean distance
    S₁= m.StiefelVector(U₁, 0)
    S₂= m.StiefelVector(U₂, 0)
    m.buildcomplement(S₁)
    m.buildcomplement(S₂)
    Δ = m.logshooting(S₁, S₂)
    A = U₁'*Δ
    B = S₁.Θ.Q'*Δ
    distancesEuclidean[i] = sqrt(norm(A)^2 + norm(B)^2)

    #Compute canonical distance
    S₁= m.StiefelVector(U₁, 1)
    S₂= m.StiefelVector(U₂, 1)
    m.buildcomplement(S₁)
    m.buildcomplement(S₂)
    Δ, Q = m.logβ(S₁, S₂, 0.5, 0)[1:2]
    A = Matrix(U₁)'*Δ
    B = Matrix(Q' * Δ) 
    distancesCanonical[i] =  sqrt(0.5*norm(A)^2 + norm(B)^2)
end

P = plot(framestyle=:box, legend=:topleft,font="Computer Modern", tickfontfamily="Computer Modern",legendfont="Computer Modern",
    legendfontsize=12,yguidefontsize=17,xguidefontsize=17, xtickfontsize = 13, ytickfontsize=13,margin = 0.3Plots.cm)
 
t = LinRange(0, π*sqrt(p), 500)

scatter!( distancesEuclidean, distancesCanonical, markersize=2, markerstrokewidth=0, label = "Samples")
plot!(t, t.*(sqrt(2)/2), linewidth = 3, label = L"\frac{\sqrt{2}}{2}d_\mathrm{E}(U, \widetilde{U})",color = :black)
plot!(t, t, linewidth = 3, label = L"d_\mathrm{E}(U, \widetilde{U})", color =:green3, linestyle=:dash)
xlabel!(L"d_\mathrm{E} (U, \widetilde{U})")
ylabel!(L"d_\mathrm{c} (U, \widetilde{U})")
script_dir = @__DIR__
path = joinpath(script_dir, "./figures/figure_2.pdf")
display(P)
savefig(P, path)