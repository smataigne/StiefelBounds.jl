using LinearAlgebra, SkewLinearAlgebra, Plots, LaTeXStrings
import .Manifold as m

#Choose N, the number of samples, and (n,p) the size of the Stiefel manifold
N = 20000
n = 5
p = 4

Frobenius = zeros(N)
distances = zeros(N)
S = zeros(2p,2p)
a = N - 1000
for i ∈ 1:N
    U  = m.randorthonormal(Float64, n, n)
    U₁ = m.Orthonormal(U.Q[:,1:p], false)
    if i ≤ a
        U₂ = m.orthogonalize!((U.Q * exp(skewhermitian!(10 * (i/N) * randn(n,n))))[:,1:p])
    else
        U₂ = m.orthogonalize!((-U.Q * exp(skewhermitian!(10 * ((i-a)/(N-a)) * randn(n,n))))[:,1:p])
    end
    Frobenius[i] = norm(Matrix(U₁ - U₂))
    S₁= m.StiefelVector(U₁, 0)
    S₂= m.StiefelVector(U₂, 0)
    Δ = m.pshooting(S₁, S₂, 1, 10)[1]
    distances[i] = norm(Δ)
end

P = plot(framestyle=:box, legend=:topleft,font="Computer Modern", tickfontfamily="Computer Modern",legendfont="Computer Modern",
    legendfontsize=12,yguidefontsize=17,xguidefontsize=17, xtickfontsize = 13, ytickfontsize=13,margin = 0.3Plots.cm)

t = Array(LinRange(0, 2*sqrt(p), 500))
ub = π/2 .* t
ub[1:250] = 2 * asin.(t[1:250] ./2)
scatter!(Frobenius./(2*sqrt(p)), min.(distances,π*sqrt(p))./(π*sqrt(p)), markersize=2, markerstrokewidth=0,label = "Samples", color=:skyblue2)
plot!(t./(2*sqrt(p)),2*sqrt(p)*asin.(t./(2*sqrt(p)))./(π*sqrt(p)), linewidth=2, label = L"\widehat{m}_\mathrm{E}",color=:black)
plot!(t./(2*sqrt(p)), ub ./(π*sqrt(p)),linewidth=2, label = L"\widehat{M}_\mathrm{E}", linestyle=:dot, color=:darkgreen)
xlabel!(L"\Vert \widetilde{U}-U \ \Vert_\mathrm{F} / 2\sqrt{p}")
ylabel!(L"d_\mathrm{E} (U, \widetilde{U})/\pi \sqrt{p}")
script_dir = @__DIR__
path = joinpath(script_dir, "./figures/figure_3.pdf")
savefig(P, path)
path2 = joinpath(script_dir, "./figures/figure_3.png")
savefig(P, path2)


