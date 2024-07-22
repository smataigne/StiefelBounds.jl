using LinearAlgebra, SkewLinearAlgebra, Plots, LaTeXStrings
import .Manifold as m

#Choose N, the number of samples, and (n,p) the size of the Stiefel manifold
N = 1000000
n = 8
p = 4

Frobenius = zeros(N)
distances = zeros(N)
S = zeros(2p,2p)
a = 9 * N ÷ 10
for i ∈ 1:N
    U  = m.randorthonormal(Float64, n, n)
    U₁ = m.Orthonormal(U.Q[:,1:p], false)
    if i ≤ a
        U₂ = m.orthogonalize!((U.Q * exp(skewhermitian!(10 * (i/N) * randn(n,n))))[:,1:p])
    else
        U₂ = m.orthogonalize!((U.Q[:,1:p] * exp(skewhermitian!(10 * ((i-a)/(N-a)) * randn(p, p)))))
    end
    Frobenius[i] = norm(Matrix(U₁ - U₂))
    S₁= m.StiefelVector(U₁, 0)
    S₂= m.StiefelVector(U₂, 0)
    Δ = m.logβ(S₁, S₂, 0.5, 1)[1]
    distances[i] = m.normβ(Δ, S₁, 0.5)
end

P = plot(framestyle=:box, legend=:topleft,font="Computer Modern", tickfontfamily="Computer Modern",legendfont="Computer Modern",
    legendfontsize=12,yguidefontsize=17,xguidefontsize=17, xtickfontsize = 13, ytickfontsize=13,margin = 0.3Plots.cm, dpi = 1000)

t = Array(LinRange(0, 2*sqrt(p), 500))

ub = π/2 .* t
ub[1:250] = 2 * asin.(t[1:250] ./2)

scatter!(Frobenius./(2*sqrt(p)), min.(distances,π*sqrt(p))./(π*sqrt(p)), markersize=2, markerstrokewidth=0,label = "Samples", color=:skyblue2)
plot!(t./(2*sqrt(p)),sqrt(2*p)*asin.(t./(2*sqrt(p)))./(π*sqrt(p)), linewidth=2, label = L"\widehat{m}_\mathrm{c}",color=:black)
plot!(t./(2*sqrt(p)), ub ./(π*sqrt(p)), linewidth=2, label = L"\widehat{M}_\mathrm{c}", linestyle=:dot, color=:darkgreen)
xlabel!(L"\Vert \widetilde{U}-U \ \Vert_\mathrm{F} / 2\sqrt{p}")
ylabel!(L"d_\mathrm{c} (U, \widetilde{U})/\pi \sqrt{p}")
script_dir = @__DIR__
#path = joinpath(script_dir, "./figures/figure_3_bis.pdf")
#savefig(P, path)
path2 = joinpath(script_dir, "./figures/figure_3_bis.png")
savefig(P, path2)