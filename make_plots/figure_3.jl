using LinearAlgebra, SkewLinearAlgebra, Plots, LaTeXStrings
import .Manifold as m

#Choose N, the number of samples, and (n,p) the size of the Stiefel manifold
N = 20000
n = 5
p = 4

Frobenius = zeros(N)
distances = zeros(N)
S = zeros(2p,2p)
for i ∈ 1:N
    U = m.randorthonormal(Float64, n, n)
    U₁ = m.Orthonormal(U.Q[:,1:p], false)
    U₂ =  m.orthogonalize!((U.Q * exp(skewhermitian!(10 * (i/N) * randn(n,n))))[:,1:p])
    
    Frobenius[i] = norm(Matrix(U₁ - U₂))
    S₁= m.StiefelVector(U₁, 0)
    S₂= m.StiefelVector(U₂, 0)
    m.buildcomplement(S₁)
    m.buildcomplement(S₂)
    Δ = m.logshooting(S₁, S₂)
    A = U₁'*Δ
    B = S₁.Θ.Q'*Δ
    distances[i] = sqrt(norm(A)^2+norm(B)^2)
end

P = plot(framestyle=:box, legend=:topleft,font="Computer Modern", tickfontfamily="Computer Modern",legendfont="Computer Modern",
    legendfontsize=12,yguidefontsize=17,xguidefontsize=17, xtickfontsize = 13, ytickfontsize=13,margin = 0.3Plots.cm)

t = LinRange(0, 2*sqrt(p),500)

scatter!(Frobenius./(2*sqrt(p)), min.(distances,π*sqrt(p))./(π*sqrt(p)), markersize=3, markerstrokewidth=0,label = "Samples")
plot!(t./(2*sqrt(p)),2*sqrt(p)*asin.(t./(2*sqrt(p)))./(π*sqrt(p)), linewidth=3, label = "Lower bound",color=:black)
plot!(t./(2*sqrt(p)), t.*π/2 ./(π*sqrt(p)),linewidth=3, label = "Upper bound", linestyle=:dash, color=:green3)
xlabel!(L"\Vert \widetilde{U}-U \ \Vert_\mathrm{F} / 2\sqrt{p}")
ylabel!(L"d_\mathrm{E} (U, \widetilde{U})/\pi \sqrt{p}")
script_dir = @__DIR__
path = joinpath(script_dir, "./figures/figure_3.pdf")
display(P)
savefig(P, path)


