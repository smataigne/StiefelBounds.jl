using LinearAlgebra, SkewLinearAlgebra, Plots, LaTeXStrings
import .Manifold as m

#Choose N, the number of samples, and (n,p) the size of the Stiefel manifold
N = 100000
n = 4
p = 3

Nsec = 200
Frobenius = zeros(Nsec)
distances1 = ones(Nsec)
distances2 = ones(Nsec)
S = zeros(2p,2p)
for i ∈ 1:N
    U  = m.randorthonormal(Float64, n, n)
    U₁ = m.Orthonormal(U.Q[:,1:p], false)
    coin = rand((-1,1))
    U₂ = m.orthogonalize!((coin * U.Q * exp(skewhermitian!(2 * (i/N) * randn(n,n))))[:,1:p])
    
    s = norm(Matrix(U₁ - U₂))
    index = Int(floor(s / (2*sqrt(p)) * Nsec)) + 1
    S₁= m.StiefelVector(U₁, 0)
    S₂= m.StiefelVector(U₂, 0)
    Δ = m.pshooting(S₁, S₂, 1, 10)[1]
    d1 = max((s * π/2 - norm(Δ))/norm(Δ), 0)
    d2 = (s > 2 ? max((s * π/2 - norm(Δ))/norm(Δ), 0) : max((2 * asin(s/2) - norm(Δ))/norm(Δ), 0))
    distances1[index] = (d1 < distances1[index] ? d1 : distances1[index])
    distances2[index] = (d2 < distances2[index] ? d2 : distances2[index])
end

P = plot(framestyle=:box, legend=:topright,font="Computer Modern", tickfontfamily="Computer Modern", guidefontfamily="Computer Modern", legendfont="Computer Modern",
    legendfontsize=12,yguidefontsize=17,xguidefontsize=17, xtickfontsize = 13, ytickfontsize=13,margin = 0.3Plots.cm, dpi = 1000)

t = LinRange(0, 2*sqrt(p), 500)

plot!(LinRange(0, 1, Nsec), distances1, label = "Thm 6.3", color=:skyblue2, linewidth = 2)
plot!(LinRange(0, 1, Nsec), distances2, label = "Thms 6.3 + 6.4", color=:darkgreen, linewidth = 2, linestyle=:dot)
xlabel!(L"\delta / 2\sqrt{p}")
ylabel!("Relative Gap "*L"\omega(\delta)")
display(P)

script_dir = @__DIR__
path = joinpath(script_dir, "./figures/figure_7.pdf")
savefig(P, path)
path2 = joinpath(script_dir, "./figures/figure_7.png")
savefig(P, path2)