using LinearAlgebra,SkewLinearAlgebra, Plots, LaTeXStrings
import .Manifold as m

#Choose (n,p) the size of the Stiefel manifold, N the number of samples
n = 5
p = 4
N = 20000

metric = 0

U₁ = m.Orthonormal(Matrix(1.0*I,n,p), false)
Θ = m.Orthonormal(Matrix(1.0*I,n,n)[:,p+1:end], false)
S₁= m.StiefelVector(U₁, metric)
S₁.Θ = Θ


A = [0  1 -1;
    -1  0  1;
     1 -1  0] * π / sqrt(3)
B = [1  1  1] * π / sqrt(3)

S = [2A -B';
      B  0]
N2 = 20
norms = zeros(N2, p)
distances = zeros(N2, p)
for k ∈ 1:p
    U₂ = copy(U₁)
    if k == 1
        for i ∈ 1:N2
            θ = (i / N2) * π
            U₂.Q[:,1] = cos(θ)*U₁.Q[:,1] + sin(θ)*S₁.Θ.Q[:,1]
            norms[i, k] = norm(U₂.Q-U₁.Q)
            distances[i, k] = θ
        end
    elseif iseven(k)
        Q = zeros(k, k)
        #Build Q
        for i ∈ 1:N2
            θ = (i / N2) * π
            for j ∈ 1:2:k
                Q[j:j+1, j:j+1] = [cos(θ) -sin(θ); sin(θ) cos(θ)]
            end
            
            U₂.Q[:,1:k] = U₁.Q[:,1:k] * Q
            norms[i, k] = norm(U₂.Q - U₁.Q)
            distances[i, k] = θ *sqrt(k)
        end
    else
        M = zeros(k+1, k+1)
        #Build Q
        for i ∈ 1:N2
            for j ∈ 1:2:(k-3)
                M[j:j+1, j:j+1] = [0 -π; π 0]
            end
            M[k-2:end, k-2:end] = S
            U₂.Q[:, 1:k] = ([U₁.Q[:,1:k] S₁.Θ.Q[:,1]] * exp((i / N2) * M))[:, 1:k]
            U₂.Q[:, k-2:k] = U₂.Q[:, k-2:k] * exp(- (i / N2) * A)
            norms[i, k] = norm(U₂.Q-U₁.Q)
            distances[i, k] = (i / N2) * π *sqrt(k)
        end
    end
end
norms ./= 2 * sqrt(p)
distances  ./= π * sqrt(p)


Frobenius = zeros(N)
gdistances = zeros(N)
S = zeros(2p,2p)
for i ∈ 1:N
    U = m.randorthonormal(Float64, n, n)
    U₃ = m.Orthonormal(U.Q[:,1:p], false)
    coin = rand((-1,1))
    U₂ =  m.orthogonalize!((coin*U.Q*exp(skewhermitian!(10*(i/N) *randn(n,n))))[:,1:p])
    
    Frobenius[i] = norm(Matrix(U₃ - U₂))
    S₃= m.StiefelVector(U₃, 0)
    S₂= m.StiefelVector(U₂, 0)
    Δ, Q = m.pshooting(S₃, S₂, 1, 10)[1:2]
    gdistances[i] = norm(Δ)
    #=
    if 2 < Frobenius[i] < 2√2 && gdistances[i] >  sqrt((π^2) + 4 * asin(sqrt(Frobenius[i]^2 - 4)/2)^2)*1.01
        print("A is ")
        A = U₃.Q'Δ
        display(A)
        print("B is ")
        display(norm(Q'Δ))
    end
    =#
end

t = Array(LinRange(0, 2*sqrt(p),500))
ub = π/2 .* t
jj = findmax(findall(x -> x < 2, t))[1]
ub[1:jj] = 2 * asin.(t[1:jj] ./2)
P = plot(framestyle=:box, legend=:topleft,font="Computer Modern", tickfontfamily="Computer Modern",legendfont="Computer Modern",
    legendfontsize=14,yguidefontsize=17,xguidefontsize=17, xtickfontsize = 13, ytickfontsize=13,margin = 0.3Plots.cm, dpi = 1000)

scatter!(Frobenius./(2*sqrt(p)), min.(gdistances,π*sqrt(p))./(π*sqrt(p)), markersize=3, markerstrokewidth=0,label = "Samples", color=:skyblue2)
plot!(t./(2*sqrt(p)), ub ./(π*sqrt(p)),linewidth=2, label = L"\widehat{M}_\mathrm{E}", linestyle=:dot, color=:darkgreen)

plot!(norms[:, 1], distances[:, 1], color = :black, label = "Minimal geodesics", linewidth = 2, linestyle=:solid)
for i ∈ 2:p
    plot!(norms[:, i], distances[:, i], color = :black, label = false, linewidth = 2, linestyle=:solid)
end
#=
for k ∈ 1:p-1
    ε = 1e-5
    δ = LinRange(2√k + ε, 2√(k+1) - ε, 100)
    u = sqrt.((k * π^2) .+ 4 .* asin.(sqrt.(δ.^2 .- 4k)./2).^2)
    plot!(δ ./ (2√p), u ./ (π√p), color =:purple, linweidth = 2, label = false)
end
=#
xlabel!(L"\Vert \widetilde{U}-U \ \Vert_\mathrm{F} / 2\sqrt{p}")
ylabel!(L"d_\mathrm{E} (U, \widetilde{U})/\pi \sqrt{p}")
script_dir = @__DIR__
path = joinpath(script_dir, "./figures/figure_10.pdf")
savefig(P, path)
path2 = joinpath(script_dir, "./figures/figure_10.png")
savefig(P, path2)
