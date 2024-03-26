using LinearAlgebra,SkewLinearAlgebra, Plots, LaTeXStrings
import .Manifold as m

#Choose (n,p) the size of the Stiefel manifold
n = 5
p = 4

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
N = 20
norms = zeros(N, p)
distances = zeros(N, p)
for k ∈ 1:p
    U₂ = copy(U₁)
    if k == 1
        for i ∈ 1:N
            θ = (i / N) * π
            U₂.Q[:,1] = cos(θ)*U₁.Q[:,1] + sin(θ)*S₁.Θ.Q[:,1]
            norms[i, k] = norm(U₂.Q-U₁.Q)
            distances[i, k] = θ
        end
    elseif iseven(k)
        Q = zeros(k, k)
        #Build Q
        for i ∈ 1:N
            θ = (i / N) * π
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
        for i ∈ 1:N
            for j ∈ 1:2:(k-3)
                M[j:j+1, j:j+1] = [0 -π; π 0]
            end
            M[k-2:end, k-2:end] = S
            U₂.Q[:, 1:k] = ([U₁.Q[:,1:k] S₁.Θ.Q[:,1]] * exp((i / N) * M))[:, 1:k]
            U₂.Q[:, k-2:k] = U₂.Q[:, k-2:k] * exp(- (i / N) * A)
            norms[i, k] = norm(U₂.Q-U₁.Q)
            distances[i, k] = (i / N) * π *sqrt(k)
        end
    end
end
norms ./= 2 * sqrt(p)
distances  ./= π * sqrt(p)

#Choose N, the number of samples.
N = 10000
Frobenius = zeros(N)
gdistances = zeros(N)
S = zeros(2p,2p)
for i ∈ 1:N
    U = m.randorthonormal(Float64, n, n)
    U₃ = m.Orthonormal(U.Q[:,1:p], false)
    U₂ =  m.orthogonalize!((U.Q*exp(skewhermitian!(10*(i/N) *randn(n,n))))[:,1:p])
    
    Frobenius[i] = norm(Matrix(U₃ - U₂))
    S₃= m.StiefelVector(U₃, 0)
    S₂= m.StiefelVector(U₂, 0)
    m.buildcomplement(S₃)
    m.buildcomplement(S₂)
    Δ = m.logshooting(S₃, S₂)
    A₂ = U₁'*Δ
    B₂ = S₁.Θ.Q'*Δ
    gdistances[i] = sqrt(norm(A₂)^2+norm(B₂)^2)
end

t = LinRange(0, 2*sqrt(p),500)
P = plot(framestyle=:box, legend=:topleft,font="Computer Modern", tickfontfamily="Computer Modern",legendfont="Computer Modern",
    legendfontsize=12,yguidefontsize=17,xguidefontsize=17, xtickfontsize = 13, ytickfontsize=13,margin = 0.3Plots.cm)

scatter!(Frobenius./(2*sqrt(p)), min.(gdistances,π*sqrt(p))./(π*sqrt(p)), markersize=3, markerstrokewidth=0,label = "Samples")
plot!(norms[:, 1], distances[:, 1], color = :black, label = "Geodesics", linewidth = 3, linestyle=:solid)
for i ∈ 2:p
    plot!(norms[:, i], distances[:, i], color = :black, label = false, linewidth = 3, linestyle=:solid)
end
plot!(t./(2*sqrt(p)), t.*π/2 ./(π*sqrt(p)),linewidth=3, label = "Upper bound", linestyle=:dash, color=:green3)
xlabel!(L"\Vert \widetilde{U}-U \ \Vert_\mathrm{F} / 2\sqrt{p}")
ylabel!(L"d_\mathrm{E} (U, \widetilde{U})/\pi \sqrt{p}")
script_dir = @__DIR__
path = joinpath(script_dir, "./figures/figure_7.pdf")
display(P)
savefig(P, path)