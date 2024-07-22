using LinearAlgebra, SkewLinearAlgebra, Plots, LaTeXStrings
import .Manifold as m

N = 10000
n = 3
p = 2

Frobenius = zeros(N)
distances = zeros(N)
β = 2
count = 1
U  = m.Orthonormal(Matrix(1.0I,n, n), false)
U₁ = m.Orthonormal(Matrix(1.0I,n, p), false)
S₁= m.StiefelVector(U₁, β)
for i ∈ 1:N-200
    global count, U, U₁, S₁
    #m.buildcomplement(S₁)
    if iseven(count)
        U₂ =  m.orthogonalize!((U.Q*exp(skewhermitian!(2 * (i / N) * randn(n,n))))[:,1:p])   
        count +=1
    else
        U₂ =  m.orthogonalize!((-U.Q*exp(skewhermitian!(3 * (i / N) * randn(n,n))))[:,1:p])
        count += 1
    end
    Frobenius[i] = norm(Matrix(U₁ - U₂))
    if Frobenius[i] > 0
        S₂= m.StiefelVector(U₂, β)
        #Δ, Q = m.logβ(S₁, S₂, β, 2)[1:2]
        Δ, Q = m.pshooting(S₁, S₂, β, 10)[1:2]
        A = U₁'*Δ
        B = Q'*Δ
        
        
        if norm(exp(S₁, Δ) - S₂.U.Q) < 0.01 && sqrt(β * norm(A)^2 + norm(B)^2) ≤ max(1,sqrt(β)) * π * sqrt(p)
            distances[i] = sqrt(β * norm(A)^2 + norm(B)^2)
        end
    end
end

for i ∈ 1:200
    global  U, U₁, S₁
    j = i+ N-200
    U₂ =  m.orthogonalize!(U₁.Q * exp(skewhermitian!(π * (i / 200) *randn(p,p))))
    Frobenius[j] = norm(Matrix(U₁ - U₂))
    S₂= m.StiefelVector(U₂, β)
    Δ, Q = m.pshooting(S₁, S₂, β, 10)[1:2]
    A = U₁'*Δ
    B = Q'*Δ
    if norm(exp(S₁, Δ) - S₂.U.Q) < 0.01 && sqrt(β * norm(A)^2 + norm(B)^2) ≤ max(1,sqrt(β)) * π * sqrt(p)
        distances[j] = sqrt(β * norm(A)^2 + norm(B)^2)
    end
end

for i ∈1:N
    if iszero(distances[i])
        Frobenius[i] = 0
    end
end

P = plot(framestyle=:box, legend=:topleft,font="Computer Modern", tickfontfamily="Computer Modern",legendfont="Computer Modern",
    legendfontsize=14,yguidefontsize=17,xguidefontsize=17, xtickfontsize = 13, ytickfontsize=13,margin = 0.3Plots.cm, dpi = 1000)

t = LinRange(0, 2 * sqrt(p), 500)
t2 = Array(LinRange(0, 2 , 400))
t3 = LinRange(2, 2 * sqrt(p), 100)
scatter!(Frobenius./(2*sqrt(p)), distances./(π*sqrt(β*p)), markersize=2, markerstrokewidth=0,label = "Samples", color=:skyblue2)
plot!(t./(2*sqrt(p)),2*min(1,sqrt(β))*sqrt(p)*asin.(t./(2*sqrt(p)))./(π*sqrt(p*β)), linewidth=2, label = L"\widehat{m}_\beta",color=:black)
plot!(t2./(2*sqrt(p)), max(1,sqrt(β))* 2 .* asin.(t2./2)./(π*sqrt(p*β)) , linewidth=2, label = false,color=:darkgreen, linestyle=:dot)
plot!(t3./(2*sqrt(p)), max(1,sqrt(β))* π/2 .* t3./(π*sqrt(p*β)) , linewidth=2, label = L"\widehat{M}_\beta",color=:darkgreen, linestyle=:dot)
xlabel!(L"\Vert \widetilde{U}-U \ \Vert_\mathrm{F} / 2\sqrt{p}")
ylabel!(L"d_\mathrm{\beta} (U, \widetilde{U})/\pi \sqrt{p}")

script_dir = @__DIR__
path = joinpath(script_dir, "./figures/figure_6bis.pdf")
savefig(P, path)
path = joinpath(script_dir, "./figures/figure_6bis.png")
savefig(P, path)


