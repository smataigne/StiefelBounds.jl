using LinearAlgebra, SkewLinearAlgebra, Plots, LaTeXStrings
import .Manifold as m

#Parameters

n = 6
p = n - 1
N = 500000
β = 0.5    #keep β = 0.5 (canonical metric)

U₁ = zeros(n, n)
U₂ = similar(U₁, n, n)
temp = similar(U₁, n, n)

Frobenius = zeros(N)
distances = zeros(N)
δeven = 0
δodd = 0
for i ∈ 1:(7N÷10)
    global  δeven, δodd
    U₁ .= exp(skewhermitian!(randn(n , n)))
    coin = (rand() > 0.8 ? -1 : 1)
    temp .= U₁
    if iseven(p)
        temp[:, 1:p] .*= coin
    else
        temp .*= coin
    end
    if coin == 1
        δeven += 1/N * 10/7
        mul!(U₂ , temp, exp(skewhermitian!(randn(n, n) * 2 * δeven)), 1, 0)
    else
        δodd += 1/N *  10/7
        mul!(U₂ , temp, exp(skewhermitian!(randn(n, n) * 2 * δodd)), 1, 0)
    end
    Frobenius[i] = norm(U₁[:, 1:p] - U₂[:, 1:p])
    distances[i] = norm(m.skewlog(U₁'U₂)) * sqrt(2) / 2
end
δ = 0
#Sample first corner
for i ∈ (8N÷10)+1:(9N÷10)
    global  δ
    U₁ .= exp(skewhermitian!(randn(n , n)))
    temp .= U₁
    θ = rand() * π
    v1 = copy(temp[:, 1]); vn = copy(temp[:, n] )
    @. temp[:, 1] = cos(θ) *  v1 + sin(θ) * vn
    @. temp[:, n] = -sin(θ) *  v1 + cos(θ) * vn
    δ += 1/N
    mul!(U₂ , temp, exp(skewhermitian!(randn(n, n) * 10 * δ)), 1, 0)
    Frobenius[i] = norm(U₁[:, 1:p] - U₂[:, 1:p])
    distances[i] = norm(m.skewlog(U₁'U₂)) * sqrt(2) / 2
end

#Sample second corner
δ = 0
for i ∈ (9N÷10)+1:N
    global  δ
    U₁ .= exp(skewhermitian!(randn(n , n)))
    temp .= U₁
    θ = rand() * π
    temp[:, 1] .*= -1
    v1 = copy(temp[:, 2]); vn = copy(temp[:, n] )
    @. temp[:, 2] = cos(θ) *  v1 + sin(θ) * vn
    @. temp[:, n] = sin(θ) *  v1 - cos(θ) * vn
    θ = rand() * π
    v1 = copy(temp[:, 3]); vn = copy(temp[:, n-1] )
    @. temp[:, 3] = cos(θ) *  v1 + sin(θ) * vn
    @. temp[:, n-1] = -sin(θ) *  v1 + cos(θ) * vn
    δ += 1/N
    mul!(U₂ , temp, exp(skewhermitian!(randn(n, n)  * δ)), 1, 0)
    Frobenius[i] = norm(U₁[:, 1:p] - U₂[:, 1:p])
    distances[i] = norm(m.skewlog(U₁'U₂)) * sqrt(2) / 2
end

#third corner sampling
for i ∈ (7N÷10)+1:(8N÷10)
    global  δ
    U₁ .= exp(skewhermitian!(randn(n , n)))
    temp .= U₁
    temp[:,1:3].*=-1
    θ = rand() * π
    v1 = copy(temp[:, 4]); vn = copy(temp[:, n] )
    @. temp[:, 4] = cos(θ) *  v1 + sin(θ) * vn
    @. temp[:, n] = sin(θ) *  v1 - cos(θ) * vn
    θ = rand() * π
    v1 = copy(temp[:, 5]); vn = copy(temp[:, n] )
    @. temp[:, 5] = cos(θ) *  v1 + sin(θ) * vn
    @. temp[:, n] = -sin(θ) *  v1 + cos(θ) * vn
    δ += 1/N
    mul!(U₂ , temp, exp(skewhermitian!(randn(n, n)  * δ)), 1, 0)
    Frobenius[i] = norm(U₁[:, 1:p] - U₂[:, 1:p])
    distances[i] = norm(m.skewlog(U₁'U₂)) * sqrt(2) / 2
end



t = LinRange(0, 2*sqrt(p)-0.01,500)
t2 = LinRange(0, 2*sqrt(p-1),500)
t3 = LinRange(2*sqrt(p-1), 2*sqrt(p),100)
P = plot(framestyle=:box, legend=:topleft,font="Computer Modern", tickfontfamily="Computer Modern",legendfont="Computer Modern",
    legendfontsize=14,yguidefontsize=17,xguidefontsize=17, xtickfontsize = 13, ytickfontsize=13,margin = 0.3Plots.cm, dpi = 1000)
t4 = Array(LinRange(0, 2 , 400))
t5 = LinRange(2, 2 * sqrt(p), 100)
scatter!(Frobenius./(2*sqrt(p)), min.(distances,π*sqrt(p))./(π*sqrt(p / 2)), markersize=2, markerstrokewidth=0,label = "Samples", color=:skyblue2)
plot!(t./(2*sqrt(p)), sqrt(2p)*asin.(t./(2*sqrt(p)))./(π*sqrt(p / 2)), linewidth=2, label = L"\widehat{m}_\mathrm{c}",color=:black)
plot!(t4./(2*sqrt(p)), 2 .* asin.(t4./2)./(π*sqrt(p*β)) , linewidth=2, label = false,color=:darkgreen, linestyle=:dot)
plot!(t5./(2*sqrt(p)), π/2 .* t5./(π*sqrt(p*β)) , linewidth=2, label = L"\widehat{M}_\mathrm{c}",color=:darkgreen, linestyle=:dot)
#Julia legend doesn't displays correctly
plot!(t2[1:2]./(2*sqrt(p)), min.(sqrt(2*(p-1))*asin.(t2[1:2]./(2*sqrt(p-1)))./(π*sqrt(p / 2)), sqrt(2(p+1))*asin.(t2[1:2]./(2*sqrt(p)))./(π*sqrt(p / 2))), label = "Thm. 5.8", color =:purple, linewidth=1, linestyle=:dashdot)
plot!(t2./(2*sqrt(p)), min.(sqrt(2*(p-1))*asin.(t2./(2*sqrt(p-1)))./(π*sqrt(p / 2)), sqrt(2(p+1))*asin.(t2./(2*sqrt(p)))./(π*sqrt(p / 2))), label = false, color =:purple, linewidth=2, linestyle=:dashdot)
plot!(t3./(2*sqrt(p)),  sqrt(2(p+1))*asin.(t3./(2*sqrt(p)))./(π*sqrt(p / 2)), label = false, color =:purple, linewidth=2, linestyle=:dashdot)
xlabel!(L"\Vert \widetilde{U}-U \ \Vert_\mathrm{F} / 2\sqrt{p}")
ylabel!(L"\sqrt{2} d_\mathrm{c} (U, \widetilde{U})/\pi \sqrt{p}")
script_dir = @__DIR__
#path = joinpath(script_dir, "./figures/figure_6.pdf")
#savefig(P, path)
path2 = joinpath(script_dir, "./figures/figure_6.png")
savefig(P, path2)


