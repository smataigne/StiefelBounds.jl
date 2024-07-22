using LinearAlgebra, SkewLinearAlgebra, Plots, LaTeXStrings

n = 4
p = 2
N = 10000
l = 15
βs = zeros(N)

f(t, Aaug, S) = sqrt(2p - 2 * tr((SkewLinearAlgebra.hermitian!(exp(t .* S) * exp(t .* Aaug)))[1:p, 1:p]))

r = zeros(N)

for i ∈ 1:N
    β = rand() * 1.4 + 0.1
    βs[i] = β
    A = skewhermitian!(rand(p, p))*10
    B = randn(n - p, p)

    nm = sqrt(β * norm(A)^2 + norm(B)^2)

    A /= nm
    B ./= nm
    S = skewhermitian([2β*A -B';B zeros(n - p, n - p)])
    Aaug = skewhermitian([(1 - 2β)*A zeros(p, n-p);zeros(n - p, n)])

    δ = 1e-6
    tmax = max(1,√β) * δ
    tmin = min(1, √β) * δ
    t = (tmax + tmin) / 2
    iter = 0
    δ_true = f(t, Aaug, S)
    while abs(δ_true - δ) > 1e-9 && iter < 50
        iter += 1
        if δ_true > δ
            tmax = t
        else
            tmin = t
        end
        t = (tmax + tmin) / 2
        δ_true = f(t, Aaug, S)
    end
    r[i] = t / δ
end

P = plot(framestyle=:box, legend=:bottomright,font="Computer Modern", tickfontfamily="Computer Modern",legendfont="Computer Modern",
    legendfontsize=12,yguidefontsize=17,xguidefontsize=17, xtickfontsize = 13, ytickfontsize=13,margin = 0.3Plots.cm, dpi = 1000)

t = LinRange(0.1, 1.5, 100)    
scatter!(βs, r , label = "Samples", markersize = 2, color=:skyblue2, markerstrokewidth = 0)
plot!(t, min.(ones(100), sqrt.(t)), label = L"\min \{1, \sqrt{\beta}\}", linestyle=:solid, color =:black, linewidth = 2)
plot!(t, max.(ones(100), sqrt.(t)), label = L"\max \{1, \sqrt{\beta}\}", linestyle=:dot, color=:darkgreen, linewidth = 2)
xlabel!(L"$\beta$")
ylabel!(L"$\lim_{\delta\rightarrow 0} r_{\beta,U,\Delta} (\delta)$")
script_dir = @__DIR__
path = joinpath(script_dir, "./figures/figure_11.pdf")
savefig(P, path)
path = joinpath(script_dir, "./figures/figure_11.png")
savefig(P, path)
