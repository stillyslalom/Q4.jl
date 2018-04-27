using Plots, LaTeXStrings
import Plots.pdf
cd(@__DIR__)
include("p4.jl")

f(x, y) = cos(π*x) - cos(π*y)
res = [4, 8, 16, 32]

function runtrials(res)
    p, err, t = [], Float64[], Float64[]
    for r in res
        xx = linspace(0, 1, r + 1)
        yy = linspace(0, 1, r ÷ 2)
        tic()
        A, b = Q4(f, xx, yy); sol = A\b
        push!(t, toc())

        u = reshape(A\b, (r+1, r÷2))
        push!(p, heatmap(xx, yy, u', aspect_ratio=:equal,
              xlab=L"x", ylab=L"y", title="$(r^2) elements", legendtitle=L"u(x,y)"))

        push!(err, norm(u - f.(xx, yy')/(1 + pi^2))/norm(u))
    end
    return p, err, t
end

@time p, err, t = runtrials(res)
p[4]
plot(p[1], p[2], p[3], p[4], layout = @layout([a b; c d]))
pdf("p4_soln")

plot(res, [err (res).^-2.0], xaxis=:log2, yaxis=:log2, legend=false,
    xlab=L"\mathrm{N_{elements}}", ylab=L"\mathrm{Error}")
pdf("p4error")

r = 4
xx = linspace(0, 1, r + 1)
yy = linspace(0, 1, r ÷ 2)
A, b = Q4(f, xx, yy); sol = A\b
u = reshape(A\b, (r+1, r÷2))

f.(xx,yy')
