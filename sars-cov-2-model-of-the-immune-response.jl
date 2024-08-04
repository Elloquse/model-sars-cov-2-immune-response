### A Pluto.jl notebook ###
# v0.19.42

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ a8a2521e-afb5-4e7b-827f-44ca10f7e29c
import Pkg; Pkg.add("PlotlyBase")

# ╔═╡ 992d5953-826c-42de-93f7-7346d341a52a
using DifferentialEquations, Plots, PlutoUI; plotly()

# ╔═╡ 9fbd8b12-aebe-420d-8bb1-9bba1722a79f
function f!(du,u,p,t)
    k_prod, k_v_degr, k_v_ab, omega_3, mu_ep, ep0, k_infect, k_incub, k_v_ctl, mu_epi,
    mu_idc, idc0, idc_prol, c_idc_prol, omega_1, c_vnp_mdc, mu_mdc, k_ctl_act, c_ctl,
    beta_10, mdc_to_ctl, mu_tn, tn0, mu_ctl, omega_2, k_b_act, c_b, mu_b, bn0, k_p_act,
    c_plasma, mu_p, k_ab_act, mu_ab = p
    du[1] = k_prod * u[4] - k_v_degr * u[1] - k_v_ab * u[12] * u[1] - omega_3 * u[1]
    du[2] = mu_ep * ep0 - mu_ep * u[2] - k_infect * u[1] * u[2]
    du[3] = k_infect * u[1] * u[2] - k_incub * u[3]
    du[4] = k_incub * u[3] - k_v_ctl * u[4] * u[9] - mu_epi * u[4]
    du[5] = mu_idc * idc0 - mu_idc * u[5] + idc_prol * u[5] * (u[1]/(u[1] + c_idc_prol)) - omega_1 * u[5] * (u[1]/(u[1] + c_vnp_mdc))
    du[6] = omega_1 * u[5] * (u[1]/(u[1] + c_vnp_mdc)) - mu_mdc * u[6]
    du[7] = k_ctl_act * u[7] * (u[6]/(u[6] + c_ctl)) - beta_10 * u[7] * (u[6]/(u[6] + mdc_to_ctl)) - mu_tn * u[7] + mu_tn * tn0
    du[8] = beta_10 * u[7] * (u[6]/(u[6] + mdc_to_ctl)) - mu_ctl * u[8]
    du[9] = omega_2 * u[8] - mu_ctl * u[9]
    du[10] = k_b_act * u[10] * (u[6]/(u[6] + c_b)) - mu_b * u[10] + mu_b * bn0 - k_p_act * u[10] * (u[6]/(u[6] + c_plasma))
    du[11] = k_p_act * u[10] * (u[6]/(u[6] + c_plasma)) - mu_p * u[11]
    du[12] = k_ab_act * u[11] - mu_ab * u[12]
    nothing
end

# ╔═╡ dfae18d3-bed7-42a2-828d-c693eebb61b1
begin
	u0 = [1000.0, 5500000.0, 0.0, 0.0, 2500000.0, 0.0, 16000.0, 0.0, 0.0, 33000.0, 0.0, 0.0]
	tspan = (0.0, 30.0)
	p = (10.0, 0.5, 0.5, 0.1, 0.8, 5500000.0, 1.5e-7, 1.65, 0.001, 0.5, 0.13, 2500000.0, 0.15, 1.0e7, 0.3, 5000000.0, 0.5, 15.0, 1.0e7, 22.0, 1.0e7, 0.008, 16000.0, 0.05, 0.015, 5.5, 1.0e7, 0.008, 33000.0, 6.0, 1.0e7, 0.2, 3.0e-5, 0.025)
end

# ╔═╡ d05b1b68-2e87-49b1-9b6d-f8231181121d
md"## ODE solution"

# ╔═╡ f32d3960-f0d3-4360-9dfd-be946a4e8abc
prob_ODE = ODEProblem(f!, u0, tspan, p)

# ╔═╡ 91cf046b-15b2-4b67-bcff-df71f100a6c9
sol_ODE = solve(prob_ODE, Rodas5(), abstol=1.0E-20, reltol=1.0E-12, dt=0.5)

# ╔═╡ 4d523a95-a671-4708-9a91-0acaa0e7855c
md" ## SDE solution"

# ╔═╡ 1019f8ed-906b-4921-ad31-553c1ffe5deb
@bind new_sol Button("New Solution")

# ╔═╡ 55afa6ed-d67a-4931-acf5-0c8e6e07f6ce
md"## Plots"

# ╔═╡ d607c2b0-86aa-4606-a62c-70af763f4b2e
@bind new_plot Button("New plot (SDE)")

# ╔═╡ 7075d6fc-fb3b-403c-ab0e-7dcc454f8a24
begin
	new_plot
	p_ODE = plot(sol_ODE.t, sol_ODE[1,:], size=(600,500), legend=false, linewidth=3, title="Solution (ODE + Noize)")
end

# ╔═╡ 40c2790d-f077-4913-8b82-44ca72568a00
@bind noize Slider(0.0:0.1:3.0, 1.0, true)

# ╔═╡ 1a332714-8c3b-4ff6-9795-8fcabf5483fb
function g!(du,u,p,t)
	
    du[1] = noize * u[1]
    du[2] = 0
    du[3] = 0
    du[4] = 0
    du[5] = 0
    du[6] = 0
    du[7] = 0
    du[8] = 0
    du[9] = 0
    du[10] = 0
    du[11] = 0
    du[12] = 0
    nothing
end

# ╔═╡ 0a329bf5-319f-42de-b478-4c04df179f37
prob_sde = SDEProblem(f!, g!, u0, tspan, p)

# ╔═╡ 011c1d3b-bb7a-4b56-aaee-29901338fb90
begin
	new_sol
	sol_SDE = solve(prob_sde, EM(), dt=0.005)
end

# ╔═╡ 1ac87c3c-f3dd-4099-9839-ff006a57871d
plot!(p_ODE, sol_SDE, noize=0.5)

# ╔═╡ Cell order:
# ╠═a8a2521e-afb5-4e7b-827f-44ca10f7e29c
# ╠═992d5953-826c-42de-93f7-7346d341a52a
# ╠═9fbd8b12-aebe-420d-8bb1-9bba1722a79f
# ╠═1a332714-8c3b-4ff6-9795-8fcabf5483fb
# ╠═dfae18d3-bed7-42a2-828d-c693eebb61b1
# ╟─d05b1b68-2e87-49b1-9b6d-f8231181121d
# ╠═f32d3960-f0d3-4360-9dfd-be946a4e8abc
# ╠═91cf046b-15b2-4b67-bcff-df71f100a6c9
# ╟─4d523a95-a671-4708-9a91-0acaa0e7855c
# ╟─1019f8ed-906b-4921-ad31-553c1ffe5deb
# ╠═0a329bf5-319f-42de-b478-4c04df179f37
# ╠═011c1d3b-bb7a-4b56-aaee-29901338fb90
# ╟─55afa6ed-d67a-4931-acf5-0c8e6e07f6ce
# ╠═7075d6fc-fb3b-403c-ab0e-7dcc454f8a24
# ╟─d607c2b0-86aa-4606-a62c-70af763f4b2e
# ╠═1ac87c3c-f3dd-4099-9839-ff006a57871d
# ╠═40c2790d-f077-4913-8b82-44ca72568a00
