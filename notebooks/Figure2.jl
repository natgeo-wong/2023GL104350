### A Pluto.jl notebook ###
# v0.19.26

using Markdown
using InteractiveUtils

# ╔═╡ 681658b0-5914-11eb-0d65-bbace277d145
begin
	using Pkg; Pkg.activate()
	using DrWatson
	
md"Using DrWatson in order to ensure reproducibility between different machines ..."
end

# ╔═╡ 6dce35fc-5914-11eb-0ce2-0d4e164e1898
begin
	@quickactivate "2023GL104350"
	using NCDatasets
	using Printf
	using StatsBase
	
	using ImageShow, PNGFiles
	using PyCall, LaTeXStrings
	pplt = pyimport("proplot")
	
	include(srcdir("sam.jl"))
	
md"Loading modules for the 2023GL104350 paper submission ..."
end

# ╔═╡ e78a75c2-590f-11eb-1144-9127b0309135
md"
# Figure 2. Varying WTG Strengths
"

# ╔═╡ 3c16bd2b-55e0-492f-baf3-aef6a998e9d6
begin
	configDGW = [0.2,0.5,1,2,5,10,20,50,100,200,500]
	configWTG1 = [
		sqrt(2),2,2*sqrt(2.5),5,5*sqrt(2),10,
		10*sqrt(2),20,20*sqrt(2.5),50,50*sqrt(2)
	]
	configWTG2 = [
		1,
		sqrt(2),2,2*sqrt(2.5),5,5*sqrt(2),10,
		10*sqrt(2),20,20*sqrt(2.5),50,50*sqrt(2)
	]
	nconDGW = length(configDGW)
md"Loading time dimension and defining the damping experiments ..."
end

# ╔═╡ 57bbb1c8-e0f8-4b49-8f55-86af376c5167
# ╠═╡ show_logs = false
begin
	pplt.close()
	fig,axs = pplt.subplots(
		ncols=6,aspect=0.3,axwidth=0.8,sharey=0,wspace=[1,6,1,6,1]
	)

	ds_rceprcp = NCDataset(datadir("precipitation","RCE-P1282km300V64.nc"))
	prcp_RCE = ds_rceprcp["precipitation"][:] / 24
	prcpRCEμ_P = mean(prcp_RCE[1001:2000,:])
	close(ds_rceprcp)

	ds_rceprcp = NCDataset(datadir("precipitation","RCE-T1282km300V64.nc"))
	prcp_RCE = ds_rceprcp["precipitation"][:] / 24
	prcpRCEμ_T = mean(prcp_RCE[1001:2000,:])
	close(ds_rceprcp)

	for conDGW in configDGW

		fnc = "DGW-P1282km300V64-$(dampingstrprnt(conDGW)).nc"
		ds_dgwprcp = NCDataset(datadir("precipitation",fnc))

		prcp  = ds_dgwprcp["precipitation"][:] / 24
		for ien = 1 : 15

			prcpii  = prcp[end-2399:end,ien]
			prcpii  = dropdims(mean(reshape(prcpii,6,:),dims=1),dims=1)
			prcpii  = prcpii[.!isnan.(prcpii)]
			prcpμ   = mean(prcpii)
			prcpσ   = zeros(2,1)

			if !isnan(prcpμ)
				prcpσ[1] = prcpμ - minimum(prcpii)
				prcpσ[2] = maximum(prcpii) - prcpμ
				if prcpμ < prcpRCEμ_P * 0.9
					clr = "yellow7"
				elseif prcpμ > prcpRCEμ_P * 1.1
					clr = "blue5"
				else
					clr = "green6"
				end
				axs[1].scatter(prcpμ,conDGW,c=clr,s=50)
				axs[1].errorbar(prcpμ,conDGW,0,prcpσ,c=clr)
			end
			
		end

		close(ds_dgwprcp)

		fnc = "DGW-T1282km300V64-$(dampingstrprnt(conDGW)).nc"
		ds_dgwprcp = NCDataset(datadir("precipitation",fnc))

		prcp  = ds_dgwprcp["precipitation"][:] / 24
		for ien = 1 : 15

			prcpii  = prcp[end-2399:end,ien]
			prcpii  = dropdims(mean(reshape(prcpii,6,:),dims=1),dims=1)
			prcpii  = prcpii[.!isnan.(prcpii)]
			prcpμ   = mean(prcpii)
			prcpσ   = zeros(2,1)

			if !isnan(prcpμ)
				prcpσ[1] = prcpμ - minimum(prcpii)
				prcpσ[2] = maximum(prcpii) - prcpμ
				if prcpμ < prcpRCEμ_T * 0.9
					clr = "yellow7"
				elseif prcpμ > prcpRCEμ_T * 1.1
					clr = "blue5"
				else
					clr = "green6"
				end
				axs[2].scatter(prcpμ,conDGW,c=clr,s=50)
				axs[2].errorbar(prcpμ,conDGW,0,prcpσ,c=clr)
			end
			
		end

		close(ds_dgwprcp)

	end

	for conWTG in configWTG2

		fnc = "TGR-T1282km300V64-$(relaxscalestrprnt(conWTG)).nc"
		ds_wtgprcp = NCDataset(datadir("precipitation",fnc))

		prcp  = ds_wtgprcp["precipitation"][:] / 24
		for ien = 1 : 15

			prcpii  = prcp[end-2399:end,ien]
			prcpii  = dropdims(mean(reshape(prcpii,6,:),dims=1),dims=1)
			prcpii  = prcpii[.!isnan.(prcpii)]
			prcpμ   = mean(prcpii)
			prcpσ   = zeros(2,1)

			if !isnan(prcpμ)
				prcpσ[1] = prcpμ - minimum(prcpii)
				prcpσ[2] = maximum(prcpii) - prcpμ
				if prcpμ < prcpRCEμ_T * 0.9
					clr = "yellow7"
				elseif prcpμ > prcpRCEμ_T * 1.1
					clr = "blue5"
				else
					clr = "green6"
				end
				axs[4].scatter(prcpμ,conWTG,c=clr,s=50)
				axs[4].errorbar(prcpμ,conWTG,0,prcpσ,c=clr)
			end
			
		end

		close(ds_wtgprcp)

		fnc = "SPC-T1282km300V64-$(relaxscalestrprnt(conWTG)).nc"
		ds_wtgprcp = NCDataset(datadir("precipitation",fnc))

		prcp  = ds_wtgprcp["precipitation"][:] / 24
		for ien = 1 : 15

			prcpii  = prcp[end-2399:end,ien]
			prcpii  = dropdims(mean(reshape(prcpii,6,:),dims=1),dims=1)
			prcpii  = prcpii[.!isnan.(prcpii)]
			prcpμ   = mean(prcpii)
			prcpσ   = zeros(2,1)

			if !isnan(prcpμ)
				prcpσ[1] = prcpμ - minimum(prcpii)
				prcpσ[2] = maximum(prcpii) - prcpμ
				if prcpμ < prcpRCEμ_T * 0.9
					clr = "yellow7"
				elseif prcpμ > prcpRCEμ_T * 1.1
					clr = "blue5"
				else
					clr = "green6"
				end
				axs[6].scatter(prcpμ,conWTG,c=clr,s=50)
				axs[6].errorbar(prcpμ,conWTG,0,prcpσ,c=clr)
			end
			
		end

		close(ds_wtgprcp)

	end

	for conWTG in configWTG1

		fnc = "TGR-P1282km300V64-$(relaxscalestrprnt(conWTG)).nc"
		ds_wtgprcp = NCDataset(datadir("precipitation",fnc))

		prcp  = ds_wtgprcp["precipitation"][:] / 24
		for ien = 1 : 15

			prcpii  = prcp[end-2399:end,ien]
			prcpii  = dropdims(mean(reshape(prcpii,6,:),dims=1),dims=1)
			prcpii  = prcpii[.!isnan.(prcpii)]
			prcpμ   = mean(prcpii)
			prcpσ   = zeros(2,1)

			if !isnan(prcpμ)
				prcpσ[1] = prcpμ - minimum(prcpii)
				prcpσ[2] = maximum(prcpii) - prcpμ
				if prcpμ < prcpRCEμ_P * 0.9
					clr = "yellow7"
				elseif prcpμ > prcpRCEμ_P * 1.1
					clr = "blue5"
				else
					clr = "green6"
				end
				axs[3].scatter(prcpμ,conWTG,c=clr,s=50)
				axs[3].errorbar(prcpμ,conWTG,0,prcpσ,c=clr)
			end
			
		end

		close(ds_wtgprcp)

		fnc = "SPC-P1282km300V64-$(relaxscalestrprnt(conWTG)).nc"
		ds_wtgprcp = NCDataset(datadir("precipitation",fnc))

		prcp  = ds_wtgprcp["precipitation"][:] / 24
		for ien = 1 : 15

			prcpii  = prcp[end-2399:end,ien]
			prcpii  = dropdims(mean(reshape(prcpii,6,:),dims=1),dims=1)
			prcpii  = prcpii[.!isnan.(prcpii)]
			prcpμ   = mean(prcpii)
			prcpσ   = zeros(2,1)

			if !isnan(prcpμ)
				prcpσ[1] = prcpμ - minimum(prcpii)
				prcpσ[2] = maximum(prcpii) - prcpμ
				if prcpμ < prcpRCEμ_P * 0.9
					clr = "yellow7"
				elseif prcpμ > prcpRCEμ_P * 1.1
					clr = "blue5"
				else
					clr = "green6"
				end
				axs[5].scatter(prcpμ,conWTG,c=clr,s=50)
				axs[5].errorbar(prcpμ,conWTG,0,prcpσ,c=clr)
			end
			
		end

		close(ds_wtgprcp)

	end
	
	axs[1].plot([1,1]*prcpRCEμ_P,[0,2000],c="grey")
	axs[3].plot([1,1]*prcpRCEμ_P,[0,2000],c="grey")
	axs[5].plot([1,1]*prcpRCEμ_P,[0,2000],c="grey")
	axs[2].plot([1,1]*prcpRCEμ_T,[0,2000],c="grey")
	axs[4].plot([1,1]*prcpRCEμ_T,[0,2000],c="grey")
	axs[6].plot([1,1]*prcpRCEμ_T,[0,2000],c="grey")

	axs[1].format(ltitle="(a) DGW",ylabel=L"$a_m$ / day$^{-1}$")
	axs[3].format(ltitle="(b) TGR",ylabel=L"$\tau$ / hr")
	axs[5].format(ltitle="(c) SPC",ylabel=L"$\tau_1$ / hr")

	for ii in 1 : 2 : 5
		axs[ii].format(
			xscale="symlog",xscale_kw=Dict("linthresh"=>0.01),ultitle="(i) RRTM",
			xlim=(0,10),xlabel=L"Precipitation Rate / mm hr$^{-1}$"
		)
	end
	
	for ii in 1 : 2
		axs[ii].format(ylim=(0.05,2000),yscale="log")
	end
	
	for ii in 3 : 6
		axs[ii].format(ylim=(0.5,200),yscale="log")
	end
	for ii in 2 : 2 : 6
		axs[ii].format(
			xscale="symlog",xscale_kw=Dict("linthresh"=>0.01),
			yticklabels=["","","",""],ytickminor=10:10:100,ultitle="(ii) Ideal",
			xlabel=L"Precipitation Rate / mm hr$^{-1}$",xlim=(0,10)
		)
	end
	
	for ax in axs
		ax.format(lrtitle="Wet",lltitle="Dry")
	end
	
	fig.savefig(plotsdir("fig2-bifurcation.png"),transparent=false,dpi=400)
	load(plotsdir("fig2-bifurcation.png"))
end

# ╔═╡ Cell order:
# ╟─e78a75c2-590f-11eb-1144-9127b0309135
# ╟─681658b0-5914-11eb-0d65-bbace277d145
# ╟─6dce35fc-5914-11eb-0ce2-0d4e164e1898
# ╟─3c16bd2b-55e0-492f-baf3-aef6a998e9d6
# ╟─57bbb1c8-e0f8-4b49-8f55-86af376c5167
