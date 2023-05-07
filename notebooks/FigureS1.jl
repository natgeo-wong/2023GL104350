### A Pluto.jl notebook ###
# v0.19.25

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
	using NumericalIntegration
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
# Figure S1. Realistic Radiation
"

# ╔═╡ 3c16bd2b-55e0-492f-baf3-aef6a998e9d6
begin
	configDGW = [1,2,5,10,20,50,100,200,500,1000]
	configWTG = [
		1,sqrt(2),2,2*sqrt(2.5),5,5*sqrt(2),10,
		10*sqrt(2),20,20*sqrt(2.5),50,50*sqrt(2),100
	]
	nconDGW = length(configDGW)
	nconWTG = length(configWTG)
	bluesDGW = pplt.get_colors("Blues",(nconDGW+2))
	bluesWTG = pplt.get_colors("Blues",(nconWTG+2))
	grnsDGW  = pplt.get_colors("Teal",(nconDGW+2))
	grnsWTG  = pplt.get_colors("Teal",(nconWTG+2))
	brwnsDGW = pplt.get_colors("Brown",(nconDGW+2))
	brwnsWTG = pplt.get_colors("Brown",(nconWTG+2))
	lgd = Dict("frame"=>false,"ncols"=>4)
md"Loading time dimension and defining the damping experiments ..."
end

# ╔═╡ 8b857b61-6607-41fb-96ae-64929f64a22f
function calcT2e(T::Real)

    tb = T - 273.15
    if tb <= 0
    	return exp(43.494 - 6545.8/(tb+278)) / (tb+868)^2
    else
    	return exp(34.494 - 4924.99/(tb+237.1)) / (tb+105)^1.57
    end

end

# ╔═╡ 44153543-3861-4435-bb53-bf57a9ab5be3
function calcT2q(T::Real,p::Real)
    e = calcT2e(T)
    return e * 0.621981 / (p - 0.378019 * e)
end

# ╔═╡ 57bbb1c8-e0f8-4b49-8f55-86af376c5167
# ╠═╡ show_logs = false
begin
	pplt.close()
	fig,axs = pplt.subplots(
			ncols=5,nrows=2,aspect=0.4,axwidth=1,
			sharex=0,sharey=0,wspace=1.5
	)


	ds_rceprcp = NCDataset(datadir("precipitation","RCE-P1282km300V64.nc"))
	ds_rce2D   = NCDataset(datadir("2D","RCE-P1282km300V64.nc"))
	ds_rce3D   = NCDataset(datadir("3D","RCE-P1282km300V64.nc"))

	prcp_RCE = ds_rceprcp["precipitation"][:] / 24
	axs[1].plot([1,1]*mean(prcp_RCE[501:1000,:]),[0,2000],c="grey")
	axs[6].plot([1,1]*mean(prcp_RCE[501:1000,:]),[0,2000],c="grey")
	
	parea = ds_rce2D["rainarea"][:]
	axs[2].plot([1,1]*mean(parea[501:1000,:]),[0,2000],c="grey")
	axs[7].plot([1,1]*mean(parea[501:1000,:]),[0,2000],c="grey")
	
	prate = ds_rce2D["rainarearate"][:] / 24
	axs[3].plot([1,1]*mean(prate[501:1000,:]),[0,2000],c="grey")
	axs[8].plot([1,1]*mean(prate[501:1000,:]),[0,2000],c="grey")
	
	pw = ds_rce2D["pw"][:]
	axs[4].plot([1,1]*mean(pw[501:1000,:]),[0,2000],c="grey")
	axs[9].plot([1,1]*mean(pw[501:1000,:]),[0,2000],c="grey")

	p  = ds_rce3D["p"][:] * 100
	ta = ds_rce3D["t"][:,1:1000,:]
	ta = dropdims(mean(ta,dims=2),dims=2)
	qs = zeros(size(ta))
	sw = zeros(10)
	for ien = 1 : 10, ip = 1 : 64
		qs[ip,ien] = calcT2q(ta[ip,ien],p[ip])
	end
	for ien = 1 : 10
		sw[ien] = integrate(reverse(p),reverse(qs[:,ien])) / 9.81
	end
	crh = mean(pw[501:1000,:]) ./ sw * 100
	axs[5].plot([1,1]*mean(crh),[0,2000],c="grey")
	axs[10].plot([1,1]*mean(crh),[0,2000],c="grey")

	close(ds_rceprcp)
	close(ds_rce2D)
	close(ds_rce3D)

	for conDGW in configDGW

		fnc = "DGW-P1282km300V64-damping$(dampingstrprnt(conDGW)).nc"
		ds_dgwprcp = NCDataset(datadir("precipitation",fnc))
		ds_dgw2D   = NCDataset(datadir("2D",fnc))
		ds_dgw3D   = NCDataset(datadir("3D",fnc))

		prcp  = ds_dgwprcp["precipitation"][:] / 24
		parea = ds_dgw2D["rainarea"][:]
		prate = ds_dgw2D["rainarearate"][:] / 24
		pw    = ds_dgw2D["pw"][:]
		for ien = 1 : 15

			prcpii = prcp[end-1000:end,ien]
			prcpii = mean(prcpii[.!isnan.(prcpii)])
			pareaii = parea[end-1000:end,ien]
			pareaii = mean(pareaii[.!isnan.(pareaii)])
			prateii = prate[end-1000:end,ien]
			prateii = mean(prateii[.!isnan.(prateii)])
			pwii = pw[end-1000:end,ien]
			pwii = mean(pwii[.!isnan.(pwii)])

			p  = ds_dgw3D["p"][:] * 100
			ta = ds_dgw3D["t"][:,end-1000:end,ien]
			ta = dropdims(mean(ta,dims=2),dims=2)
			qs = zeros(size(ta))
			for ip = 1 : 64
				qs[ip] = calcT2q(ta[ip],p[ip])
			end
			sw = integrate(reverse(p),reverse(qs)) / 9.81
			crhii = mean(pwii) ./ sw * 100
			
			if (ien >= 6) && (ien <= 10)
				clr = "yellow7"
			else
				clr = "blue5"
			end

			if ien <= 5
				axs[1].plot(prcpii,conDGW,marker=".",c="k",ms=4)
				axs[2].plot(pareaii,conDGW,marker=".",c="k",ms=4)
				axs[3].plot(prateii,conDGW,marker=".",c="k",ms=4)
				axs[4].plot(pwii,conDGW,marker=".",c="k",ms=4)
				axs[5].plot(crhii,conDGW,marker=".",c="k",ms=4)
			else
				axs[1].scatter(prcpii,conDGW,c=clr,alpha=0.1,s=50)
				axs[2].scatter(pareaii,conDGW,c=clr,alpha=0.1,s=50)
				axs[3].scatter(prateii,conDGW,c=clr,alpha=0.1,s=50)
				axs[4].scatter(pwii,conDGW,c=clr,alpha=0.1,s=50)
				axs[5].scatter(crhii,conDGW,c=clr,alpha=0.1,s=50)
			end
			
		end

		close(ds_dgwprcp)
		close(ds_dgw2D)
		close(ds_dgw3D)

	end

	for conWTG in configWTG

		fnc = "WTG-P1282km300V64-$(relaxscalestrprnt(conWTG)).nc"
		ds_wtgprcp = NCDataset(datadir("precipitation",fnc))
		ds_wtg2D   = NCDataset(datadir("2D",fnc))
		ds_wtg3D   = NCDataset(datadir("3D",fnc))

		prcp  = ds_wtgprcp["precipitation"][:] / 24
		parea = ds_wtg2D["rainarea"][:]
		prate = ds_wtg2D["rainarearate"][:] / 24
		pw    = ds_wtg2D["pw"][:]
		for ien = 1 : 15

			prcpii = prcp[end-1000:end,ien]
			prcpii = mean(prcpii[.!isnan.(prcpii)])
			pareaii = parea[end-1000:end,ien]
			pareaii = mean(pareaii[.!isnan.(pareaii)])
			prateii = prate[end-1000:end,ien]
			prateii = mean(prateii[.!isnan.(prateii)])
			pwii = pw[end-1000:end,ien]
			pwii = mean(pwii[.!isnan.(pwii)])

			p  = ds_wtg3D["p"][:] * 100
			ta = ds_wtg3D["t"][:,end-1000:end,ien]
			ta = dropdims(mean(ta,dims=2),dims=2)
			qs = zeros(size(ta))
			for ip = 1 : 64
				qs[ip] = calcT2q(ta[ip],p[ip])
			end
			sw = integrate(reverse(p),reverse(qs)) / 9.81
			crhii = mean(pwii) ./ sw * 100
			
			if (ien >= 6) && (ien <= 10)
				clr = "yellow7"
			else
				clr = "blue5"
			end

			if ien <= 5
				axs[6].plot(prcpii,conWTG,marker=".",c="k",ms=4)
				axs[7].plot(pareaii,conWTG,marker=".",c="k",ms=4)
				axs[8].plot(prateii,conWTG,marker=".",c="k",ms=4)
				axs[9].plot(pwii,conWTG,marker=".",c="k",ms=4)
				axs[10].plot(crhii,conWTG,marker=".",c="k",ms=4)
			else
				axs[6].scatter(prcpii,conWTG,c=clr,alpha=0.1,s=50)
				axs[7].scatter(pareaii,conWTG,c=clr,alpha=0.1,s=50)
				axs[8].scatter(prateii,conWTG,c=clr,alpha=0.1,s=50)
				axs[9].scatter(pwii,conWTG,c=clr,alpha=0.1,s=50)
				axs[10].scatter(crhii,conWTG,c=clr,alpha=0.1,s=50)
			end
			
		end

		close(ds_wtgprcp)
		close(ds_wtg2D)
		close(ds_wtg3D)

	end
	
	axs[1].format(
		ylabel=L"$a_m$ / day$^{-1}$",
		xscale="symlog",xscale_kw=Dict("linthresh"=>0.01),
		xlim=(0,10),ultitle="(a)",
		ltitle=L"(1) Damped Gravity Wave Implementation | Sensitivity to $a_m$",
	)
	
	axs[2].format(
		xscale="symlog",xscale_kw=Dict("linthresh"=>0.01),
		xlim=(0,1),ultitle="(b)",
	)
	
	axs[3].format(
		xscale="symlog",xscale_kw=Dict("linthresh"=>1),
		xlim=(0,10),ultitle="(c)",
	)
	
	axs[4].format(xlim=(0,75),ultitle="(d)")
	axs[5].format(xlim=(0,100),ultitle="(e)")
	
	for ii in 1 : 5
		axs[ii].format(ylim=(0.5,2000),yscale="log")
	end
	for ii in 2 : 5
		axs[ii].format(yticklabels=["","","",""],ytickminor=10:10:100)
	end
	
	axs[6].format(
		ylabel=L"$\tau$ / hr",
		xscale="symlog",xscale_kw=Dict("linthresh"=>0.01),
		xlim=(0,10),xlabel=L"Domain $P$ / mm hr$^{-1}$",ultitle="(a)",
		ltitle=L"(2) Temperature Gradient Relaxation Implementation | Sensitivity to $\tau$"
	)
	
	axs[7].format(
		xscale="symlog",xscale_kw=Dict("linthresh"=>0.01),
		xlim=(0,1),xlabel="Rain Area Fraction",ultitle="(b)",
	)
	
	axs[8].format(
		xscale="symlog",xscale_kw=Dict("linthresh"=>1),
		xlim=(0,10),xlabel=L"Rain Area $P$ / mm hr$^{-1}$",ultitle="(c)",
	)
	
	axs[9].format(xlim=(0,75),xlabel="PWV / mm",ultitle="(d)")
	axs[10].format(xlim=(0,100),xlabel="CRH / %",ultitle="(e)")
	
	for ii in 6 : 10
		axs[ii].format(ylim=(0.5,200),yscale="log")
	end
	for ii in 7 : 10
		axs[ii].format(yticklabels=["","","",""],ytickminor=10:10:100)
	end
	
	for ax in axs
		ax.format(lrtitle="Wet",lltitle="Dry")
	end
	
	fig.savefig(plotsdir("figS1-bifurcation.png"),transparent=false,dpi=400)
	load(plotsdir("figS1-bifurcation.png"))
end

# ╔═╡ Cell order:
# ╟─e78a75c2-590f-11eb-1144-9127b0309135
# ╟─681658b0-5914-11eb-0d65-bbace277d145
# ╠═6dce35fc-5914-11eb-0ce2-0d4e164e1898
# ╟─3c16bd2b-55e0-492f-baf3-aef6a998e9d6
# ╠═8b857b61-6607-41fb-96ae-64929f64a22f
# ╠═44153543-3861-4435-bb53-bf57a9ab5be3
# ╟─57bbb1c8-e0f8-4b49-8f55-86af376c5167
