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
	using DSP
	using NCDatasets
	using PlutoUI
	using Printf
	using StatsBase
	
	using ImageShow, PNGFiles
	using PyCall, LaTeXStrings
	pplt = pyimport("proplot")
	
	include(srcdir("sam.jl"))
	
md"Loading modules for the 2023GL104350 project..."
end

# ╔═╡ e78a75c2-590f-11eb-1144-9127b0309135
md"
# Figure 3. Time Series, Power Spectrum
"

# ╔═╡ a63de98c-5b35-11eb-0a8f-b7a1ebd441b6
begin
	configDGW = [0.01,0.02,0.05,0.1,0.2,0.5,1,2,5,10,20,50,100,200,500,1000]
	configWTG = [
		0.01,0.02,0.02*sqrt(2.5),0.05,0.05*sqrt(2),
		0.1,0.2,0.2*sqrt(2.5),0.5,0.5*sqrt(2),
		1,sqrt(2),2,2*sqrt(2.5),5,5*sqrt(2),
		10,10*sqrt(2),20,20*sqrt(2.5),50,50*sqrt(2),100
	]
	nconDGW = length(configDGW)
	nconWTG = length(configWTG)
    blues_DGW = pplt.get_colors("Blues",(nconDGW+2))
	blues_WTG = pplt.get_colors("Blues",(nconWTG+2))
	lgd_DGW = Dict("frame"=>false,"ncols"=>2)
	lgd_WTG = Dict("frame"=>false,"ncols"=>3)
	md"Loading time dimension and defining the damping experiments ..."
end

# ╔═╡ 864f1d31-c629-412b-825a-98fe9398b591
begin
	pplt.close()
	fts,ats = pplt.subplots(nrows=2,aspect=2,axwidth=3.5)

	for ic in 1 : nconDGW

		fnc = "DGW-S1284km300V64-damping$(dampingstrprnt(configDGW[ic])).nc"
		ds_dgwprcp = NCDataset(datadir("precipitation",fnc))

		tdgw    = ds_dgwprcp["time"][:]; tdgw = reshape(tdgw,8,:)
		tdgw    = dropdims(mean(tdgw,dims=1),dims=1)
		prcpdgw = ds_dgwprcp["precipitation"][:] / 24
		
		
		for ien = 1 : 15

			prcpii = prcpdgw[:,ien]
			prcpii = reshape(prcpii,8,:)
			prcpii = dropdims(mean(prcpii,dims=1),dims=1)

			if ien == 1
				constr = @sprintf("%.1e",configDGW[ic])
				ats[1].plot(
					tdgw,prcpii,color=blues_DGW[ic+1],
					label=(L"$a_m =$" * " $(constr)" * L" day$^{-1}$"),
					legend="r",legend_kw=lgd_DGW
				)
			else
				ats[1].plot(tdgw,prcpii,color=blues_DGW[ic+1])
			end
			
		end

		close(ds_dgwprcp)

	end

	for ic in 1 : nconWTG

		fnc = "WTG-S1284km300V64-$(relaxscalestrprnt(configWTG[ic])).nc"
		ds_wtgprcp = NCDataset(datadir("precipitation",fnc))

		twtg    = ds_wtgprcp["time"][:]; twtg = reshape(twtg,8,:)
		twtg    = dropdims(mean(twtg,dims=1),dims=1)
		prcpwtg = ds_wtgprcp["precipitation"][:] / 24
		
		
		for ien = 1 : 15

			prcpii = prcpwtg[:,ien]
			prcpii = reshape(prcpii,8,:)
			prcpii = dropdims(mean(prcpii,dims=1),dims=1)

			if ien == 1
				constr = @sprintf("%.1e",configWTG[ic])
				ats[2].plot(
					twtg,prcpii,color=blues_WTG[ic+1],
					label=(L"$\tau =$" * " $(constr) hr"),
					legend="r",legend_kw=lgd_WTG
				)
			else
				ats[2].plot(twtg,prcpii,color=blues_WTG[ic+1])
			end
			
		end

		close(ds_wtgprcp)

	end

	ds_rceprcp = NCDataset(datadir("precipitation","RCE-S1284km300V64.nc"))

	t_RCE = ds_rceprcp["time"][:]
	prcp_RCE = ds_rceprcp["precipitation"][:] / 24
	ats[1].plot(t_RCE,prcp_RCE[:,1],c="k",label="RCE",legend="r")
	ats[2].plot(t_RCE,prcp_RCE[:,1],c="k",label="RCE",legend="r")
	ats[1].plot(t_RCE,prcp_RCE[:,2:10],c="k")
	ats[2].plot(t_RCE,prcp_RCE[:,2:10],c="k")

	ats[1].format(ultitle="(a) Precipitation Time-Series (DGW)")
	ats[2].format(ultitle="(b) Precipitation Time-Series (TGR)")

	close(ds_rceprcp)
	
	for ax in ats
		ax.format(
			xlim=(000,250),#yscale="symlog",yscale_kw=Dict("linthresh"=>0.01),
			ylim=(0,0.5),
			ylabel=L"Rainfall Rate / mm hr$^{-1}$",xlabel="Days"
		)
	end
	
	fts.savefig(plotsdir("fig3-timeseries.png"),transparent=false,dpi=400)
	load(plotsdir("fig3-timeseries.png"))
end

# ╔═╡ 2a3c3054-90ad-4f2b-a4d2-bdfa45983377
begin
	signalpower_DGW = zeros(401,nconDGW,15)
	signalfreq_DGW  = zeros(401,nconDGW,15)
	totalmember_DGW = zeros(1,nconDGW)
	for icon = 1 : nconDGW

		fnc = "DGW-S1284km300V64-damping$(dampingstrprnt(configDGW[icon])).nc"
		ds_dgwprcp = NCDataset(datadir("precipitation",fnc))
		prcpdgw = ds_dgwprcp["precipitation"][:]
		close(ds_dgwprcp)
		
		nmem = 0
		for imem = 1 : 15
			prcp = prcpdgw[:,imem]
			if sum(.!isnan.(prcp)) == 2000
				nmem += 1
				prcp = prcp[(end-799):end] .- mean(prcp[(end-799):end])
				pdg = periodogram(prcp,fs=8)
				signalpower_DGW[:,icon,imem] .= pdg.power
				signalfreq_DGW[:,icon,imem]  .= pdg.freq
			end
		end
		totalmember_DGW[icon] = nmem
	end
	signalpower_DGW = dropdims(sum(signalpower_DGW,dims=3),dims=3)
	signalpower_DGW = signalpower_DGW ./ totalmember_DGW
	signalfreq_DGW  = dropdims(sum(signalfreq_DGW ,dims=3),dims=3)
	signalfreq_DGW  = signalfreq_DGW  ./ totalmember_DGW
	signalfreq_DGW  = dropdims(mean(signalfreq_DGW[:,2:end],dims=2),dims=2)
end

# ╔═╡ 6fc11aaf-0b10-468b-8943-4675715dc11f
begin
	signalpower_TGR = zeros(401,nconWTG,15)
	signalfreq_TGR  = zeros(401,nconWTG,15)
	totalmember_TGR = zeros(1,nconWTG)
	for icon = 1 : nconWTG

		fnc = "WTG-S1284km300V64-$(relaxscalestrprnt(configWTG[icon])).nc"
		ds_wtgprcp = NCDataset(datadir("precipitation",fnc))
		prcpwtg = ds_wtgprcp["precipitation"][:]
		close(ds_wtgprcp)
		
		nmem = 0
		for imem = 1 : 15
			prcp = prcpwtg[:,imem]
			if sum(.!isnan.(prcp)) == 2000
				nmem += 1
				prcp = prcp[(end-799):end] .- mean(prcp[(end-799):end])
				pdg = periodogram(prcp,fs=8)
				signalpower_TGR[:,icon,imem] .= pdg.power
				signalfreq_TGR[:,icon,imem]  .= pdg.freq
			end
		end
		totalmember_TGR[icon] = nmem
	end
	signalpower_TGR = dropdims(sum(signalpower_TGR,dims=3),dims=3)
	signalpower_TGR = signalpower_TGR ./ totalmember_TGR
	signalfreq_TGR  = dropdims(sum(signalfreq_TGR ,dims=3),dims=3)
	signalfreq_TGR  = signalfreq_TGR  ./ totalmember_TGR
	signalfreq_TGR  = dropdims(mean(signalfreq_TGR,dims=2),dims=2)
end

# ╔═╡ 6a21129c-9a25-4e32-8fc8-34106ee9279e
begin
	pplt.close(); fps,aps = pplt.subplots(ncols=2,aspect=2,axwidth=3,sharey=0)

	lvls = 10:10:100

	c = aps[1].contourf(
		1 ./signalfreq_DGW[2:end],configDGW,
		(signalpower_DGW)'[:,2:end],
		levels=lvls,extend="both",cmap="fire"
	)
	aps[1].format(
		ylabel=L"$a_m$ / day$^{-1}$",ylim=(0.01,100),
		ultitle="(c) Power Spectral Density (DGW)"
	)

	c = aps[2].contourf(
		1 ./signalfreq_TGR[2:end],configWTG,
		(signalpower_TGR)'[:,2:end],
		levels=lvls,extend="both",cmap="fire"
	)
	aps[2].format(
		ylabel=L"$\tau$ / hr",ylim=(0.01,100),
		ultitle="(d) Power Spectral Density (TGR)",ytickloc="right"
	)

	for ax in aps
		ax.format(
			xscale="log",xlim=(0.2,50),xlabel="Period / Days",
			yscale="log",yformatter="log"
		)
	end

	fps.colorbar(c,length=0.75,label="Power / dB")
	fps.savefig(plotsdir("fig3-powerspectrum.png"),transparent=false,dpi=400)
	load(plotsdir("fig3-powerspectrum.png"))
end

# ╔═╡ Cell order:
# ╟─e78a75c2-590f-11eb-1144-9127b0309135
# ╟─681658b0-5914-11eb-0d65-bbace277d145
# ╟─6dce35fc-5914-11eb-0ce2-0d4e164e1898
# ╟─a63de98c-5b35-11eb-0a8f-b7a1ebd441b6
# ╟─864f1d31-c629-412b-825a-98fe9398b591
# ╟─2a3c3054-90ad-4f2b-a4d2-bdfa45983377
# ╟─6fc11aaf-0b10-468b-8943-4675715dc11f
# ╟─6a21129c-9a25-4e32-8fc8-34106ee9279e
