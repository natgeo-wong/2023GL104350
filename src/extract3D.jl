using Dates
using DrWatson
using Logging
using NCDatasets
using Printf
using Statistics

function extract3D(
    schname :: String,
    expname :: String,
    runname :: String;
    nt      :: Int = 2000,
    tperday :: Int = 1,
    nmember :: Int = 15,
)

    wwtg = zeros(64,nt,nmember) * NaN
    tabs = zeros(64,nt,nmember) * NaN
    shum = zeros(64,nt,nmember) * NaN
    z = zeros(64,nmember) * NaN
    p = zeros(64,nmember) * NaN
    t = 0 : nt
    t = collect(t[1:(end-1)] .+ t[2:end]) / 2
    t = t / tperday

    for ids = 1 : nmember

        ensemble = @sprintf("%02d",ids)
        @info "$(Dates.now()) - Opening $(datadir(
            "archive","$schname","$expname","$runname","OUT_STAT",
            "$(schname)_ExploreWTGSpace-$(expname)-member$(ensemble).nc"
        ))"
        ods = NCDataset(datadir(
            "archive","$schname","$expname","$runname","OUT_STAT",
            "$(schname)_ExploreWTGSpace-$(expname)-member$(ensemble).nc"
        ))
        nit = ods.dim["time"]
        z[:,ids] .= ods["z"][:]
        p[:,ids] .= ods["p"][:]
        wwtg[:,1:nit,ids] .= ods["WWTG"][:]
        tabs[:,1:nit,ids] .= ods["TABS"][:]
        shum[:,1:nit,ids] .= ods["QV"][:] ./ 1000
        close(ods)

    end

    fnc = datadir("3D","$(schname)-$(expname)-$(runname).nc")
    if isfile(fnc); rm(fnc,force=true) end

    nds = NCDataset(fnc,"c",attrib = Dict(
        "Conventions" => "CF-1.6",
        "history"     => "Created on $(Dates.now()) using NCDatasets.jl",
        "comments"    => "Creating NetCDF files in the same format that data is saved on the Climate Data Store"
    ))

    nds.dim["time"] = nt
    nds.dim["level"] = 64
    nds.dim["ensemble"] = nmember

    nctime = defVar(nds,"time",Float64,("time",),attrib = Dict(
        "units"     => "days after model-day 0",
        "full_name" => "Day"
    ))

    ncz = defVar(nds,"z",Float64,("level",),attrib = Dict(
        "units"     => "m",
        "full_name" => "height"
    ))

    ncp = defVar(nds,"p",Float64,("level",),attrib = Dict(
        "units"     => "hPa",
        "full_name" => "pressure_level"
    ))

    ncwwtg = defVar(nds,"wwtg",Float64,("level","time","ensemble"),attrib = Dict(
        "units"     => "m s**-1",
        "long_name" => "weak_temperature_gradient_vertical_velocity",
        "full_name" => "WTG Vertical Velocity"
    ))

    nctabs = defVar(nds,"t",Float64,("level","time","ensemble"),attrib = Dict(
        "units"     => "K",
        "long_name" => "temperature",
        "full_name" => "Atmospheric Temperature"
    ))

    ncshum = defVar(nds,"q",Float64,("level","time","ensemble"),attrib = Dict(
        "units"     => "kg kg**-1",
        "long_name" => "specifc_humidity",
        "full_name" => "Specific Humidity"
    ))

    nctime[:] = t
    ncz[:] = dropdims(mean(z,dims=2),dims=2)
    ncp[:] = dropdims(mean(p,dims=2),dims=2)
    ncwwtg[:]  = wwtg
    nctabs[:]  = tabs
    ncshum[:]  = shum

    close(nds)

end