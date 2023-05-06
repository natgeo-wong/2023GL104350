using Dates
using DrWatson
using Logging
using NCDatasets
using Printf

function extract2D(
    schname :: String,
    expname :: String,
    runname :: String;
    nt      :: Int = 2000,
    tperday :: Int = 1,
    nmember :: Int = 15,
)

    pa  = zeros(nt,nmember) * NaN
    pw  = zeros(nt,nmember) * NaN
    pr  = zeros(nt,nmember) * NaN
    olr = zeros(nt,nmember) * NaN
    prc = zeros(nt,nmember) * NaN
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
        prc[1:nit,ids] .= ods["PREC"][:]
        pa[1:nit,ids]  .= ods["AREAPREC"][:]
        olr[1:nit,ids] .= ods["LWNT"][:]
        pw[1:nit,ids]  .= ods["PW"][:]
        pr[1:nit,ids]  .= prc[1:nit,ids] ./ pa[1:nit,ids]
        ra = @view pr[1:nit,ids]
        ra[isnan.(ra)] .= 0; ra[ra.==Inf] .= 0
        close(ods)

    end

    fnc = datadir("2D","$(schname)-$(expname)-$(runname).nc")
    if isfile(fnc); rm(fnc,force=true) end

    nds = NCDataset(fnc,"c",attrib = Dict(
        "Conventions" => "CF-1.6",
        "history"     => "Created on $(Dates.now()) using NCDatasets.jl",
        "comments"    => "Creating NetCDF files in the same format that data is saved on the Climate Data Store"
    ))

    nds.dim["time"] = nt
    nds.dim["ensemble"] = nmember

    nctime = defVar(nds,"time",Float64,("time",),attrib = Dict(
        "units"     => "days after model-day 0",
        "full_name" => "Day"
    ))

    ncolr = defVar(nds,"olr",Float64,("time","ensemble"),attrib = Dict(
        "units"     => "W m**-2",
        "long_name" => "outgoing_longwave_radiation",
        "full_name" => "Outgoing Longwave Radiation"
    ))

    ncpw = defVar(nds,"pw",Float64,("time","ensemble"),attrib = Dict(
        "units"     => "kg m**-2",
        "long_name" => "precipitable_water",
        "full_name" => "Precipitable Water"
    ))

    ncra = defVar(nds,"rainarea",Float64,("time","ensemble"),attrib = Dict(
        "units"     => "0-1",
        "long_name" => "fraction_of_area_rain",
        "full_name" => "Fraction of Area with Rainfall"
    ))

    nctime[:] = t
    ncolr[:]  = olr
    ncpw[:]   = pw
    ncra[:]   = pr

    close(nds)

end

function extract2D(
    schname :: String,
    expname :: String;
    nt      :: Int = 2000,
    tperday :: Int = 1,
    nmember :: Int = 15,
)

    pa  = zeros(nt,nmember) * NaN
    pw  = zeros(nt,nmember) * NaN
    pr  = zeros(nt,nmember) * NaN
    olr = zeros(nt,nmember) * NaN
    prc = zeros(nt,nmember) * NaN
    t = 0 : nt
    t = collect(t[1:(end-1)] .+ t[2:end]) / 2
    t = t / tperday

    for ids = 1 : nmember

        ensemble = @sprintf("%02d",ids)
        @info "$(Dates.now()) - Opening $(datadir(
            "archive","$schname","$expname","OUT_STAT",
            "$(schname)_ExploreWTGSpace-$(expname)-member$(ensemble).nc"
        ))"
        ods = NCDataset(datadir(
            "archive","$schname","$expname","OUT_STAT",
            "$(schname)_ExploreWTGSpace-$(expname)-member$(ensemble).nc"
        ))
        nit = ods.dim["time"]
        prc[1:nit,ids] .= ods["PREC"][:]
        pa[1:nit,ids]  .= ods["AREAPREC"][:]
        olr[1:nit,ids] .= ods["LWNT"][:]
        pw[1:nit,ids]  .= ods["PW"][:]
        pr[1:nit,ids]  .= prc[1:nit,ids] ./ pa[1:nit,ids]
        ra = @view pr[1:nit,ids]
        ra[isnan.(ra)] .= 0; ra[ra.==Inf] .= 0
        close(ods)

    end

    fnc = datadir("2D","$(schname)-$(expname).nc")
    if isfile(fnc); rm(fnc,force=true) end

    nds = NCDataset(fnc,"c",attrib = Dict(
        "Conventions" => "CF-1.6",
        "history"     => "Created on $(Dates.now()) using NCDatasets.jl",
        "comments"    => "Creating NetCDF files in the same format that data is saved on the Climate Data Store"
    ))

    nds.dim["time"] = nt
    nds.dim["ensemble"] = nmember

    nctime = defVar(nds,"time",Float64,("time",),attrib = Dict(
        "units"     => "days after model-day 0",
        "full_name" => "Day"
    ))

    ncolr = defVar(nds,"olr",Float64,("time","ensemble"),attrib = Dict(
        "units"     => "W m**-2",
        "long_name" => "outgoing_longwave_radiation",
        "full_name" => "Outgoing Longwave Radiation"
    ))

    ncpw = defVar(nds,"pw",Float64,("time","ensemble"),attrib = Dict(
        "units"     => "kg m**-2",
        "long_name" => "precipitable_water",
        "full_name" => "Precipitable Water"
    ))

    ncra = defVar(nds,"rainarea",Float64,("time","ensemble"),attrib = Dict(
        "units"     => "0-1",
        "long_name" => "fraction_of_area_rain",
        "full_name" => "Fraction of Area with Rainfall"
    ))

    nctime[:] = t
    ncolr[:]  = olr
    ncpw[:]   = pw
    ncra[:]   = pr

    close(nds)

end