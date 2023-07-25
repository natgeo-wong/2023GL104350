using Dates
using DrWatson
using Logging
using NCDatasets
using Printf

function extractprecip(
    schname :: String,
    expname :: String,
    runname :: String;
    nt      :: Int = 2000,
    tperday :: Int = 1,
    nmember :: Int = 15,
)

    prcp = zeros(nt,nmember) * NaN
    t = 0 : nt
    t = collect(t[1:(end-1)] .+ t[2:end]) / 2
    t = t / tperday

    for ids = 1 : nmember

        ensemble = @sprintf("%02d",ids)
        @info "$(Dates.now()) - Opening $(datadir(
            "$schname","$expname","$runname","OUT_STAT",
            "$(schname)_ExploreWTGSpace-$(expname)-member$(ensemble).nc"
        ))"
        fnc = datadir(
            "$schname","$expname","$runname","OUT_STAT",
            "$(schname)_ExploreWTGSpace-$(expname)-member$(ensemble).nc"
        )
        if isfile(fnc)
            ods = NCDataset(fnc)
            try
                nit = ods.dim["time"]
                prcp[1:nit,ids] .= ods["PREC"][:]
            catch
                @warn "Unable to extract precipitation data from $(fnc)"
            end
            close(ods)
        else
            @warn "No file exists at $(fnc), please run this configuration again ..."
        end

    end

    fnc = datadir("precipitation","$(schname)-$(expname)-$(runname).nc")
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

    ncvar = defVar(nds,"precipitation",Float64,("time","ensemble"),attrib = Dict(
        "units"     => "mm day**-1",
        "long_name" => "daily_precipitation_rate",
        "full_name" => "Daily Precipitation Rate"
    ))

    nctime[:] = t
    ncvar[:]  = prcp

    close(nds)

end

function extractprecip(
    schname :: String,
    expname :: String;
    nt      :: Int = 2000,
    tperday :: Int = 1,
    nmember :: Int = 15,
)

    prcp = zeros(nt,nmember) * NaN
    t = 0 : nt
    t = collect(t[1:(end-1)] .+ t[2:end]) / 2
    t = t / tperday

    for ids = 1 : nmember

        ensemble = @sprintf("%02d",ids)
        @info "$(Dates.now()) - Opening $(datadir(
            "$schname","$expname","OUT_STAT",
            "$(schname)_ExploreWTGSpace-$(expname)-member$(ensemble).nc"
        ))"
        ods = NCDataset(datadir(
            "$schname","$expname","OUT_STAT",
            "$(schname)_ExploreWTGSpace-$(expname)-member$(ensemble).nc"
        ))
        nit = ods.dim["time"]
        prcp[1:nit,ids] .= ods["PREC"][:]
        close(ods)

    end

    fnc = datadir("precipitation","$(schname)-$(expname).nc")
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

    ncvar = defVar(nds,"precipitation",Float64,("time","ensemble"),attrib = Dict(
        "units"     => "mm day**-1",
        "long_name" => "daily_precipitation_rate",
        "full_name" => "Daily Precipitation Rate"
    ))

    nctime[:] = t
    ncvar[:]  = prcp

    close(nds)

end