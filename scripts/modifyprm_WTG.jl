using DrWatson
@quickactivate "2023GL104350"
using Printf

include(srcdir("sam.jl"))

wtgscheme = "WTG"
expname = "P1282km300V64"
tprm  = projectdir("exp","tmp.prm")
plist = [
    1,sqrt(2),2,2*sqrt(2.5),5,5*sqrt(2),
    10,10*sqrt(2),20,20*sqrt(2.5),50,50*sqrt(2),100
]

for powerii in plist
    conii = relaxscalestrprnt(powerii)
    mkpath(projectdir("exp","prm",wtgscheme,expname,conii))
    for imember = 1 : 15
        oprm  = projectdir("exp","prm",wtgscheme,expname,"relaxscale01d0","member01.prm")
        nprm  = projectdir("exp","prm",wtgscheme,expname,conii,"member$(@sprintf("%02d",imember)).prm")
        open(tprm,"w") do fprm
            open(oprm,"r") do rprm
                s = read(rprm,String)
                s = replace(s,"member01"=>"member$(@sprintf("%02d",imember))")
                s = replace(s,"nensemble     = 1,"=>"nensemble     = $(imember),")
                s = replace(s,
                    "tau_wtg      =   1.000000"=>
                    "tau_wtg      = $(@sprintf("%10f",powerii))"
                )
                write(fprm,s)
            end
        end
        mkpath(projectdir("exp","prm",wtgscheme,expname,conii))
        mv(tprm,nprm,force=true)
    end
end
