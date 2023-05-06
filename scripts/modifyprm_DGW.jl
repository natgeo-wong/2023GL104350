using DrWatson
@quickactivate "2023GL104350"
using Printf

include(srcdir("sam.jl"))

schname = "DGW"
expname = "P1282km300V64"
tprm  = projectdir("exp","tmp.prm")
plist1 = [1,2,5,10,20,50,100,200,500,1000]

for powerii in plist1
    conii = "damping$(dampingstrprnt(powerii))"
    mkpath(projectdir("exp","prm",schname,expname,conii))
    for imember = 1 : 15
        oprm  = projectdir("exp","prm",schname,expname,"damping001","member01.prm")
        nprm  = projectdir("exp","prm",schname,expname,conii,"member$(@sprintf("%02d",imember)).prm")
        open(tprm,"w") do fprm
            open(oprm,"r") do rprm
                s = read(rprm,String)
                s = replace(s,"member01"=>"member$(@sprintf("%02d",imember))")
                s = replace(s,"nensemble     = 1,"=>"nensemble     = $(imember),")
                s = replace(s,"am_wtg  = 1., "=>"am_wtg  = $(powerii), ")
                write(fprm,s)
            end
        end
        mkpath(projectdir("exp","prm","DGW",expname,conii))
        mv(tprm,nprm,force=true)
    end
end
