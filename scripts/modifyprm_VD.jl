using DrWatson
@quickactivate "2023GL104350"
using Printf

# include(srcdir("sam.jl"))

schname = "VDT"
tprm  = projectdir("exp","tmp.prm")
expvec = [
    "H0d00F1d00","H0d50F1d00","H1d00F1d00",
    "H1d00F0d50","H1d00F0d40","H1d00F0d30",
    "H1d00F0d20","H1d00F0d10","H1d00F0d00"
]

for expname in expvec
    mkpath(projectdir("exp","prm",schname,expname))
    for imember = 1 : 15
        oprm  = projectdir("exp","prm",schname,expname,"member01.prm")
        nprm  = projectdir("exp","prm",schname,expname,"member$(@sprintf("%02d",imember)).prm")
        open(tprm,"w") do fprm
            open(oprm,"r") do rprm
                s = read(rprm,String)
                s = replace(s,"S1284km300V64"=>expname)
                s = replace(s,"member01"=>"member$(@sprintf("%02d",imember))")
                s = replace(s,"nensemble     = 1,"=>"nensemble     = $(imember),")
                write(fprm,s)
            end
        end
        mv(tprm,nprm,force=true)
    end
end
