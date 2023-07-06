using DrWatson
@quickactivate "2023GL104350"
using Logging
using Printf

include(srcdir("sam.jl"))

schname = "DGW"
expname = "P1282km300V64"
tprm  = projectdir("exp","tmp.prm")
plist1 = [0,0.2,0.5,1,2,5,10,20,50,100,200,500]
plist2 = [0,0.02,0.05,0.1,0.2,0.5,1,2,5,10,20,50,100,200,500]

for powerii in plist2
    conii = dampingstrprnt(powerii)
    mkpath(projectdir("exp","prm",schname,expname,conii))
    for imember = 1 : 15
        oprm  = projectdir("scripts","modifysam","prm","$(schname)_$(expname).prm")
        nprm  = projectdir("exp","prm",schname,expname,conii,"member$(@sprintf("%02d",imember)).prm")
        open(tprm,"w") do fprm
            open(oprm,"r") do rprm
                s = read(rprm,String)
                s = replace(s,"[xx]" => @sprintf("%02d",imember))
                s = replace(s,"[en]" => "$(imember)")
                s = replace(s,"[tau]" => @sprintf("%7e",1))
                if !iszero(powerii)
                    s = replace(s,"[bool]" => "true")
                    s = replace(s,"[am]" => @sprintf("%7e",powerii))
                else
                    s = replace(s,"[bool]" => "false")
                    s = replace(s,"[am]" => @sprintf("%7e",1))
                end
                write(fprm,s)
            end
        end
        mkpath(projectdir("exp","prm","DGW",expname,conii))
        mv(tprm,nprm,force=true)
        @info "Creating new prm file for $schname $expname $conii ensemble member $imember"
    end
end
