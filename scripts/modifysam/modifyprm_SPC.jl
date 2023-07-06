using DrWatson
@quickactivate "2023GL104350"
using Logging
using Printf

include(srcdir("sam.jl"))

expname = "P1282km300V64"
tprm  = projectdir("exp","tmp.prm")
plist1 = [
    0,sqrt(2),2,2*sqrt(2.5),5,5*sqrt(2),
    10,10*sqrt(2),20,20*sqrt(2.5),50,50*sqrt(2),
]
plist2 = vcat(
    plist1/100,1,
    plist1
)

for powerii in plist1
    conii = relaxscalestrprnt(powerii)
    mkpath(projectdir("exp","prm","SPC",expname,conii))
    for imember = 1 : 15
        mstr = @sprintf("%02d",imember)
        oprm = projectdir("run","modifysam","prmtemplates","SPC_$(expname).prm")
        nprm = projectdir("exp","prm","SPC",expname,conii,"member$(mstr).prm")
        open(tprm,"w") do fprm
            open(oprm,"r") do rprm
                s = read(rprm,String)
                s = replace(s,"[xx]" => mstr)
                s = replace(s,"[en]" => "$(imember)")
                s = replace(s,"[am]" => @sprintf("%7e",1))
                if !iszero(powerii)
                    s = replace(s,"[bool]" => "true")
                    s = replace(s,"[tau]" => @sprintf("%7e",powerii))
                else
                    s = replace(s,"[bool]" => "false")
                    s = replace(s,"[tau]" => @sprintf("%7e",1))
                end
                write(fprm,s)
            end
        end
        mkpath(projectdir("exp","prm","SPC",expname,conii))
        mv(tprm,nprm,force=true)
        @info "Creating new prm file for SPC $expname $conii ensemble member $imember"
    end
end
