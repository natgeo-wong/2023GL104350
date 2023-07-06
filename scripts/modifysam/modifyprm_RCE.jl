using DrWatson
@quickactivate "2023GL104350"
using Printf

expname = "P1282km300V64"
tprm  = projectdir("exp","tmp.prm")

for imember = 1 : 10
    mstr = @sprintf("%02d",imember)
    oprm  = projectdir("exp","prmtemplates","RCE","DGW_$(expname).prm")
    nprm  = projectdir("exp","prm","RCE",expname,"member$(mstr).prm")
    open(tprm,"w") do fprm
        open(oprm,"r") do rprm
            s = read(rprm,String)
            s = replace(s,"[xx]" => mstr)
            s = replace(s,"[en]" => "$(imember)")
            write(fprm,s)
        end
    end
    mv(tprm,nprm,force=true)
end
