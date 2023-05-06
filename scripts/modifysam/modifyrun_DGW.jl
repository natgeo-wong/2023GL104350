using DrWatson
@quickactivate "ExploreWTGSpace"
using Printf

include(srcdir("sam.jl"))

schname = "DGW"
expname = "P1282km300V64"
pwrvec1 = [1,2,5,10,20,50,100,200,500,1000]
pwrvec2 = [1,2,5,10,14,20,32,50,71,100,200,500,1000]
pwrvec  = pwrvec1

mrun = projectdir("run","modelrun.sh")
brun = projectdir("run","Build.csh")

open(mrun,"r") do frun
    s = read(frun,String)
    for pwrii in pwrvec

        pwrname = "damping$(dampingstrprnt(pwrii))"

        for ensembleii in 1 : 15

            ensname = "ensemble$(@sprintf("%02d",ensembleii))"
            nrun = projectdir("run",schname,expname,pwrname,"$(ensname).sh")

            open(nrun,"w") do wrun
                sn = replace(s ,"[email]"=>"")
                sn = replace(sn,"[user]"=>"")
                sn = replace(sn,"[project]"=>"ExploreWTGSpace")
                sn = replace(sn,"[experiment]"=>"$(expname)")
                sn = replace(sn,"[config]"=>"$(pwrname)")
                if ensembleii < 6
                    sn = replace(sn,"[sndname]"=>"$(expname)")
                elseif (ensembleii > 5) && (ensembleii < 11)
                    sn = replace(sn,"[sndname]"=>"$(expname)_hot")
                else
                    sn = replace(sn,"[sndname]"=>"$(expname)_cld")
                end
                sn = replace(sn,"[lsfname]"=>"noforcing")
                sn = replace(sn,"[schname]"=>"$(schname)")
                sn = replace(sn,"member[xx]"=>"member$(@sprintf("%02d",ensembleii))")
                write(wrun,sn)
            end

        end

    end
end

open(brun,"r") do frun
    s = read(frun,String)
    for pwrii in pwrvec

        pwrname = "damping$(dampingstrprnt(pwrii))"
        nrun = projectdir("run",schname,expname,pwrname,"Build.csh")

        open(nrun,"w") do wrun
            sn = replace(s ,"[user]"=>"")
            sn = replace(sn,"[schname]"=>schname)
            sn = replace(sn,"[expname]"=>"$(expname)")
            sn = replace(sn,"[runname]"=>"$(pwrname)")
            write(wrun,sn)
        end

    end
end