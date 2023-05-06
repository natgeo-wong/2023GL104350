using DrWatson
@quickactivate "ExploreWTGSpace"
using Printf

include(srcdir("sam.jl"))

mrun = projectdir("run","modelrun.sh")
brun = projectdir("run","Build.csh")

open(mrun,"r") do frun
    s = read(frun,String)
    for hii in 0 : 0.1 : 1
        expname = "H$(@sprintf("%3.1f",hii))F1d0"
        expname = replace(expname,"."=>"d")

        for ensembleii in 1 : 15

            ensname = "ensemble$(@sprintf("%02d",ensembleii))"
            nrun = projectdir("run","DCS",expname,"$(ensname).sh")

            open(nrun,"w") do wrun
                sn = replace(s ,"[email]"=>"")
                sn = replace(sn,"[user]"=>"")
                sn = replace(sn,"[project]"=>"ExploreWTGSpace")
                sn = replace(sn,"[experiment]"=>"$(expname)")
                sn = replace(sn,"[config]"=>"")
                if ensembleii < 6
                    sn = replace(sn,"[sndname]"=>"S1284km300V64")
                elseif (ensembleii > 5) && (ensembleii < 11)
                    sn = replace(sn,"[sndname]"=>"S1284km300V64_hot")
                else
                    sn = replace(sn,"[sndname]"=>"S1284km300V64_cld")
                end
                sn = replace(sn,"[lsfname]"=>"noforcing")
                sn = replace(sn,"[schname]"=>"DCS")
                sn = replace(sn,"member[xx]"=>"member$(@sprintf("%02d",ensembleii))")
                write(wrun,sn)
            end

        end

    end
    for fii in 0 : 0.1 : 0.9
        expname = "H1d0F$(@sprintf("%3.1f",fii))"
        expname = replace(expname,"."=>"d")

        for ensembleii in 1 : 15

            ensname = "ensemble$(@sprintf("%02d",ensembleii))"
            nrun = projectdir("run","DCS",expname,"$(ensname).sh")

            open(nrun,"w") do wrun
                sn = replace(s ,"[email]"=>"")
                sn = replace(sn,"[user]"=>"")
                sn = replace(sn,"[project]"=>"ExploreWTGSpace")
                sn = replace(sn,"[experiment]"=>"$(expname)")
                sn = replace(sn,"[config]"=>"")
                if ensembleii < 6
                    sn = replace(sn,"[sndname]"=>"S1284km300V64")
                elseif (ensembleii > 5) && (ensembleii < 11)
                    sn = replace(sn,"[sndname]"=>"S1284km300V64_hot")
                else
                    sn = replace(sn,"[sndname]"=>"S1284km300V64_cld")
                end
                sn = replace(sn,"[lsfname]"=>"noforcing")
                sn = replace(sn,"[schname]"=>"DCS")
                sn = replace(sn,"member[xx]"=>"member$(@sprintf("%02d",ensembleii))")
                write(wrun,sn)
            end

        end

    end
end

open(brun,"r") do frun
    s = read(frun,String)
    for hii in 0 : 0.1 : 1
        expname = "H$(@sprintf("%3.1f",hii))F1d0"
        expname = replace(expname,"."=>"d")
        nrun = projectdir("run","DCS",expname,"Build.csh")

        open(nrun,"w") do wrun
            sn = replace(s ,"[user]"=>"")
            sn = replace(sn,"[schname]"=>"DCS")
            sn = replace(sn,"[expname]"=>"$(expname)")
            sn = replace(sn,"[runname]"=>"")
            write(wrun,sn)
        end
    end
    for fii in 0 : 0.1 : 0.9
        expname = "H1d0F$(@sprintf("%3.1f",fii))"
        expname = replace(expname,"."=>"d")
        nrun = projectdir("run","DCS",expname,"Build.csh")

        open(nrun,"w") do wrun
            sn = replace(s ,"[user]"=>"")
            sn = replace(sn,"[schname]"=>"DCS")
            sn = replace(sn,"[expname]"=>"$(expname)")
            sn = replace(sn,"[runname]"=>"")
            write(wrun,sn)
        end
    end
end