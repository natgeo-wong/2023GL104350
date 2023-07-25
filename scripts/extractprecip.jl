using DrWatson
@quickactivate "2023GL104350"

include(srcdir("extractprecip.jl"))

schname = "DGW"
expname = "P1282km300V64"

DGW1 = [0,0.2,0.5,1,2,5,10,20,50,100,200,500]
DGW2 = [0,0.02,0.05,0.1,0.2,0.5,1,2,5,10,20,50,100,200,500]
TGR1 = [
    0,sqrt(2),2,2*sqrt(2.5),5,5*sqrt(2),
    10,10*sqrt(2),20,20*sqrt(2.5),50,50*sqrt(2),
]
TGR2 = vcat(TGR1/100,1,TGR1)

pwrvec  = DGW1

for pwr in pwrvec
    extractprecip(schname,expname,dampingstrprnt(pwr))
end
