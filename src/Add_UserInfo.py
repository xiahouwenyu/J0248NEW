import ROOT

file = ROOT.TFile.Open("./residual_all.root", "Updates")
for bin in range(6):
    
    tdata=file.Get("nHit%02d"%int(bin)).data
    tbkg=file.Get("nHit%02d"%int(bin)).bkg

    obj1=ROOT.TParameter(int)("Nside",1024)
    obj2=ROOT.TParameter(int)("Scheme",0)

    tdata.GetUserInfo().Add(obj1)
    tdata.GetUserInfo().Add(obj2)
    tbkg.GetUserInfo().Add(obj1)
    tbkg.GetUserInfo().Add(obj2)

file.Write("")
file.Close()