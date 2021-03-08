import ROOT

def singleton(class_):
  instances = {}
  def getinstance(*args, **kwargs):
    if class_ not in instances:
        instances[class_] = class_(*args, **kwargs)
    return instances[class_]
  return getinstance

@singleton
class color():
  pass

color.data           = ROOT.kBlack
color.DY             = ROOT.kBlue-9
color.TTJets         = 7
color.nonprompt      = ROOT.kBlue-9
color.singleTop      = 40
color.TTX_rare       = ROOT.kRed-10
color.TTH            = ROOT.kRed-10
color.TWZ            = ROOT.kRed
color.TTZ            = 91 #ROOT.kBlack
color.TTW            = ROOT.kOrange + 4
color.TTTT           = ROOT.kGreen + 4
color.TTZtoLLNuNu    = 91
color.signal         = ROOT.kOrange
color.TTZtoQQ        = 91
color.TTG            = ROOT.kRed
color.TTG_signal     = 91
color.TZQ            = ROOT.kRed - 9
color.WJetsToLNu     = ROOT.kRed-10
color.WJets          = ROOT.kRed-10
color.diBoson        = ROOT.kOrange
color.multiBoson     = ROOT.kOrange
color.ZZ             = ROOT.kGreen+3
color.WZ             = ROOT.kAzure+6 #51
color.WZ_amc         = ROOT.kAzure+6 #51
color.WW             = ROOT.kOrange-7
color.VV             = 30
color.WG             = ROOT.kOrange-5
color.triBoson       = ROOT.kOrange+1
color.rare_noZZ      = 8
color.WZZ            = ROOT.kYellow
color.WWG            = ROOT.kYellow-5
color.QCD            = 46
