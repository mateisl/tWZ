file = open("mvabdt.sh","w")
variables = ["randZWdphidRnonZjet0bjet", "randZWdphibjet", "randZWdphijet", "randZWdphijet0", "randZWdphidR", "randZWdphidR", "randZWdphidRnonZ", "randZWdphi", "highrankjet0bjet", "highrankjet01bjet", "tommy", "highrankbjetRnonZ", "randZWdphilnonZfwdjet", "randZWdphifwdjet", "highrankjet0btagbjet"]

for var in variables: 
    file.write("cmssw python train_TWZ_3l.py --mva bdt --variables "+var+" --NTrees 450 --maxdepth 8 --ncuts 25 --overwrite \n")
    
    ntrees = 100  
    
    for i in range(5):
        ncuts = 10
        for j in range(5):
            for k in range(1,4):
                file.write("cmssw python train_TWZ_3l.py --mva bdt --variables "+var+" --NTrees "+str(ntrees)+" --maxdepth "+str(k)+" --ncuts "+str(ncuts)+"\n")
            ncuts = ncuts*2
        ntrees = (ntrees*10) / 2
             
file.close()
