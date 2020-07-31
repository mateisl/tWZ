file = open("mva.sh","w")
file.write("cmssw python train_TWZ_3l.py --mva bdt --variables randZWdphilnonZ --NTrees 25 --maxdepth 1 --ncuts 50 --overwrite \n")
for i in range(1,4): 
    NTrees = 25
    for j in range(5):
        ncuts = 50 
        for k in range(4): 
            file.write("cmssw python train_TWZ_3l.py --mva bdt --variables randZWdphilnonZ --NTrees "+str(NTrees)+" --maxdepth "+str(i)+" --ncuts "+str(ncuts)+"\n") 
            ncuts *= 2        
        NTrees *= 2 

file.close()
