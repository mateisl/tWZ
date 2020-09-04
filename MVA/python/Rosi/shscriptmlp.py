
pen("mvamlp.sh","w")

variables = ["randZWdphidRnonZjet0", "randZWdphijet0", "randZWdphidRnonZ"]#, "randZWdphi"] , "highrankjet0bjet", "highrankjet01bjet", "tommy", "highrankbjetRnonZ", "randZWdphilnonZfwdjet", "randZWdphifwdjet", "highrankjet0btagbjet"]
val  = [0.3,0.5,0.6,0.8,1]
rang = [5, 7] 
ran  = [-1,-3,5]

for var in variables: 
    file.write("cmssw python train_TWZ_3l.py --mva mlp --variables "+var+" --layer Np1 --sampling 1 --epoch 1 --overwrite \n")
    for i in val: 
        for j in val: 
            file.write("cmssw python train_TWZ_3l.py --mva mlp --variables "+var+" --layer Np5 --sampling "+str(i)+" --epoch "+str(j)+"\n")
            file.write("cmssw python train_TWZ_3l.py --mva mlp --variables "+var+" --layer Np7 --sampling "+str(i)+" --epoch "+str(j)+"\n")
            file.write("cmssw python train_TWZ_3l.py --mva mlp --variables "+var+" --layer NcN --sampling "+str(i)+" --epoch "+str(j)+"\n")
            file.write("cmssw python train_TWZ_3l.py --mva mlp --variables "+var+" --layer NcNm1 --sampling "+str(i)+" --epoch "+str(j)+"\n")
            file.write("cmssw python train_TWZ_3l.py --mva mlp --variables "+var+" --layer NcNm1c1 --sampling "+str(i)+" --epoch "+str(j)+"\n")
            file.write("cmssw python train_TWZ_3l.py --mva mlp --variables "+var+" --layer NcNm3 --sampling "+str(i)+" --epoch "+str(j)+"\n")
            file.write("cmssw python train_TWZ_3l.py --mva mlp --variables "+var+" --layer NcNp5c1 --sampling "+str(i)+" --epoch "+str(j)+"\n")
            file.write("cmssw python train_TWZ_3l.py --mva mlp --variables "+var+" --layer NcNp5 --sampling "+str(i)+" --epoch "+str(j)+"\n")
#            for layerone in rang : 
#                if layerone > 0: layerone = "p"+str(layerone)
#                firstlayer = "N"+str(layerone)
#                onelayer   = firstlayer.replace("-","m") 
#                file.write("cmssw python train_TWZ_3l.py --mva mlp --variables "+var+" --layer "+onelayer+" --sampling "+str(i)+" --epoch "+str(j)+"\n")
#                file.write("cmssw python train_TWZ_3l.py --mva mlp --variables "+var+" --layer "+onelayer+"cN --sampling "+str(i)+" --epoch "+str(j)+"\n")
#                for layertwo in ran :
#                    if layertwo > 0: layertwo = "p"+str(layertwo)
#                    secondlayer = "N"+str(layertwo)
#                    twolayers   = secondlayer.replace("-","m") 
#                    file.write("cmssw python train_TWZ_3l.py --mva mlp --variables "+var+" --layer Nc"+twolayers+" --sampling "+str(i)+" --epoch "+str(j)+"\n")
#                    file.write("cmssw python train_TWZ_3l.py --mva mlp --variables "+var+" --layer "+onelayer+"c"+twolayers+" --sampling "+str(i)+" --epoch "+str(j)+"\n")
#                    for layerthree in range(1): #,2):
#                        thirdlayer = str(layerthree)
#                        alllayers = onelayer + "c" + twolayers +"c"+ thirdlayer
#    
#                        file.write("cmssw python train_TWZ_3l.py --mva mlp --variables "+var+" --layer NcNc"+thirdlayer+" --sampling "+str(i)+" --epoch "+str(j)+"\n")
#                        file.write("cmssw python train_TWZ_3l.py --mva mlp --variables "+var+" --layer "+onelayer+"cNc"+thirdlayer+" --sampling "+str(i)+" --epoch "+str(j)+"\n")
#                        file.write("cmssw python train_TWZ_3l.py --mva mlp --variables "+var+" --layer Nc"+twolayers+"c"+thirdlayer+" --sampling "+str(i)+" --epoch "+str(j)+"\n")
#                        file.write("cmssw python train_TWZ_3l.py --mva mlp --variables "+var+" --layer "+alllayers+" --sampling "+str(i)+" --epoch "+str(j)+"\n")         
#            
file.close()
