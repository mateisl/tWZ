file = open("mvamlp.sh","w")
file.write("cmssw python train_TWZ_3l.py --mva mlp --variables randZWdphilnonZfwdjet --layer Np1 --sampling 1 --epoch 1 --overwrite \n")

val  = [0.5,1]
rang = [-5,-3,-1,1,3,5,7]
ran  = [-3,-1,1,5]
#rang = range(-5,5)
#rang.remove(0) 
for i in val: 
    for j in val: 
        file.write("cmssw python train_TWZ_3l.py --mva mlp --variables randZWdphilnonZfwdjet --layer N --sampling "+str(i)+" --epoch "+str(j)+"\n")
        file.write("cmssw python train_TWZ_3l.py --mva mlp --variables randZWdphilnonZfwdjet --layer NcN --sampling "+str(i)+" --epoch "+str(j)+"\n")
        for layerone in rang : 
            if layerone > 0: layerone = "p"+str(layerone)
            firstlayer = "N"+str(layerone)
            onelayer   = firstlayer.replace("-","m") 
            file.write("cmssw python train_TWZ_3l.py --mva mlp --variables randZWdphilnonZfwdjet --layer "+onelayer+" --sampling "+str(i)+" --epoch "+str(j)+"\n")
            file.write("cmssw python train_TWZ_3l.py --mva mlp --variables randZWdphilnonZfwdjet --layer "+onelayer+"cN --sampling "+str(i)+" --epoch "+str(j)+"\n")
            for layertwo in ran :
                if layertwo > 0: layertwo = "p"+str(layertwo)
                secondlayer = "N"+str(layertwo)
                twolayers   = secondlayer.replace("-","m") 
                file.write("cmssw python train_TWZ_3l.py --mva mlp --variables randZWdphilnonZfwdjet --layer Nc"+twolayers+" --sampling "+str(i)+" --epoch "+str(j)+"\n")
                file.write("cmssw python train_TWZ_3l.py --mva mlp --variables randZWdphilnonZfwdjet --layer "+onelayer+"c"+twolayers+" --sampling "+str(i)+" --epoch "+str(j)+"\n")
                for layerthree in range(1,2):
                    thirdlayer = str(layerthree)
                    alllayers = onelayer + "c" + twolayers +"c"+ thirdlayer

                    file.write("cmssw python train_TWZ_3l.py --mva mlp --variables randZWdphilnonZfwdjet --layer NcNc"+thirdlayer+" --sampling "+str(i)+" --epoch "+str(j)+"\n")
                    file.write("cmssw python train_TWZ_3l.py --mva mlp --variables randZWdphilnonZfwdjet --layer "+onelayer+"cNc"+thirdlayer+" --sampling "+str(i)+" --epoch "+str(j)+"\n")
                    file.write("cmssw python train_TWZ_3l.py --mva mlp --variables randZWdphilnonZfwdjet --layer Nc"+twolayers+"c"+thirdlayer+" --sampling "+str(i)+" --epoch "+str(j)+"\n")
                    file.write("cmssw python train_TWZ_3l.py --mva mlp --variables randZWdphilnonZfwdjet --layer "+alllayers+" --sampling "+str(i)+" --epoch "+str(j)+"\n")         
        
file.close()
