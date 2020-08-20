file = open("mvamlp.sh","w")
file.write("cmssw python train_TWZ_3l.py --mva mlp --variables randZWdphilnonZ --layer Np1 --sampling 1 --epoch 1 --overwrite \n")

val  = [0.2,0.5,0.8,1]
rang = range(-5,5)
rang.remove(0) 
for layerone in rang : 
    if layerone > 0: layerone = "p"+str(layerone)
    firstlayer = "N"+str(layerone)
    onelayer   = firstlayer.replace("-","m") 
    for layertwo in rang :
        if layertwo > 0: layertwo = "p"+str(layertwo)
        secondlayer = "N"+str(layertwo)
        twolayers   = secondlayer.replace("-","m") 
        for layerthree in range(1,5):
            thirdlayer = str(layerthree)
            alllayers = onelayer + "c" + twolayers +"c"+ thirdlayer
            for i in val:
                for j in val: 
                    file.write("cmssw python train_TWZ_3l.py --mva mlp --variables randZWdphilnonZ --layer N --sampling "+str(i)+" --epoch "+str(j)+"\n")
                    file.write("cmssw python train_TWZ_3l.py --mva mlp --variables randZWdphilnonZ --layer "+onelayer+" --sampling "+str(i)+" --epoch "+str(j)+"\n")

                    file.write("cmssw python train_TWZ_3l.py --mva mlp --variables randZWdphilnonZ --layer Nc"+twolayers+" --sampling "+str(i)+" --epoch "+str(j)+"\n")
                    file.write("cmssw python train_TWZ_3l.py --mva mlp --variables randZWdphilnonZ --layer NcN --sampling "+str(i)+" --epoch "+str(j)+"\n")
                    file.write("cmssw python train_TWZ_3l.py --mva mlp --variables randZWdphilnonZ --layer "+onelayer+"cN --sampling "+str(i)+" --epoch "+str(j)+"\n")
                    file.write("cmssw python train_TWZ_3l.py --mva mlp --variables randZWdphilnonZ --layer "+onelayer+"c"+twolayers+" --sampling "+str(i)+" --epoch "+str(j)+"\n")

                    file.write("cmssw python train_TWZ_3l.py --mva mlp --variables randZWdphilnonZ --layer NcNc"+thirdlayer+" --sampling "+str(i)+" --epoch "+str(j)+"\n")
                    file.write("cmssw python train_TWZ_3l.py --mva mlp --variables randZWdphilnonZ --layer "+onelayer+"cNc"+thirdlayer+" --sampling "+str(i)+" --epoch "+str(j)+"\n")
                    file.write("cmssw python train_TWZ_3l.py --mva mlp --variables randZWdphilnonZ --layer Nc"+twolayers+"c"+thirdlayer+" --sampling "+str(i)+" --epoch "+str(j)+"\n")
                    file.write("cmssw python train_TWZ_3l.py --mva mlp --variables randZWdphilnonZ --layer "+alllayers+" --sampling "+str(i)+" --epoch "+str(j)+"\n")         

file.close()
