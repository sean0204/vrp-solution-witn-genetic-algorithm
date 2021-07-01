import random, math, sys
import matplotlib.pyplot as plt # 畫圖
from copy import deepcopy
from tqdm import *  # 進度條
from scipy.spatial.distance import pdist
from scipy.spatial.distance import squareform
import numpy as np
import pandas as pd

DEBUG = False

geneNum = 100 # 種群數量
generationNum = 1500  # 迭代次數

CENTER = 0  # 配送中心

HUGE = 999999
VARY = 0.05  # 變異機率
n = 90  # 客戶點數量
k = 8   # 車輛數量
Q = 264 # 額定載重量, t
q=np.array([128.8,192,192,192,264,264,264,264,264,264,264,264,264,264,264,264,264,416,416])
carC=[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21]
carx=[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21]
dis = 160  # 續航里程, km
costPerKilo = 10  # 油價
epu = 20  # 早到逞罰成本
lpu = 30  # 晚到逞罰成本
speed = 40  # 速度，km/h

Tkmax=480
M=99999


h=[0,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20]
h=[]
h.append(0)
for i in range(90):
    h.append(20)




def getdata(filename):
    
    data=pd.read_csv(filename)
    dArray=data.values[:,1:]
     
    return dArray


def getalldata():
    # a_ik=getdata('a_ik.csv') #車輛能到達之站點
    a_k=getdata('a_k.csv')   #車輛之燃油效率
    c_ij=getdata('c_ij.csv') #配送路徑
    d_i=getdata('d_i.csv')   #站點需求量
    e_i=getdata('e_i.csv')   #各站點最早到達時間
    l_i=getdata('l_i.csv')   #各站點最晚到達時間
    Q_k=getdata('Q_k.csv')   #車容量
    t_ij=getdata('t_ij.csv') #行駛路徑之運輸時間

    a_k=a_k.flatten()
    d_i=d_i.flatten()
    e_i=e_i.flatten()
    l_i=l_i.flatten()
    Q_k=Q_k.flatten()
    
    c_ij=c_ij/1000
    t_ij=t_ij/60
    
    return a_k ,c_ij ,d_i ,e_i ,l_i ,Q_k ,t_ij
a_k ,distB ,t ,eh ,lh ,q ,t_ij =getalldata()

class Gene:                          #初始化
    def __init__(self, name='Gene', data=None):
        self.name = name
        self.length = n + 1
        self.usedcar=[]
        self.subpathT=[]
        self.subpathDist=[]
        self.subpathPos=[]
        self.pathtime=[]
        self.twoTimeIndex=[]
        self.carSubpath=[]
        self.TimeIneed=[]
        self.subPTNeed=[]
        self.genefit=0
        self.pathtimeadd20=[]
        self.costlist=[]
        self.possible=[]
        if data is None:
            self.data = self._getGene(self.length)
        else:
            self.data = data
        self.fit = self.getFit()
        self.chooseProb = 0  # 选择概率
        
        
    # randomly choose a gene
    def _generate(self, length):  #亂碼產生初帶
        data = [i for i in range(1, length)]
        random.shuffle(data)
        data.insert(0, CENTER)
        data.append(CENTER)
        return data

    # insert zeors at proper positions
    def _insertZeros(self, data):  #將亂碼依照需求量插入0
        sum = 0
        newData = []
        usedcar=[]
        # random.shuffle(carx)
        # x=len(carx)-1
        # Q=q[carx[x]]
        global distB
        Q=270
        for index, pos in enumerate(data):
            sum += t[pos]
            if sum > Q:
                newData.append(CENTER)
                sum = t[pos]
                # x-=1
                # Q=q[carx[x]]
            newData.append(pos)
        
        return newData
    
    # return a random gene with proper center assigned
    def _getGene(self, length):
        data = self._generate(length)
        data = self._insertZeros(data)
        # data,usedcar,subpathT,subpathDist,subpathPos = self._insertZeros(data)
        return data
    # ,usedcar,subpathT,subpathDist,subpathPos

    # return fitness
    def getFit(self):#計算適應度
        fit = distCost = timeCost = overloadCost = overTimeCost=pathcost = 0
        culposT = culdist = x =0
        sumtime =sumtimeadd20 =0
        usedcar=[]
        subpathT=[]
        pathDist=[]
        pathPos=[]
        subpathPos=[]
        pathtime=[]
        pathtimeadd20=[]
        global distB
        for i,pos in  enumerate(self.data):
            if i==0:
                pathPos.append(pos)
                continue
            else:
                culposT+=t[pos]
                culdist+=distB[self.data[i-1]][self.data[i]]
                pathPos.append(pos)
                sumtime+=(t_ij[self.data[i-1]][self.data[i]])
                sumtimeadd20+=(t_ij[self.data[i-1]][self.data[i]])+20
                if(self.data[i]==0):
                    subpathT.append(culposT)
                    culposT=t[pos]
                    pathDist.append(culdist)
                    culdist=0
                    subpathPos.append(pathPos)
                    pathPos=[]
                    pathPos.append(0)
                    pathtime.append(sumtime)
                    sumtime=0
                    pathtimeadd20.append(sumtimeadd20)
                    sumtimeadd20=0
        self.subpathT=subpathT
        self.subpathDist=pathDist
        self.subpathPos=subpathPos
        self.pathtime=pathtime 
        self.pathtimeadd20=pathtimeadd20
        

        
     
        timeSpent = 0
        for i, pos in enumerate(self.data):
            # skip first center
            if i == 0:
                continue
            # new car
            elif pos == CENTER:
                timeSpent = 0
            # update time spent on road
            timeSpent += (t_ij[self.data[i-1]][self.data[i]])
            # arrive early
            # if timeSpent < eh[pos]:
            #     timeCost += ((eh[pos] - timeSpent) * epu)
            #     timeSpent = eh[pos]
            # # arrive late
            if timeSpent > lh[pos]:
                timeCost += (timeSpent - lh[pos])
            if timeSpent > 540:
                timeCost+=(timeSpent*HUGE)
            # update time
            timeSpent += h[pos]
        # overload cost and out of fuel cost
        
        choosetime=[]
        sumtwoT=[]
        twoTimeI=[]
     
        choosetime=np.array(pathtimeadd20)
        choosetime.sort
        choosetimeBig=choosetime[choosetime>270]
        choosetimeLess=choosetime[choosetime<=270]
      
        culchTime=0
        pathIndex=[]
        pathTimeindex=[]
        need=[]
        Tneed=[]
        maxT=[]
        subPTindex=[]
        if(len(choosetimeLess))==1:
            twoTimeI.append(choosetimeLess[-1])
            find=pathtimeadd20.index(choosetimeLess[-1])
            pathIndex.append(find)
            pathNeed=subpathT[find]
            need.append(pathNeed)
            pathTimeindex.append(pathIndex)
            maxT.append(pathNeed)
        else:   
            for i in range(len(choosetimeLess)):
                culchTime+=choosetimeLess[i]
                if culchTime<540:
                    sumtwoT.append(choosetimeLess[i])
                    find=pathtimeadd20.index(choosetimeLess[i])
                    pathIndex.append(find)
                    pathNeed=subpathT[find]
                    need.append(pathNeed)
                    if i==len(choosetimeLess)-1:
                        twoTimeI.append(sumtwoT)
                        if len(pathIndex)!=1:
                            pathTimeindex.append(pathIndex)
                            Tneed.append(need)
                            maxT.append(max(need))
                        else:
                            subPTindex.append(pathIndex)
                            pathTimeindex.append(subPTindex)
                            subPTindex=[]
                            Tneed.append(need)
                            maxT.append(max(need))
                elif culchTime>540:
                    twoTimeI.append(sumtwoT)
                    sumtwoT=[]
                    culchTime=choosetimeLess[i]
                    sumtwoT.append(choosetimeLess[i])
                    pathTimeindex.append(pathIndex)
                    Tneed.append(need)
                    maxT.append(max(need))
                    need=[]
                    pathIndex=[]
                    find=pathtimeadd20.index(choosetimeLess[i])
                    pathIndex.append(find)
                    pathNeed=subpathT[find]
                    need.append(pathNeed)
                    if i==len(choosetimeLess)-1:
                        twoTimeI.append(choosetimeLess[i])
                        subPTindex.append(find)
                        pathTimeindex.append(subPTindex)
                        subPTindex=[]
                        Tneed.append(pathNeed)
                        maxT.append(pathNeed)
        for i in range(len(choosetimeBig)):
            twoTimeI.append(choosetimeBig[i])
            find=pathtimeadd20.index(choosetimeBig[i])
            subPTindex.append(find)
            pathTimeindex.append(subPTindex)
            pathNeed=subpathT[find]
            Tneed.append(pathNeed)
            maxT.append(pathNeed)
            subPTindex=[]
        self.twoTimeIndex=twoTimeI
        self.carSubpath=pathTimeindex
        self.TimeIneed=Tneed
        

        self.subPTNeed=maxT
        fitcar=[]
        for i in range(len(maxT)):
            near=(np.abs(q -maxT[i] )).argmin()
            if q[near]<maxT[i]:
                near+=1
                if carC[near]==carC[-1]:
                    fitcar.append(carC[near])
                else:
                    while carC[near] in fitcar:
                            near+=1
                            if carC[near]==carC[-1]:
                                fitcar.append(carC[near])
                                break
                    fitcar.append(carC[near])
                    
            else:
                while carC[near] in fitcar:
                            near+=1
                            if carC[near]==carC[-1]:
                                fitcar.append(carC[near])
                                break
                fitcar.append(carC[near])
        self.usedcar=fitcar
        
        pathcost=0

        for i in range(len(pathTimeindex)):
            carFuel=a_k[fitcar[i]]            
            for j in range(len(pathTimeindex[i])):
                subDist=pathDist[pathTimeindex[i][j]]
                pathcost+=subDist*carFuel*19.2

            

        load = 0
        distAfterCharge = 0
        Q_i=0
        for i, pos in enumerate(self.data):
            # skip first center
            if i == 0:
                continue
            # charge here
            if pos > n:
                distAfterCharge = 0
            # at center, re-load
            elif pos == CENTER:
                load = 0
                distAfterCharge = 0
            # normal
            else:
                load += t[pos]
               
                if Q_i<2:
                    Q=416
                    overloadCost += (HUGE * (load > Q))
                    Q_i+=1
                else:
                    Q=264
                    overloadCost += (HUGE * (load > Q))
        
        carcost=len(fitcar)*966.667

        sumTotaltime=sum(pathtimeadd20)
        fit = pathcost + timeCost + overloadCost + carcost + sumTotaltime
        
        costlist=[]
        costlist.append(fit)
        costlist.append(pathcost)
        costlist.append(sumTotaltime)
        costlist.append(carcost)
        costlist.append(timeCost)
        self.genefit=costlist
        
        
        return 1/fit

    def updateChooseProb(self, sumFit):#更新個體適應度
        self.chooseProb = self.fit / sumFit

    def moveRandSubPathLeft(self):#選擇要交配的一緞子路徑
        path = random.randrange(k-1)  # choose a path index
        index = self.data.index(CENTER, path+1) # move to the chosen index
        # move first CENTER
        locToInsert = 0
        self.data.insert(locToInsert, self.data.pop(index))
        index += 1
        locToInsert += 1
        # move data after CENTER)
        while self.data[index] != CENTER:
            self.data.insert(locToInsert, self.data.pop(index))
            index += 1
            locToInsert += 1
            
    # plot this gene in a new window
    # def plot(self):
    #     Xorder = [X[i] for i in self.data]
    #     Yorder = [Y[i] for i in self.data]
    #     plt.plot(Xorder, Yorder, c='black', zorder=1)
    #     plt.scatter(X, Y, zorder=2)
    #     plt.scatter([X[0]], [Y[0]], marker='o', zorder=3)
    #     plt.title(self.name)
    #     for i in range(28):
    #         plt.annotate(i, (X[i], Y[i]))
    #     plt.show()       


def getSumFit(genes):#計算總適應度
    sum = 0
    for gene in genes:
        sum += gene.fit
    return sum


# return a bunch of random genes
def getRandomGenes(size):
    genes = []
    for i in range(size):
        genes.append(Gene("Gene "+str(i)))
    return genes


# 计算适应度和
# def getSumFit(genes):
#     sumFit = 0
#     for gene in genes:
#         sumFit += gene.fit
#     return sumFit


# 更新选择概率
def updateChooseProb(genes):
    sumFit = getSumFit(genes)
    for gene in genes:
        gene.updateChooseProb(sumFit)

def updateUsedcar(genes):#更新選擇車輛
    for gene in genes:
        gene.chooseCar()

# 计算累计概率
def getSumProb(genes):
    sum = 0
    for gene in genes:
        sum += gene.chooseProb
    return sum


# 选择复制，选择前 1/3
def choose(genes):
    num = int(geneNum/6) * 2    # 选择偶数个，方便下一步交叉
    # sort genes with respect to chooseProb
    key = lambda gene: gene.chooseProb
    genes.sort(reverse=True, key=key)
    # return shuffled top 1/3
    return genes[0:num]


# 交叉一对
def crossPair(gene1, gene2, crossedGenes):
    gene1.moveRandSubPathLeft()
    gene2.moveRandSubPathLeft()

    newGene1 = []
    newGene2 = []
    # copy first paths
    centers = 0
    firstPos1 = 1
    for pos in gene1.data:
        firstPos1 += 1
        centers += (pos == CENTER)
        newGene1.append(pos)
        if centers >= 2:
            break
    centers = 0
    firstPos2 = 1
    for pos in gene2.data:
        firstPos2 += 1
        centers += (pos == CENTER)
        newGene2.append(pos)
        if centers >= 2:
            break
    
    # copy data not exits in father gene
    for pos in gene2.data:
        if pos not in newGene1:
            newGene1.append(pos)
    for pos in gene1.data:
        if pos not in newGene2:
            newGene2.append(pos)
    # add center at end
    newGene1.append(CENTER)
    newGene2.append(CENTER)

    key = lambda gene: gene.fit
    possible = []
    # while gene1.data[firstPos1] != CENTER:
    for i in range(7):
        newGene = newGene1.copy()
        # newGene.insert(firstPos1+random.randint(5,8), CENTER)
        newGene.insert(firstPos1, CENTER)  
        # newGene = Gene(data=newGene.copy())
        possible.append(newGene)
        firstPos1 += 1
        

        
    P0index=[]
    cultime=0
    sumneed=0
    for i in range(len(possible)):
        x=0
        for j in possible[i]:
            if j==0:
              P0index.append(x)
            x+=1
        k=P0index[len(P0index)-2]+1
        lenk=len(possible[i])
        for k in range(lenk):
            cultime += (t_ij[possible[i][k-1]][possible[i][k]])+20
            sumneed+=t[possible[i][k]]
            if possible[i][k]==0:
                cultime=t_ij[0][possible[i][k]]
                sumneed=t[possible[i][k]]
                continue
            elif  cultime>540 or sumneed>270:
                # print('%d============%d'%(cultime,sumneed))
                possible[i].insert(k-1,CENTER)
                cultime=t_ij[0][possible[i][k]]+20
                sumneed=t[possible[i][k]]
                lenk=len(possible[i])
                
        P0index=[]

    P0index=[]
    cultime=0
    sumneed=0
    for i in range(len(possible)):
        for j in possible[i]:
            if j==0:
              P0index.append(x)
        k=P0index[len(P0index)-2]+1
        lenk=len(possible[i])

        for k in range(lenk):
            cultime += (t_ij[possible[i][k-1]][possible[i][k]])+20
            sumneed+=t[possible[i][k]]

            if possible[i][k]==0:
                cultime=t_ij[0][possible[i][k]]
                sumneed=t[possible[i][k]]
                continue
            elif  cultime>500 or sumneed>270:
                # print('%d============%d'%(cultime,sumneed))
                possible[i].insert(k-1,CENTER)
                cultime=t_ij[0][possible[i][k]]+20
                sumneed=t[possible[i][k]]
                
                
        P0index=[]
    for i in range(len(possible)):
        possible[i] =  Gene(data=possible[i].copy())


    
    possible.sort(reverse=True, key=key)
    # print(possible[0].data)
    if possible is None:
        print('bug')

    
    crossedGenes.append(possible[0])
    # key = lambda gene: gene.fit
    possible = []
    # while gene2.data[firstPos2] != CENTER:
    for i in range(7):
        newGene = newGene2.copy()
        # newGene.insert(firstPos2+random.randint(5,8), CENTER)
        newGene.insert(firstPos2, CENTER)
        possible.append(newGene)
        firstPos2 += 1

    
    if possible == None:
        print('bug')
    
    P0index=[]
    cultime=0
    sumneed=0
    for i in range(len(possible)):
        x=0
        for j in possible[i]:
            if j==0:
              P0index.append(x)
            x+=1
        k=P0index[len(P0index)-2]+1
        lenk=len(possible[i])
        for k in range(lenk):
            cultime += (t_ij[possible[i][k-1]][possible[i][k]])+20
            sumneed+=t[possible[i][k]]
            if possible[i][k]==0:
                cultime=t_ij[0][possible[i][k]]
                sumneed=t[possible[i][k]]
                continue
            elif  cultime>540 or sumneed>270:
                # print('%d============%d'%(cultime,sumneed))
                possible[i].insert(k-1,CENTER)
                cultime=t_ij[0][possible[i][k]]+20
                sumneed=t[possible[i][k]]
                lenk=len(possible[i])
      
        P0index=[]
    P0index=[]
    cultime=0
    sumneed=0
    for i in range(len(possible)):
        x=0
        for j in possible[i]:
            if j==0:
              P0index.append(x)
            x+=1
        k=P0index[len(P0index)-2]+1
        lenk=len(possible[i])
        for k in range(lenk):
            cultime += (t_ij[possible[i][k-1]][possible[i][k]])+20
            sumneed+=t[possible[i][k]]
            if possible[i][k]==0:
                cultime=t_ij[0][possible[i][k]]
                sumneed=t[possible[i][k]]
                continue
            elif  cultime>500 or sumneed>270:
                # print('%d============%d'%(cultime,sumneed))
                possible[i].insert(k-1,CENTER)
                cultime=t_ij[0][possible[i][k]]+20
                sumneed=t[possible[i][k]]
                lenk=len(possible[i])
    for i in range(len(possible)):
        possible[i] =  Gene(data=possible[i].copy())
                
        P0index=[]
    key = lambda gene: gene.fit
    possible.sort(reverse=True, key=key)
    a = []
 
     
        


    crossedGenes.append(possible[0])
    
   

# 交叉
def cross(genes):
    crossedGenes = []
    for i in range(0, len(genes), 2):
        crossPair(genes[i], genes[i+1], crossedGenes)  
    # for i in range(len(crossedGenes)):
    #     print('aaaa',crossedGenes[i].data)

    return crossedGenes



# 合并
def mergeGenes(genes, crossedGenes):
    # sort genes with respect to chooseProb
    key = lambda gene: gene.chooseProb
    genes.sort(reverse=True, key=key)
    pos = geneNum - 1
    
    
    for gene in crossedGenes:
        genes[pos] = gene
        pos -= 1
    # print('bbb',genes[-1].data)


    return  genes



# 变异一个
def varyOne(gene):
    varyNum = 10
    variedGenes = []
    for i in range(varyNum):
        p1, p2 = random.choices(list(range(1,len(gene.data)-2)), k=2)
        newGene = gene.data.copy()
        newGene[p1], newGene[p2] = newGene[p2], newGene[p1] # 交换
        variedGenes.append(Gene(data=newGene.copy()))
    key = lambda gene: gene.fit
    variedGenes.sort(reverse=True, key=key)

    return variedGenes[0]


# 变异
def vary(genes):
    for index, gene in enumerate(genes):
        # 精英主义，保留前三十
        if index < 30:
            continue
        if random.random() < VARY:
            genes[index] = varyOne(gene)
    # print('bbb',genes[-1].data)

    return genes

def printGeneData(self):
    totalcost=timecost=fuelcost=0
        
    for i in range(len(self.carSubpath)):
        print('車輛編號',self.usedcar[i])
        print('車輛行走路徑',self.carSubpath[i])
        carFuel=a_k[self.usedcar[i]]
        for j in range(len(self.carSubpath[i])):
            subDist=self.subpathDist[self.carSubpath[i][j]]
            fuelcost+=subDist*carFuel*19.2
            timecost+=self.pathtimeadd20[self.carSubpath[i][j]]
        print('車輛路徑成本',fuelcost)
        print('車輛行駛時間',timecost)
        totalcost=fuelcost+timecost+966.667
        print('車輛總成本',totalcost)
        totalcost=fuelcost=timecost=0


if __name__ == "__main__" and not DEBUG:
    genes = getRandomGenes(geneNum) # 初始种群

    #迭代
    for i in tqdm(range(generationNum)):
        updateChooseProb(genes)
        sumProb = getSumProb(genes)
        chosenGenes = choose(deepcopy(genes))   # 选择
        crossedGenes = cross(chosenGenes)   # 交叉
        genes = mergeGenes(genes, crossedGenes) # 复制交叉至子代种群
        genes = vary(genes) # under construction


    
    # sort genes with respect to chooseProb
    key = lambda gene: gene.fit
    genes.sort(reverse=True, key=key)   # 以fit對種群排序
    print('last',genes[len(genes)-1].data)
    print('\r\n')

    print('fit:', genes[0].genefit)
    
    print('cost',genes[len(genes)-1].genefit)

    

if DEBUG:
    print("START")
    gene = Gene()
    print(gene.data)
    gene.moveRandSubPathLeft()
    print(gene.data)


    print("FINISH")

