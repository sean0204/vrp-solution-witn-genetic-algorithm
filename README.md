# vrp-solution-witn-genetic-algorithm
solve VRP with genetic algorithm based on python

def getfit():

第一個for 
 取基因碼並從起始計算到基因結束
 將此段基因碼計算
 子路徑站點
 子路徑需求量
 子路徑距離
 子路徑時間

第二for 
 計算基因是否符合時窗限制條件
 如果為否將適應值乘極大值
 此基因將在篩選時被淘汰

if(len(choosetimeLess))==1到self.usedcard
 挑選車輛行走路徑
 多段路徑可以使用同一台車運送時
 儲存車輛編號與車輛運送之子路徑

第三FOR
 各別子路徑與車輛運送之路徑
 計算個別車輛燃油運送之子路徑距離的成本

第四FOR
 計算基因子路徑是否符合車輛運送容量
 如果為否將成本乘極大值
 篩選時將被淘汰

def crossPair:
將兩段基因碼進行交配
並且生成多個新子代
從子代之中選擇最好的基因並且更新成本

def crsoo:
取母體前三分之一進行交配生成子代

def mergeGenes:
將產生的子代從母體最後一個進行更新(更新母體後三分之一的基因資料)
