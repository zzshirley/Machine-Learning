#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 20 15:25:02 2018

@author: xiaotong
"""

from __future__ import division
from itertools import *
import xlrd

def loadData(url,sequences):
    data=xlrd.open_workbook(url)
    #sqlfile = open("data_test.txt", "w")#将读取的数据写进文件中
    sh=data.sheet_by_index(0)#获得excels的工作单元
    numbers={'start':'a','checkgrade':'b','explaingrade':'c','retest':'d','reselfeva':'e','lookover':'f','ask':'g','desion':'h'}
    num_row=sh.nrows
    for i in range(1,num_row-1):       
        sequences.append(numbers[sh.row(i)[5].value])
    return sequences

def print_dptable(M):
    print( "    ")
    for i in range(len(M)): print("%12s" % ("%d" % i))
    print()
    for j in M[0].keys():
        print ("%7s: " %j,)
        for t in range(len(M)):
            print ("%7s" % ("%.10f" % M[t][j]))
        print()

class HMM:
    
    def __init__(self,obs,state,start_p,trans_p,emit_p):
        
        self.obs=obs
        self.states=states
        self.start_p=start_p
        self.trans_p=trans_p
        self.emit_p=emit_p
    
    def print_txt(self,M0,M1,M2):
        print('trans_p\n',M0)
        print("emit_p\n",M1)
        print("start_p\n",M2)
        sqlfile = open("data_test.txt", "w")#将读取的数据写进文件中
        
        
        sqlfile.writelines("trans_p\n")
        for i in self.states:
            for t in M0[i]:
                #sqlfile.writelines("--")
                sqlfile.write(str(M0[i][t]))#float类型不可以用write/writelines方法
                sqlfile.writelines(" ")
            sqlfile.writelines("\n")
        sqlfile.writelines("\n================")
        
        sqlfile.writelines("\nemit_p\n")
        for i in self.states:
            for t in M1[i]:
                sqlfile.write(str(M1[i][t]))#float类型不可以用write/writelines方法
                sqlfile.writelines(" ")
            sqlfile.writelines("\n")
        sqlfile.writelines("\n================")
        
        sqlfile.writelines("\nstart_p\n")
        for i in self.states:
            #sqlfile.writelines("--")
            sqlfile.write(str(M2[i]))#float类型不可以用write/writelines方法
            sqlfile.writelines(" ")
        sqlfile.writelines("\n")
        
    def print_start_p(self,M):
        print('=============\n',M)
        sqlfile = open("data.txt", "w")#将读取的数据写进文件中
        for i in self.states:
            #sqlfile.writelines("--")
            sqlfile.write(str(M[i]))#float类型不可以用write/writelines方法
            sqlfile.writelines(" ")
        sqlfile.writelines("\n")
        sqlfile.writelines("\n================")
    

    def forward(self):
        self.alpha = [{}]
        self.prob_alpha = 0
        for i in self.states:#计算初始alpha值
            self.alpha[0][i] = self.start_p[i] * self.emit_p[i][self.obs[0]]
        #print(self.alpha[0])
        for t in range(1,len(self.obs)):
            self.alpha.append({})
            for i in self.states:
                self.alpha[t][i]=self.emit_p[i][self.obs[t]]*sum(self.alpha[t-1][j]*self.trans_p[j][i] for j in self.states)
        self.prob_alpha = sum(self.alpha[t][j] for j in self.states)
        
        return (self.alpha, self.prob_alpha)
    
    def backward(self):
        self.beta=[]
        self.prob_beta=0
        for t in range(len(self.obs)):
            self.beta.append({})
            for i in self.states:
                if t==len(self.obs)-1:
                    self.beta[t][i]=1
                else:
                    self.beta[t][i]=0
                    
        for t in reversed(range(len(self.obs))):
            for i in self.states:
                if t in range(len(self.obs)-1):
                    self.beta[t][i] = sum(self.trans_p[i][j] * self.emit_p[j][self.obs[t+1]] *self.beta[t+1][j] for j in self.states)
        
        self.prob_beta=sum(self.start_p[i]*self.emit_p[i][self.obs[0]]*self.beta[0][i] for i in self.states)
        
        return(self.beta,self.prob_beta)
        
    def Gamma(self):
        self.gamma0=[]
        for t in range(len(self.obs)-1):
            self.gamma0.append({})
            for i in self.states:
                self.gamma0[t][i]=(self.alpha[t][i]*self.beta[t][i])/sum(self.alpha[t][j]*self.beta[t][j] for j in self.states)
        return self.gamma0
    
    def Yita(self):
        self.yita0=[]
        for t in range(len(self.obs)-1):
            self.yita0.append({})
            for i in self.states:
                self.yita0[t][i]={}
                for j in self.states:
                    numerator=self.alpha[t][i]*self.trans_p[i][j]*self.emit_p[j][self.obs[t+1]]*self.beta[t+1][j]#分子
                    denominator=sum(sum(self.alpha[t][i]*self.trans_p[i][j]* self.beta[t+1][j]*self.emit_p[j][self.obs[t+1]] for j in self.states) for i in self.states)
                    self.yita0[t][i][j]=numerator/denominator
        return self.yita0
    
    def Baum_Welch(self,iters=50):
        
        def iterations():
            
            new_pi={}
            new_emit={}
            new_trans={}
            
            HMM.Yita()
            HMM.Gamma()
            
            for i in self.states:
                new_pi[i]=round(self.gamma0[0][i],2)
                
            for i in self.states:
                new_trans[i]={}
                for j in self.states:
                    new_trans[i][j]=round(sum(self.yita0[t][i][j] for t in range(len(self.obs)-1))/sum(self.gamma0[t][i] for t in range(len(self.obs)-1)),2)
                    
            new_emit={}
            for j in self.states:
                new_emit[j]={}
                for k in self.obs:
                    num=0
                    den=0
                    for t in range(len(self.obs)-1):
                        if k==self.obs[t]:
                            num+=self.gamma0[t][j]
                        den+=self.gamma0[t][j]
                    new_emit[j][k]=round(num/den,2)
            self.start_p = new_pi
            self.trans_p = new_trans
            self.emit_p = new_emit
            
        for i in range(iters):
            HMM.forward()
            HMM.backward()
            iterations()
            
    def execute(self,*args):
        for arg in args:
            if arg == 'forward':
                print ('Forward algoritms')
                self.forward()
                #self.print_dptable(self.alpha)
                #print (self.prob_alpha, observations)
                #print('\n')
            if arg =='backward':
                print('Backward algoritms')
                self.backward()
                print_dptable(self.beta)
                print (self.prob_beta, observations)
                print('\n')
            if arg == 'Baum_Welch':
                print ('Baum-Welch algoritms')
                self.forward()
                self.backward()
                self.Baum_Welch()
                self.print_txt(self.trans_p,self.emit_p,self.start_p)
                print('\n')


if __name__ == "__main__":
    

    S = "/Users/xiaotong/Documents/新测评系统数据/10-24/10-24协商数据.xlsx"#读取的文件路径
    states = ('S1', 'S2','S3')
    
    sequences= []

    start_probability = {'S1' : 0.4, 'S2' : 0.3,'S3':0.3}

    transition_probability = {
            'S1' : {'S1' : 0.3, 'S2' : 0.3, 'S3' : 0.4},
            'S2' : {'S1' : 0.5, 'S2' : 0.2, 'S3' : 0.3},
            'S3' : {'S1' : 0.33, 'S2' : 0.33, 'S3' : 0.34}
            }

    emission_probability = {
            'S1' : {'a' : 0.1, 'b' : 0.15, 'c' : 0.125, 'd' : 0.05, 'e' : 0.2, 'f' : 0.125, 'g' : 0.125, 'h' : 0.125},
            'S3' : {'a' : 0.15, 'b' : 0.1, 'c' : 0.125, 'd' : 0.125, 'e' : 0.2, 'f' : 0.005, 'g' : 0.125, 'h' : 0.125},
            'S2' : {'a' : 0.225, 'b' : 0.225, 'c' : 0.05, 'd' : 0.1, 'e' : 0.125, 'f' : 0.1, 'g' : 0.05, 'h' : 0.125}
            }

    observations=loadData(S,sequences)
    
    HMM = HMM(observations, states, start_probability, transition_probability, emission_probability)

    HMM.execute('Baum_Welch')
