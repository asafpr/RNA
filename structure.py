'''
compute the accessibility of nucleotides in a sequence 
'''
import random
import string
import os,sys

char_set = string.ascii_uppercase + string.digits

def average(vlist,winlen,totlen):
    '''
    computes the average of each segmant in the length of winlen
    '''
    retlist=[]
    for i in range(totlen-winlen+1):
        tot=0
        for j in range(i,i+winlen):
            if j in vlist: tot+=vlist[j]
        retlist.append(tot/float(winlen))
    return retlist

def computeACC(seq,winlen=10,temp=37,RNAfold='RNAfold'):
    '''
    Use RNAfold -p to caculate the ensemble of secondary structures then
    open the dot.ps and read these values, sum them for each nucleotide
    and takes a 1-this number to have the accessibility
    '''
    name=''.join(random.sample(char_set,6))
    cmd = 'echo ">'+name+'\n'+seq+'" | '+RNAfold+' -p -T '+str(temp)+' >& /dev/null'
    os.system(cmd)
    probs={}
    #name of file should be name_dp.ps
    try:
        fin = open(name+'_dp.ps')
    except:
        return None
    for line in fin:
        if line.rstrip()[-4:]=='ubox' and line[0]!='%':
            pos1,pos2,prob=line.split(' ')[0:3]
            p1=int(pos1)-1
            p2=int(pos2)-1
            if p1 not in probs: probs[p1]=0
            if p2 not in probs: probs[p2]=0            
            probs[p1]+=float(prob)**2
            probs[p2]+=float(prob)**2
    fin.close()
    os.unlink(name+'_dp.ps')
    os.unlink(name+'_ss.ps')
    avprobs = average(probs,winlen,len(seq))
    for i in range(len(avprobs)):
        avprobs[i]=1-avprobs[i]
    return avprobs

    
    
def main(argv):
    import getopt
    optlist,argv = getopt.getopt(argv,'w:t:')
    winlen=10
    temp=37
    for opt,val in optlist:
        if opt =='-w':
            winlen=int(val)
        elif opt=='-t':
            temp=int(val)
    print computeACC(argv[0],winlen=winlen,temp=temp)


if __name__=='__main__':
    main(sys.argv[1:])

    
