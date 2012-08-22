""" Parses the output of Blast2N and outputs the differences between the two sequences, which can be easily read and used to convert a vcf
"""

import sys, pickle

if __name__ == '__main__':
    f=open(sys.argv[1],'rU')
    wdata=open(sys.argv[1]+'.changes','w')
    wshifts=open(sys.argv[1]+'.shifts','w')
    # Dictionary of positions and their bases
    pos={}
    diff = {}
    shifts={}
    shift=0
    old_margin = 0
    new_margin = 0
    old=[]
    changes=[]
    new=[]

    while True:
        try:
            old_line = f.readline().split('  ')[-2].lstrip(' ').rstrip(' ')
        except:
            break
        old.extend(old_line)
        changes.extend(f.readline()[14:14+len(old_line)])#Can't split becasue it's possible for there to be two spaces in the content, which we don't want to split
        new.extend(f.readline().split('  ')[-2].lstrip(' ').rstrip(' '))
        _ = f.readline()
    
    for n, c in enumerate(changes):
        if c != '|':
            if old[n] == '-':
                old_margin -= 1
                shift+=1
                shifts[n+1+old_margin]=shift
                if old[n+1]=='-':
                    pass
                else:
                    i=-1
                    while old[n+i] == '-':
                        i-=1
                    oldbases = old[n+i]
                    newbases = ''.join([new[n+x] for x in range(i,1)])
                    pos[n + old_margin] = n + 1 + new_margin
                    diff[n + old_margin] = [oldbases.upper(), newbases.upper()]
            elif new[n] == '-':
                new_margin -= 1
                shift-=1
                shifts[n+1+old_margin]=shift
                if new[n+1]=='-':
                    pass
                else:
                    i=-1
                    while new[n+i] == '-':
                        i-=1
                    newbases = new[n+i]
                    oldbases = ''.join([old[n+x] for x in range(i,1)])
                    pos[n + old_margin] = n + 1 + new_margin
                    diff[n + old_margin] = [oldbases.upper(), newbases.upper()]
            else:
                pos[n + 1 + old_margin] = n + 1 + new_margin
                diff[n + 1 + old_margin] = [old[n].upper(), new[n].upper()]

    #Outputs positions where indels affect the positioning between the two sequences
    for key, value in sorted(shifts.items()):
        wshifts.write(str(key)+"\t"+str(value)+"\n")

    #Outputs old position, new position, old ref, new ref
    for key, value in sorted(pos.items()):
        wdata.write(str(key) + "\t" + str(value) + "\t" + "\t".join(diff[key])+"\n") 


