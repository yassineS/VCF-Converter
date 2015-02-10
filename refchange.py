"""This script convert a vcf file variants coordinates from one reference to another
"""
import sys
#Command should look like: python refchange.py <blast file> <vcf file>

# COORD_CHANGES for hg19 to rsrs
COORD_CHANGES = {}

# Define the bad loci here
BAD_LOCI = []

def change_coord(pos, coord_map = COORD_CHANGES):
    change = 0
    for key, value in sorted(coord_map.items()):
        if key < pos:
            change = value
        else: break
    return change


if __name__ == '__main__':
    c =  open(sys.argv[1]+'.changes','rU')
    s =  open(sys.argv[1]+'.shifts','rU')
    vcf= open(sys.argv[2],'rU')
    w =  open(sys.argv[2]+'.conv', 'w')
    
    # Makes COORD_CHANGES
    for line in s:
        line = line.strip().split('\t')
        COORD_CHANGES[int(line[0])]=int(line[1])

    d={}
    # Parse the readable blast format
    for line in c:
        line=line.strip().split('\t')
        d[int(line[0])]=[line[1],line[2],line[3]]

    debug = 0

    changes=[]

    def new_snp(pos, numppl, d):
        s='chrM'#TODO: find a way to figure out from vcf
        s+='\t'+str(d[pos][0])
        s+='\t.'
        s+='\t'+d[pos][2]+'\t'+d[pos][1]
        s+='\t.\t.\t.\t.'
        for i in range(numppl):
            s+='\t1/1'
        s+='\n'
        return s

    def new_snp_pos(d,c):
        l=sorted(list(d.keys()))
        for i in l:
            if i not in c:
                return i
        return 'Done'

    for line in vcf:
        if line[0] == '#':
            if line.split('\t')[0]=='#CHROM':
                numppl=len(line.split('\t')[9:])
            w.write(line)
        else:
            line=line.strip().split('\t')
            pos=int(line[1])
            while(new_snp_pos(d,changes) < pos):
                if d[new_snp_pos(d,changes)][2] != '-':
                    w.write(new_snp(new_snp_pos(d,changes),numppl,d))
                changes.append(new_snp_pos(d,changes))

            newpos=int(line[1])+change_coord(int(line[1]))
            if pos in d.keys():
                changes.append(pos)
                oldref=d[int(line[1])][1].upper()
                newref=d[int(line[1])][2].upper()
                ppl=[]
                if newref == line[4].upper():
                    alt=oldref
                    for p in line[9:]:
                        p=p.split(':')
                        freqs=p[1].split(',')[1]+','+p[1].split(',')[0]
                        if p[0] == '0/0':
                            ppl.append('1/1:'+freqs+':'+':'.join(p[2:]))
                        elif p[0] == '1/1':
                            ppl.append('0/0:'+freqs+':'+':'.join(p[2:]))
                        else:#1/0 or ./., neither of which need changing
                            ppl.append(':'.join(p))
                elif len(line[4].split(','))==1:
                    alt=line[4]+','+oldref
                    for p in line[9:]:
                        p=p.split(':')
                        freqs='0,'+p[1].split(',')[1]+','+p[1].split(',')[0]
                        if p[0] == '0/0':
                            ppl.append('2/2:'+freqs+':'+':'.join(p[2:]))
                        elif p[0] in ['1/0','0/1']:
                            ppl.append('1/2:'+freqs+':'+':'.join(p[2:]))
                        else:#1/1 or ./., neither of which need changing
                            ppl.append(':'.join(p))
                elif newref in line[4].split(','):
                    switch = line[4].split(',').index(newref)
                    alt = ','.join([oldref if x==newref else x for x in line[4].split(',')])
                    def switch_indexes(n,s):
                        if n == s:
                            return 0
                        elif n == 0:
                            return s
                        return n
                    for p in line[9:]:
                        p=p.split(':')
                        freqs=','.join([p[1].split(',')[switch_indexes(n,switch+1)] for n,x in enumerate(p[1].split(','))])
                        f='/'.join([str(switch_indexes(x,switch+1)) for x in p[0].split('/')])
                        ppl.append(f+':'+freqs+':'+':'.join(p[2:]))
                else:
                    alt=line[4]+','+oldref
                    for p in line[9:]:
                        p=p.split(':')
                        freqs='0,'+p[1].split(',')[1]+','+p[1].split(',')[2]+','+p[1].split(',')[0]
                        f='/'.join([x if x!='0' else '3' for x in p[0].split('/')])
                        ppl.append(f+':'+freqs+':'+':'.join(p[2:]))

                w.write(line[0]+'\t'+str(newpos)+'\t'+line[2]+'\t'+newref+'\t'+alt+'\t'+'\t'.join(line[5:9])+'\t'+'\t'.join(ppl)+'\n')
    while new_snp_pos(d,changes) != 'Done':
        w.write(new_snp(new_snp_pos(d,changes), numppl,d))
        changes.append(new_snp_pos(d,changes))
