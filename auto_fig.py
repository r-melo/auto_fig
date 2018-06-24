from difflib import SequenceMatcher
import matplotlib.pyplot as plt
import matplotlib.pylab as pl
import numpy as np
import glob
import csv

def prepara(vec,ordem):
    for i in range(len(vec)):
        vec[i] = [vec[i][j] for j in ordem]
    for i,col in enumerate(vec[2:]):
        vec[i+2] = list(map(float,col))
    vec[1] = list(map(int,vec[1]))
    return vec

def fazfig(name):

    file = []
    with open('averageRelative_' + name, newline='') as csvfile:
        reader = csv.reader(csvfile, delimiter='\t')
        for row in reader:
            file += [row]
    cab = file[0]
    data = file[1:]
    file = []
    with open('pvalue_' + name, newline='') as csvfile:
        reader = csv.reader(csvfile, delimiter='\t')
        for row in reader:
            file += [row]
    pvnames = file[0]
    pval = file[1:]
    file = []
    with open(name, newline='') as csvfile:
        reader = csv.reader(csvfile, delimiter='\t')
        for row in reader:
            file += [row]
    trans = file[2:]
    #Proteins = cab[0]
    #position = cab[1]
    catsize = []
    catnames = []
    itercab = iter(cab[2:])
    for indcat in range(len(cab[2:])//2):
        cat = next(itercab)
        catname = cat.split()
        catsize += [len(catname) -2]
        if len(catname) > 3:
            string1 = catname[1]
            string2 = catname[2]
            match = SequenceMatcher(None, string1, string2).find_longest_match(1, len(string1)-1, 1, len(string2)-1)
            catnames += [string1[match.a: match.a + match.size]]
        else:
            catnames += catname[1][1:-1]
        cat = next(itercab)
        
    datacols = [list(x) for x in zip(*data)]
    pvalcols = [list(x) for x in zip(*pval)]
    transcols = [list(x) for x in zip(*trans)]
    
    ordenado = list(map(int,datacols[1]))
    ordenado.sort()
    ordem = []
    for val in ordenado:
        ordem += [i for i,x in enumerate(datacols[1]) if int(x) == val]

    transcols = prepara(transcols, ordem)
    pvalcols = prepara(pvalcols, ordem)
    datacols = prepara(datacols, ordem)
    for i, col in enumerate(datacols[3::2]):
        datacols[2 * (i+1) + 1] = list(map(lambda x: x/np.sqrt(catsize[i]), col))

    x = datacols[1]
    proteins = datacols[0]
    datacols = datacols[2:]
    transcols = transcols[2:]
    pvalcols = pvalcols[2:]
    
    n = len(catnames)
    colors = pl.cm.nipy_spectral(np.linspace(0,1,n))

    plt.figure(figsize=(14,14))
    plt.subplot(3,1,1)
    for i,name in enumerate(catnames):
        for j in range(catsize[i]):
            soma = sum(catsize[:i])
            #print(i, 2+soma+j,len(transcols))
            plt.plot(x,transcols[soma+j],label=catnames[i],linewidth=0.5, color=colors[i])
            plt.xticks([])
    plt.xlim((1,x[-1]))
    plt.ylabel('Transcriptograms', fontsize=14)
    plt.yticks(fontsize=12)
    plt.subplot(3,1,2)
    for i, col in enumerate(datacols[::2]):
        plt.plot(x,col,label=catnames[i],linewidth=0.5, color=colors[i])
        plt.fill_between(x, [x - y for x, y in zip(col,datacols[2*(i+1)-1])],
                         [x + y for x, y in zip(col,datacols[2*(i+1)-1])], alpha=0.4, color=colors[i])
    plt.xticks([])
    plt.ylabel('Averages', fontsize=14)
    plt.yticks(fontsize=12)
    plt.xlim((1,x[-1]))
    plt.legend(bbox_to_anchor=(1.25, 2.03), fontsize=14)
    plt.subplot(3,1,3)
    for i, col in enumerate(pvalcols):
        plt.plot(x,col,'k',linewidth=0.5)
    plt.xlim((1,x[-1]))
    plt.yscale('log')
    plt.ylabel('p-value', fontsize=14)
    plt.yticks(fontsize=12)
    plt.xticks(fontsize=12)
    plt.subplots_adjust(hspace=0.01)
    plt.savefig(name[:-3] + 'png', dpi=600)
#    plt.show()


for name in glob.glob('raio30*txt'):
    fazfig(name)

