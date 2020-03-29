import numpy as np
import matplotlib.pyplot as plt

def createChart(numbers, labels, outfile):
    y_pos = np.arange(len(labels))
    colors=[]
    # create a color palette
    palette = plt.get_cmap('Set1')
    for i, label in enumerate(labels):
        colors.append(palette(i))
    bars = plt.bar(y_pos, numbers, color=colors)
    plt.xticks(y_pos, labels, color='blue',rotation=90)
    for rect in bars:
        height = rect.get_height()
        txt = "{:.2f}".format(float((height/numbers[0])*100))+"%"
        #place text above the bar
        plt.text(rect.get_x() + rect.get_width()/2.0, height, txt, ha='center', va='bottom',rotation=45)
    plt.tight_layout()
    plt.savefig(outfile)

def createChartPrc(numbers, labels, prc, outfile):
    y_pos = np.arange(len(labels))
    colors=[]
    prcs=[]
    # create a color palette
    palette = plt.get_cmap('Set1')
    for i, label in enumerate(labels):
        colors.append(palette(i))
    bars = plt.bar(y_pos, numbers, color=colors)
    plt.xticks(y_pos, labels, color='blue',rotation=90)
    for i, rect in enumerate(bars):
        height = rect.get_height()
        #txt = "{:.2f}".format(float((height/numbers[0])*100))+"%"
        #place text above the bar
        plt.text(rect.get_x() + rect.get_width()/2.0, height, prc[i], ha='center', va='bottom',rotation=45)
    plt.tight_layout()
    plt.savefig(outfile)
