import matplotlib.pyplot as plt
import csv
import sys

"""
Mandatory Arguments...
argv[2] filename := file with csv extension that contains rows of residuals. Will plot this on y axis.
argv[3] first row:= first row in file to plot, (zero-indexed)
argv[4] last row := last row to plot, followed by 

args(list of args):= label for each corresponding row, i.e (float, posit32, bfloat, posit16, etc...),
    len(args) must equal last row - first row. 

Optional Arguments (kwargs)...

traffic    := file with csv extension that contains fraction bit discrepancies.
matrixname := name of matrix
"""
def plot(filename, start, stop, *args, **kwargs):

    start, stop = int(start), int(stop)    

    if ("traffic" in kwargs):
        fig, ax = plt.subplots(2, 1)
        plt.subplots_adjust(hspace=.6)
        cp = ax[0]
        ax[1].set_facecolor('#e6e9ed')
    else:
        fig, ax = plt.subplots(1, 1)
        cp = ax
    cp.set_facecolor('#e6e9ed')

    fig.set_facecolor('#ebeef2')
    
    cp.set_yscale('log')

    pcolors = [('orange', '.'), ('black', '*'), ('blue', '^')]

    if 'matrixname' in kwargs:
        fig.suptitle(kwargs['matrixname'])
    
    with open(filename, 'r') as file:
        for plot in enumerate(csv.reader(file, delimiter=',')):
            tag     = plot[0]
            if tag not in range(start, stop+1): continue
            
            color     = pcolors[(tag % len(pcolors))][0]
            style     = pcolors[(tag % len(pcolors))][1]
            residuals = list(map(float, plot[1]))
            
            cp.plot(range(len(plot[1])), residuals, label=args[tag-start], color=color, marker=None)

        cp.set_title ('Convergence', fontsize=16)
        cp.set_xlabel('Iterations' , fontsize=14)
        cp.set_ylabel('Residual'   , fontsize=14)  

    if 'traffic' in kwargs:
        with open(kwargs['traffic']) as file:
            data          = csv.reader(file, delimiter=',')
            #next(data)
            #next(data)
            discrepancies = list(map(float, next(data)))[13:]
            occurrences   = list(map(float, next(data)))[13:]
            
            ax[1].bar(discrepancies, height=occurrences, color='gray')
            ax[1].set_title ('Extra bits of precision for posit'      )
            ax[1].set_xlabel('Bits')
            ax[1].set_ylabel('Percent of cases'       )
            ax[1].xaxis.set_ticks(discrepancies) #[0:len(discrepancies):2])

    fig.legend(prop={'size':16})
    plt.show()

def plotPrecision(data): 
    print(data)
    fig, ax = plt.subplots(1, 1)
    cp = ax
    cp.set_facecolor('#e6e9ed')

    fig.set_facecolor('#ebeef2')
    
    cp.set_yscale('log')
    cp.set_xscale('log')

    pcolors = [('blue'), ('orange'), ('gray'), ('red'), ('green')]

    fig.suptitle('Precision Comparison for Posit32 vs. Float32')

    for tag, file in enumerate(data):
        with open(file, 'r') as f:
            x, y = [], []
            for line in f:
                x.append(float(line.split(',')[0]))
                y.append(float(line.split(',')[1]))
            cp.plot(x, y, label=data[file], color=pcolors[tag%len(pcolors)], marker=None)

    #cp.set_title ('Precision distribution for Posit32 vs. Float32', fontsize=16)
    cp.set_xlabel('Numerical Value' , fontsize=14)
    cp.set_ylabel('Precision'   , fontsize=14)  

    fig.legend(prop={'size':12})
    plt.show()

"""
argv[2] := file that contains plot data.
"""
def plotTraffic(filename, title):
    fig, ax = plt.subplots(1, 1)

    fig.set_facecolor('#ebeef2')
    ax.set_facecolor('#e6e9ed')

    fig.suptitle(title, fontsize=20)

    with open(filename) as file:
        data          = csv.reader(file, delimiter=',')
        discrepancies = list(map(float, next(data)))[13:]
        occurrences   = list(map(lambda x: 100*float(x), next(data)))[13:]
        discrepancies = list(map(float, next(data)))[13:]
        occurrences   = list(map(lambda x: 100*float(x), next(data)))[13:]
        
        ax.bar(discrepancies, height=occurrences, color='gray')
        ax.set_xlabel('Bits of Precision Over Float32', fontsize=16)
        ax.set_ylabel('Percent of Cases'       , fontsize=16)
        ax.xaxis.set_ticks(discrepancies)#discrepancies[0:len(discrepancies):2])

    plt.show()

def barPlot(filename, ylabel):
    fig, ax = plt.subplots(1, 1)

    fig.set_facecolor('#ebeef2')
    ax .set_facecolor('#e6e9ed')
    ax .set_yscale('log')

    #fig.suptitle(title, fontsize=20)

    with open(filename) as file:
        data      = csv.reader(file, delimiter=',')
        types     = list(next(data))
        residuals = list(map(float, next(data)))
        
        ax.bar(range(len(residuals)), height=residuals, color='gray', tick_label=types)
        ax.set_xlabel('Numeric Type', fontsize=16)
        ax.set_ylabel(ylabel        , fontsize=16)
        ax.tick_params(axis='x',rotation=30)
    plt.tight_layout()
    plt.show()

def plotAdvantage(filename, matrixname=""): 

    fig, ax = plt.subplots(1, 1)
    ax.set_facecolor('#e6e9ed')

    fig.set_facecolor('#ebeef2')
    
    with open(filename, 'r') as file:
        for plot in csv.reader(file, delimiter=','):
            color         = 'blue'
            bit_advantages = list(map(float, plot))
            ax.plot(range(len(plot)), bit_advantages, color=color, marker=None)

        if matrixname != "": title = 'Posit Bit Advantage for ' + matrixname
        else               : title = 'Posit Bit Advantage'

        ax.set_title (title, fontsize=16)
        ax.set_xlabel('Iteration' , fontsize=14)
        ax.set_ylabel('Advantage' , fontsize=14)  
    plt.show()


"""
first argument of 0 -> CG plot.
first argument of 1 -> traffic plot.
first argument of 2 -> advantage plot. 

First argument identifies the type of plot. for second argument and above, refer to comments above. 
"""
if __name__ == "__main__":
    
    if sys.argv[1] == '0':
        plot(sys.argv[2], sys.argv[3], *[arg for arg in sys.argv[4:] if '=' not in arg], \
            **dict(arg.split('=') for arg in sys.argv[4:] if '=' in arg)) 
    elif sys.argv[1] == '1':
        plotTraffic(sys.argv[2], 'Posit Advantage for Matrix Market (ES=3)')
    elif sys.argv[1] == '2':
        plotAdvantage(sys.argv[2], sys.argv[3] if len(sys.argv) > 3 else "")
    plotPrecision({'es2precisionRelative.txt':'Posit(32, 2)', 'positprecisionRelative.txt':'Posit(32, 3)', \
        'floatprecisionRelative.txt':'Float32'})
    

    

