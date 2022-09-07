import numpy as np
import matplotlib.pyplot as plt

def read_data(file:str):
    with open(file, "r") as f:
        lines = f.readlines()
        lines = [l.rstrip("\n").lstrip().split() for l in lines]
        lines = [[float(r) for r in l] for l in lines]
        lines = np.array(lines)
    return lines

if __name__ == "__main__":
    read_file = "c.out" 
    data = read_data(read_file)

    c = np.zeros((100,100))
    for d in data:
        index1, index2 = int(d[0]), int(d[1])
        c[index1,index2] = d[-1]

    dataP = read_data("hit.out")
    p = []
    for d in dataP:
        if d[-1] == 1:
            p.append([d[0],d[1]])
    p = np.array(p)


    fig = plt.figure(dpi=300)
    ax  = fig.add_subplot(111)
    levels = np.linspace(0,0.5,100)
    c = ax.contourf(c.T, cmap='jet', levels=levels)
    ax.scatter(p[:,0], p[:,1],c='r')
    fig.colorbar(c)
    plt.savefig("test.png")


