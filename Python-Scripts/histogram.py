#!/usr/bin/env python
import matplotlib.pyplot as plt
import numpy as np

with open("../3_kmer", "r") as ins:
    array = []
    
    for line in ins:
        word = []
        word = line.split()
        array.append(int(word[1]))
plt.hist(array, bins=500)
plt.yscale('log', nonposy='clip')
plt.title("Number of reads associated with kmer [K=3]")
plt.show()

