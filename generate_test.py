import numpy as np

scaffolds = [3,8,4,2,10,5,5,2]
#scaffolds = [1,2,3]
print("Cumsum of test scaffolds:")
print(np.cumsum(scaffolds))
groups = [(0,3),(2,6)]
#groups = []

def gen_test():
    matrix = np.ones(shape=(sum(scaffolds),sum(scaffolds)))

    intra_const = 100

    for ind,val in enumerate(scaffolds):
        st = np.cumsum(scaffolds)[ind] - val
        end = np.cumsum(scaffolds)[ind]

        matrix[st:end,st:end] += intra_const

    for group in groups:
        for ind1,val1 in enumerate(group):
            for ind2,val2 in enumerate(group):
                st1 = np.cumsum(scaffolds)[val1] - scaffolds[val1]
                end1 = np.cumsum(scaffolds)[val1]
                st2 = np.cumsum(scaffolds)[val2] - scaffolds[val2]
                end2 = np.cumsum(scaffolds)[val2]
                print (st1,end1,st2,end2)
                matrix[st1:end1,st2:end2] = intra_const

    borders = np.cumsum(scaffolds)
    return  matrix,scaffolds,[0]*len(scaffolds),[str(i) for i in range(len(scaffolds))]
