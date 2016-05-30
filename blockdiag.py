import numpy as np

class Queue:
    def __init__(self):
        self.items = []
    def isEmpty(self):
        return self.items == []
    def enqueue(self, item):
        self.items.insert(0,item)
    def dequeue(self):
        return self.items.pop()
    def size(self):
        return len(self.items)


def blockdiag2blocks(matr):
    block_cnt = 0
    elem_cnt=matr.shape[0]
    print matr.shape, elem_cnt
    elem_checked = [False for i in xrange(elem_cnt)]
    
    elem_connections = [[] for i in xrange(elem_cnt)]
    block_first_elem = []
    block_elements = []

    for i in xrange(elem_cnt):
        for j in xrange(elem_cnt):
            if matr[i,j]!=0:
                elem_connections[i].append(j)
    queue = Queue()
    for i in xrange(elem_cnt):
        if elem_checked[i]:
            continue
        else:
            queue.enqueue(i)
            elem_checked[i] = True
            block_cnt+=1
            block_first_elem.append(i)
            block_elements.append([i])
            while not queue.isEmpty():
                current = queue.dequeue()
                for j in elem_connections[current]:
                    if not elem_checked[j]:
                        queue.enqueue(j)
                        block_elements[-1].append(j)
                        elem_checked[j] = True
    # make blocks
    matr_out = []
    for subgraph in block_elements:
        w = np.zeros((len(subraph),len(subgraph)))
        for i in xrange(len(subraph)):
            for j in xrange(len(subraph)):
                w[i,j] = matr[subraph[i],subraph[j]]
        matr_out.append(w)
    print elem_cnt, block_cnt, block_first_elem
    #print block_elements
    print [len(m_i) for m_i in block_elements]
    return matr_out, block_elements
    
    
    
    