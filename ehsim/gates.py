import numpy as np

zero_ket = np.array([1, 0], dtype=np.complex_)
one_ket = np.array([0, 1], dtype=np.complex_)
plus_ket = np.array([1 / np.sqrt(2), 1 / np.sqrt(2)], dtype=np.complex_)
minus_ket = np.array([1 / np.sqrt(2), -1 / np.sqrt(2)], dtype=np.complex_)

X = {
    'gate':'X',
    'rules': [{
        'match':[zero_ket],
        'replace':[one_ket]
    },
    {
        'match':[one_ket],
        'replace':[zero_ket]
    }]
}

H = {
    'gate':'H',
    'rules': [{
        'match':[zero_ket],
        'replace':[plus_ket]
    },
    {
        'match':[one_ket],
        'replace':[minus_ket]
    }]
}

SWAP = {
    'gate':'SWAP',
    'rules': [{
        'match':[zero_ket, one_ket],
        'replace':[one_ket, zero_ket]
    },
    {
        'match':[one_ket, zero_ket],
        'replace':[zero_ket, one_ket]
    }]
} 

CX = {
    'gate':'CX',
    'rules': [{
        'match':[one_ket, one_ket],
        'replace':[one_ket, zero_ket]
    },
    {
        'match':[one_ket, zero_ket],
        'replace':[one_ket, one_ket]
    }]
}

CCX = {
    'gate':'CCX',
    'rules': [{
        'match':[one_ket, one_ket, one_ket],
        'replace':[one_ket, zero_ket]
    },
    {
        'match':[one_ket, one_ket, zero_ket],
        'replace':[one_ket, one_ket]
    }]
}