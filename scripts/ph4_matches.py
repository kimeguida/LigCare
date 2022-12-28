ph4_matching = {
    
    # (ligand, cavity), (ligand, cavity)
    # ligand A, B
    # cavity: A', B'

    # (A, A') , (B, B')

    #A, B = CA, CA
    (('CA', 'CA'), ('CA', 'CA')): 0.5,

    #A, B = CA, CZ
    (('CA', 'CA'), ('CZ', 'CZ')): 1,

    #A, B = CA, O
    (('CA', 'CA'), ('O', 'O')): 1,

    #A, B = CA, OD1
    (('CA', 'CA'), ('OD1', 'OD1')): 1,

    #A, B = CA, OG
    (('CA', 'CA'), ('OG', 'OG')): 1,

    #A, B = CA, N
    (('CA', 'CA'), ('N', 'N')): 1,

    #A, B = CA, NZ
    (('CA', 'CA'), ('NZ', 'NZ')): 1,



    #A, B = CZ, CA
    (('CZ', 'CZ'), ('CA', 'CA')): 1,

    #A, B = CZ, CZ
    (('CZ', 'CZ'), ('CZ', 'CZ')): 1,

    #A, B = CZ, O
    (('CZ', 'CZ'), ('O', 'O')): 1,

    #A, B = CZ, OD1
    (('CZ', 'CZ'), ('OD1', 'OD1')): 1,

    #A, B = CZ, OG
    (('CZ', 'CZ'), ('OG', 'OG')): 1,

    #A, B = CZ, N
    (('CZ', 'CZ'), ('N', 'N')): 1,

    #A, B = CZ, NZ
    (('CZ', 'CZ'), ('NZ', 'NZ')): 1,


    
    #A, B = O, CA
    (('O', 'O'), ('CA', 'CA')): 1,

    #A, B = O, CZ
    (('O', 'O'), ('CZ', 'CZ')): 1,

    #A, B = O, O
    (('O', 'O'), ('O', 'O')): 1,

    #A, B = O, OD1
    (('O', 'O'), ('OD1', 'OD1')): 1,

    #A, B = O, OG
    (('O', 'O'), ('OG', 'OG')): 1,

    #A, B = O, N
    (('O', 'O'), ('N', 'N')): 1,

    #A, B = O, NZ
    (('O', 'O'), ('NZ', 'NZ')): 1,


    
    #A, B = OD1, CA
    (('OD1', 'OD1'), ('CA', 'CA')): 1,

    #A, B = OD1, CZ
    (('OD1', 'OD1'), ('CZ', 'CZ')): 1,

    #A, B = OD1, O
    (('OD1', 'OD1'), ('O', 'O')): 1,

    #A, B = OD1, OD1
    (('OD1', 'OD1'), ('OD1', 'OD1')): 1,

    #A, B = OD1, OG
    (('OD1', 'OD1'), ('OG', 'OG')): 1,

    #A, B = OD1, N
    (('OD1', 'OD1'), ('N', 'N')): 1,

    #A, B = OD1, NZ
    (('OD1', 'OD1'), ('NZ', 'NZ')): 1,


    
    #A, B = OG, CA
    (('OG', 'OG'), ('CA', 'CA')): 1,

    #A, B = OG, CZ
    (('OG', 'OG'), ('CZ', 'CZ')): 1,

    #A, B = OG, O
    (('OG', 'OG'), ('O', 'O')): 1,

    #A, B = OG, OD1
    (('OG', 'OG'), ('OD1', 'OD1')): 1,

    #A, B = OG, OG
    (('OG', 'OG'), ('OG', 'OG')): 1,

    #A, B = OG, N
    (('OG', 'OG'), ('N', 'N')): 1,

    #A, B = OG, NZ
    (('OG', 'OG'), ('NZ', 'NZ')): 1,



    #A, B = N, CA
    (('N', 'N'), ('CA', 'CA')): 1,

    #A, B = N, CZ
    (('N', 'N'), ('CZ', 'CZ')): 1,

    #A, B = N, O
    (('N', 'N'), ('O', 'O')): 1,

    #A, B = N, OD1
    (('N', 'N'), ('OD1', 'OD1')): 1,

    #A, B = N, OG
    (('N', 'N'), ('OG', 'OG')): 1,

    #A, B = N, N
    (('N', 'N'), ('N', 'N')): 1,

    #A, B = N, NZ
    (('N', 'N'), ('NZ', 'NZ')): 1,



    #A, B = NZ, CA
    (('NZ', 'NZ'), ('CA', 'CA')): 1,

    #A, B = NZ, CZ
    (('NZ', 'NZ'), ('CZ', 'CZ')): 1,

    #A, B = NZ, O
    (('NZ', 'NZ'), ('O', 'O')): 1,

    #A, B = NZ, OD1
    (('NZ', 'NZ'), ('OD1', 'OD1')): 1,

    #A, B = NZ, OG
    (('NZ', 'NZ'), ('OG', 'OG')): 1,

    #A, B = NZ, N
    (('NZ', 'NZ'), ('N', 'N')): 1,

    #A, B = NZ, NZ
    (('NZ', 'NZ'), ('NZ', 'NZ')): 1,



}


"""
    #A, B = DU, CA
    (('DU', 'DU'), ('CA', 'CA')): 0.5,
    (('DU', 'CA'), ('CA', 'CA')): 0.5,
    (('DU', 'DU'), ('CA', 'DU')): 0.5,
    (('DU', 'CA'), ('CA', 'DU')): 0.5,
    #A, B = DU, CZ
    (('DU', 'DU'), ('CZ', 'CZ')): 1,
    (('DU', 'CA'), ('CZ', 'CZ')): 1,
    (('DU', 'DU'), ('CZ', 'CA')): 0.5,
    (('DU', 'CA'), ('CZ', 'CA')): 0.5,
    #A, B = DU, O
    (('DU', 'DU'), ('O', 'O')): 1,
    (('DU', 'CA'), ('O', 'O')): 1,
    (('DU', 'DU'), ('O', 'OD1')): 0.75,
    (('DU', 'CA'), ('O', 'OD1')): 0.75,
    (('DU', 'DU'), ('O', 'OG')): 0.75,
    (('DU', 'CA'), ('O', 'OG')): 0.75,
    #A, B = DU, OD1
    (('DU', 'DU'), ('OD1', 'OD1')): 1,
    (('DU', 'CA'), ('OD1', 'OD1')): 1,
    (('DU', 'DU'), ('OD1', 'O')): 0.75,
    (('DU', 'CA'), ('OD1', 'O')): 0.75,
    (('DU', 'DU'), ('OD1', 'OG')): 0.75,
    (('DU', 'CA'), ('OD1', 'OG')): 0.75,
    #A, B = DU, OG
    (('DU', 'DU'), ('OG', 'OG')): 1,
    (('DU', 'CA'), ('OG', 'OG')): 1,
    (('DU', 'DU'), ('OG', 'O')): 0.75,
    (('DU', 'CA'), ('OG', 'O')): 0.75,
    (('DU', 'DU'), ('OG', 'OD1')): 0.75,
    (('DU', 'CA'), ('OG', 'OD1')): 0.75,
    #A, B = DU, N
    (('DU', 'DU'), ('N', 'N')): 1,
    (('DU', 'CA'), ('N', 'N')): 1,
    (('DU', 'DU'), ('N', 'NZ')): 0.75,
    (('DU', 'CA'), ('N', 'NZ')): 0.75,
    #A, B = DU, NZ
    (('DU', 'DU'), ('NZ', 'NZ')): 1,
    (('DU', 'CA'), ('NZ', 'NZ')): 1,
    (('DU', 'DU'), ('NZ', 'N')): 0.75,
    (('DU', 'CA'), ('NZ', 'N')): 0.75,
    #A, B = DU, DU
    (('DU', 'DU'), ('DU', 'DU')): 0.5,
    (('DU', 'CA'), ('DU', 'DU')): 0.5,
    (('DU', 'DU'), ('DU', 'CA')): 0.5,
    (('DU', 'CA'), ('DU', 'CA')): 0.5,
    """

    

ph4_pairs = {

    ('CA', 'CA'),
    ('CA', 'CZ'),
    ('CA', 'O'),
    ('CA', 'OD1'),
    ('CA', 'OG'),
    ('CA', 'N'),
    ('CA', 'NZ'),
    ('CA', 'DU'),

    ('CZ', 'CA'),
    ('CZ', 'CZ'),
    ('CZ', 'O'),
    ('CZ', 'OD1'),
    ('CZ', 'OG'),
    ('CZ', 'N'),
    ('CZ', 'NZ'),
    ('CZ', 'DU'),

    ('O', 'CA'),
    ('O', 'CZ'),
    ('O', 'O'),
    ('O', 'OD1'),
    ('O', 'OG'),
    ('O', 'N'),
    ('O', 'NZ'),
    ('O', 'DU'),

    ('OD1', 'CA'),
    ('OD1', 'CZ'),
    ('OD1', 'O'),
    ('OD1', 'OD1'),
    ('OD1', 'OG'),
    ('OD1', 'N'),
    ('OD1', 'NZ'),
    ('OD1', 'DU'),

    ('OG', 'CA'),
    ('OG', 'CZ'),
    ('OG', 'O'),
    ('OG', 'OD1'),
    ('OG', 'OG'),
    ('OG', 'N'),
    ('OG', 'NZ'),
    ('OG', 'DU'),

    ('N', 'CA'),
    ('N', 'CZ'),
    ('N', 'O'),
    ('N', 'OD1'),
    ('N', 'OG'),
    ('N', 'N'),
    ('N', 'NZ'),
    ('N', 'DU'),


    ('NZ', 'CA'),
    ('NZ', 'CZ'),
    ('NZ', 'O'),
    ('NZ', 'OD1'),
    ('NZ', 'OG'),
    ('NZ', 'N'),
    ('NZ', 'NZ'),
    ('NZ', 'DU'),

    ('DU', 'CA'),
    ('DU', 'CZ'),
    ('DU', 'O'),
    ('DU', 'OD1'),
    ('DU', 'OG'),
    ('DU', 'N'),
    ('DU', 'NZ'),
    ('DU', 'DU'),
}

frequencies = {
    "CA": 0.3817,
    "CZ": 0.1110,
    "O": 0.0659,
    "OG": 0.0867,
    "OD1": 0.0593,
    "N": 0.1376,
    "NZ": 0.0579,
    "DU": 0.0999
}

frequencies = {
    "CA": 1,
    "CZ": 1,
    "O": 1,
    "OG": 1,
    "OD1": 1,
    "N": 1,
    "NZ": 1,
    "DU": 1
}

frequencies = {
    "CA": 1,
    "CZ": 0.33,
    "O": 0.33,
    "OG": 0.33,
    "OD1": 0.33,
    "N": 0.33,
    "NZ": 0.33,
    "DU": 1
}


ph4_scores = {k: 1/(frequencies[k[0]]*frequencies[k[1]]) for 
            k in ph4_pairs
}


if __name__ == '__main__':
    print(ph4_scores)

    print('CA, CA', ph4_scores[('CA', 'CA')])
    print('CZ, CZ', ph4_scores[('CZ', 'CZ')])
    print('CA, CZ', ph4_scores[('CA', 'CZ')])

    print('N, NZ', ph4_scores[('N', 'NZ')])
    print('N, N', ph4_scores[('N', 'N')])
    print('NZ, NZ', ph4_scores[('NZ', 'NZ')])

    print('O, O', ph4_scores[('O', 'O')])
    print('OD1, OD1', ph4_scores[('OD1', 'OD1')])
    print('OG, OG', ph4_scores[('OG', 'OG')])
    print('O, OG', ph4_scores[('O', 'OG')])
    print('O, OD1', ph4_scores[('O', 'OD1')])
    print('OD1, OG', ph4_scores[('OD1', 'OG')])
