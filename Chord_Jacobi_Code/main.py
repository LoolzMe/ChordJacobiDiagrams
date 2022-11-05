import chords


diags = chords.circ_diag(3, True)

fourT = chords.gen_4T(4)

while True:
    try:
        obj = next(fourT)
    except StopIteration:
        break
    
    print(obj)
    