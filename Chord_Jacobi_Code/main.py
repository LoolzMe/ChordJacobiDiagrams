import chords

m = 3
fourT = chords.gen_4T(m)
diagrams = chords.circ_diag(m)
print(len(diagrams))

for diag in diagrams:
    print(diag)

# for T in fourT:
#     print(T)

diagrams = list(map(chords.Ch_diag, diagrams))

LS = chords.ChordsLinearSpace(diagrams)
print(LS.getMatrix())

LS.reduceByRelation(fourT)

print(LS.getBasis())
print(LS.getMatrix())
print(LS.getCorrDict())


diagrams_new = []
for diag in LS.getBasis():
    diagrams_new.append(diag.seq)

print(diagrams_new)