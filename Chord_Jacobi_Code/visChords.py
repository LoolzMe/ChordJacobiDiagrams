import cv2
import numpy as np 
import chords

def chunks(lst, n):
    """Yield successive n-sized chunks from lst."""
    for i in range(0, len(lst), n):
        yield lst[i:i + n]


width = 640
height = 480


window = np.zeros((height, width, 3)) + 255

m = 4
diagrams = chords.circ_diag(m)
fourT = chords.gen_4T(m)

diagrams = list(map(chords.Ch_diag, diagrams))

LS = chords.ChordsLinearSpace(diagrams)

LS.reduceByRelation(fourT)

diagrams = LS.getBasis()
diagrams_new = []

for diag in diagrams:
    diagrams_new.append(diag.seq)

diagrams = diagrams_new

n_diagrams = len(diagrams)

if len(diagrams) > 24:
    diagrams = chunks(diagrams, 24)
else:
    diagrams = [diagrams]

R = 30
space = 20
current_pos = [0, 0]
diagram_color = (255, 0, 0)
line_color = (0, 0, 255)
diagram_thickness = 2

def drawDiagram(diag):
    global window, current_pos
    reg_pol = []
    n = len(diag)

    window = cv2.circle(window, center=tuple(current_pos), radius=R, color=diagram_color, thickness=diagram_thickness)
    
    for phi in range(n):
        reg_pol.append((int(current_pos[0] + R * np.cos(2 * np.pi * phi / n)), int(current_pos[1] + R * np.sin(2 * np.pi * phi / n))))
    
    bg_chs = {}

    for index, ch in enumerate(diag):
        if ch in bg_chs:
            window = cv2.line(window, pt1=reg_pol[bg_chs[ch]], pt2=reg_pol[index], color=line_color, thickness=diagram_thickness)
            bg_chs.pop(ch)
        else:
            bg_chs[ch] = index


current_pos = [130, 50]


for index, obj in enumerate(diagrams):
    current_pos[0] = 130
    current_pos[1] = 50
    if index != 0:
        if cv2.waitKey(0):
            window = np.zeros((height, width, 3)) + 255
            current_pos = [130, 50]

    window = cv2.putText(window, "{}".format(n_diagrams), (20, 35), cv2.FONT_HERSHEY_TRIPLEX, fontScale=1, color=diagram_color, thickness=diagram_thickness//2)
    for index, diag in enumerate(obj):
        drawDiagram(diag)
        current_pos[0] += 9*R // 2
        if index % 4 == 3:
            current_pos[0] = 130
            current_pos[1] += 2*R + 10


    cv2.imshow("Diagrams", window)

cv2.waitKey(0)