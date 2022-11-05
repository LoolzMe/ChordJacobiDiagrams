import cv2
import numpy as np 
import chords


# To navigate through the 4T relation equations press any button on your keyboard
# On English layout by the way


width = 640
height = 480

window = np.zeros((height, width, 3)) + 255

fourT = chords.gen_4T(4)

R = 30
space = 20
current_pos = [0, 0]
diagram_color = (255, 0, 0)
line_color = (0, 0, 255)
diagram_thickness = 2

def drawDiagram(coef, diag):
    global window
    reg_pol = []
    n = len(diag)

    window = cv2.circle(window, center=tuple(current_pos), radius=R, color=diagram_color, thickness=diagram_thickness)
    if coef > 0: 
        coef = '+' + str(coef)
    else:
        coef = str(coef)
    window = cv2.putText(window, coef, (current_pos[0] - 5*R//2, current_pos[1] + 10), cv2.FONT_HERSHEY_TRIPLEX, fontScale=1, color=diagram_color, thickness=diagram_thickness//2)

    for phi in range(n):
        reg_pol.append((int(current_pos[0] + R * np.cos(2 * np.pi * phi / n)), int(current_pos[1] + R * np.sin(2 * np.pi * phi / n))))
    
    bg_chs = {}

    for index, ch in enumerate(diag):
        if ch in bg_chs:
            window = cv2.line(window, pt1=reg_pol[bg_chs[ch]], pt2=reg_pol[index], color=line_color, thickness=diagram_thickness)
            bg_chs.pop(ch)
        else:
            bg_chs[ch] = index


current_pos = [100, 30]

for obj in fourT:
    current_pos[0] = 100
    if current_pos[1] + 2*R >= height:
        if cv2.waitKey(0) & 0xFF == ord('d'):
            window = np.zeros((height, width, 3)) + 255
            current_pos = [100, 30]

    for coef, diag in obj:
        drawDiagram(coef, diag)
        current_pos[0] += 9*R // 2
    
    current_pos[1] += 2*R + 10

    cv2.imshow("Diagrams", window)

#drawDiagram(+1, [1, 1, 2, 2, 3, 4, 3, 4])




#cv2.waitKey(0) 