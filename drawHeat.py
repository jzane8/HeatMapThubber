import numpy as np
import matplotlib.pyplot as plt
import sys
import math
# from wand.image import Image
from PIL import Image
# MAGICK_HOME = "C:/Program Files/ImageMagick-6.9.10-Q16"

xArr = []
yArr = []
rArr = []

inputFile = open("/Users/richardzane/Documents/EAGLE/ulps/HeatULP/output.txt","r")
for l in inputFile:
    print(l)
    lS = l.split()
    if lS[0] == "Area:":
        bY = lS[1]
        bX = lS[2]
        if bY[0] == '-':
            print("hello")
            bY = float(bY[1:])*(-1.0)
        if bX[0] == '-':
            bX = float(bX[1:])*(-1.0)

        tY = float(lS[3])
        tX = float(lS[4])
    else:
        yArr.append(float(lS[0]) - bY)
        xArr.append(float(lS[1]) - bX)
        tR = math.sqrt(float(lS[2])/3.14)
        rArr.append(tR)

print(sys.argv)

# plate size, mm
# w = h = 20.
w = tX - bX
h = tY - bY

for i in range(len(xArr)):
    xArr[i] = w - xArr[i]

print(w)
print(h)
# intervals in x-, y- directions, mm
dx = dy = 0.1
# Thermal diffusivity of steel, mm2.s-1
# D = 4.
# Thermal diffusivity of epoxy resin, aka pcb material
Dpcb = 0.13
#thermal diffusivity of thubber
Dthub = 8.65

Tcool, Thot = 300, 350

nx, ny = int(w/dx), int(h/dy)

dx2, dy2 = dx*dx, dy*dy
# dt = dx2 * dy2 / (2 * D * (dx2 + dy2))
dtPcb = dx2 * dy2 / (2 * Dpcb * (dx2 + dy2))
dtThub = dx2 * dy2 / (2 * Dthub * (dx2 + dy2))

u0 = Tcool * np.ones((nx, ny))
u = np.empty((nx, ny))

# Initial conditions - ring of inner radius r, width dr centred at (cx,cy) (mm)

r2Arr = []
pArr = []

for rad in rArr:
    r2Arr.append(rad**2)

for i in range(nx):
    for j in range(ny):
        for k in range(len(rArr)):
            thisP = (i*dx-xArr[k])**2 + (j*dy-yArr[k])**2
            pArr.append(thisP)
            if thisP < r2Arr[k]:
                u0[i,j] = Thot

def do_timestep(u0, u, mat):
    if mat == "pcb":
        D = Dpcb
        dt = dtPcb
    elif mat == "thubber":
        D = Dthub
        dt = dtThub
    # Propagate with forward-difference in time, central-difference in space
    u[1:-1, 1:-1] = u0[1:-1, 1:-1] + D * dt * (
          (u0[2:, 1:-1] - 2*u0[1:-1, 1:-1] + u0[:-2, 1:-1])/dx2
          + (u0[1:-1, 2:] - 2*u0[1:-1, 1:-1] + u0[1:-1, :-2])/dy2 )
    u0 = u.copy()
    return u0, u

# Number of timesteps
nsteps = 501
# Output 4 figures at these timesteps
timeLim = 2000
mfig = [500]
fignum = 0
fig = plt.figure()
for m in range(nsteps):
    u0, u = do_timestep(u0, u, "pcb")
    if m in mfig:
        time = m*dtPcb*1000
        fignum += 1
        print(m, fignum, time)
        ax = fig.add_subplot()
        im = ax.imshow(u.copy(), cmap=plt.get_cmap('hot'), vmin=Tcool,vmax=Thot)
        ax.set_axis_off()
        # ax.set_title('{:.1f} ms' .format(time))
        if time >= timeLim:
            break
fig.subplots_adjust(right=1,left=0,top=1,bottom=0)
# cbar_ax = fig.add_axes([0.9, 0.15, 0.03, 0.7])
# cbar_ax.set_xlabel('$T$ / K', labelpad=20)
# fig.colorbar(im, cax=cbar_ax)
fig.savefig("/Users/richardzane/Documents/EAGLE/ulps/HeatULP/heatmap.png")
Image.open("/Users/richardzane/Documents/EAGLE/ulps/HeatULP/heatmap.png").save("/Users/richardzane/Documents/EAGLE/ulps/HeatULP/heatmap.bmp")
img = Image.open("/Users/richardzane/Documents/EAGLE/ulps/HeatULP/heatmap.bmp")
newimg = img.convert(mode='P', colors=16)
newimg.save("/Users/richardzane/Documents/EAGLE/ulps/HeatULP/heatmap.bmp")
# with Image(filename="heatimg.png") as img:
#     with img.convert('bmp') as converted:
#         newImg = converted.quantize(256,'rgb',0,False,False)
#         converted.save(filename='heatbmp.bmp')

plt.show()
