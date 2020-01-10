import numpy as np
import matplotlib.pyplot as plt
import sys
import math
from PIL import Image

xArr = []
yArr = []
rArr = []

inputFile = open("/Users/richardzane/Documents/EAGLE/ulps/HeatMapThubber/files/output.txt","r")
for l in inputFile:
    print(l)
    lS = l.split()
    if lS[0] == "Area:":
        bY = lS[1]
        bX = lS[2]
        print(bY[1:])
        print(bX)
        if bY[0] == '-':
            print("hello")
            bY = float(bY[1:])*(-1.0)
            bX = float(bX)
        elif bX[0] == '-':
            bX = float(bX[1:])*(-1.0)
            bY = float(bY)
        tY = float(lS[3])
        tX = float(lS[4])
    else:
        yArr.append(float(lS[0]) - bY)
        xArr.append(float(lS[1]) - bX)
        tR = math.sqrt(float(lS[2])/3.14)
        rArr.append(tR)

print(sys.argv)

w = tX - bX
h = tY - bY

for i in range(len(xArr)):
    xArr[i] = w - xArr[i]

print(w)
print(h)
# intervals in x-, y- directions, mm
dx = dy = 0.1

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


x = [1,2,3,4,5]
print(x[1:-1])

for rad in rArr:
    r2Arr.append(rad**2)

for i in range(nx):
    for j in range(ny):
        for k in range(len(rArr)):
            thisP = (i*dx-xArr[k])**2 + (j*dy-yArr[k])**2
            pArr.append(thisP)
            if thisP < r2Arr[k]:
                u0[i,j] = Thot
                
sDivideX = int(w/2)
sDivideY = int(h/2)

def do_timestep(u0, u):

        
#    if u0[m,m] <= sDivideX:
#        u0, u = do_timestep(u0, u, "thubber")
#    else:
#        u0, u = do_timestep(u0, u, "pcb")
    
    
   
    # Propagate with forward-difference in time, central-difference in space
    D = Dpcb
    dt = dtPcb
    u[1:sDivideX-1, 1:sDivideY-1] = u0[1:sDivideX-1, 1:sDivideY-1] + D * dt * ((u0[2:sDivideX, 1:sDivideY-1] - 2*u0[1:sDivideX-1, 1:sDivideY-1] +u0[:sDivideX-2, 1:sDivideY-1])/dx2 + (u0[1:sDivideX-1, 2:sDivideY] - 2*u0[1:sDivideX-1, 1:sDivideY-1] + u0[1:sDivideX-1, :sDivideY-2])/dy2)
    D = Dthub
    dt = dtThub
    u[sDivideX+1:-1, sDivideY+1:-1] = u0[sDivideX+1:-1, sDivideY+1:-1] + D * dt * ((u0[sDivideX+2:, sDivideY+1:-1] - 2*u0[sDivideX+1:-1, sDivideY+1:-1] +u0[sDivideX:-2, sDivideY+1:-1])/dx2 + (u0[sDivideX+1:-1, sDivideY+2:] - 2*u0[sDivideX+1:-1, sDivideY+1:-1] + u0[sDivideX+1:-1, sDivideY:-2])/dy2)
    
#    u[1:-1, 1:-1] = u0[1:-1, 1:-1] + D * dt * ((u0[2:, 1:-1] - 2*u0[1:-1, 1:-1] +u0[:-2, 1:-1])/dx2 + (u0[1:-1, 2:] - 2*u0[1:-1, 1:-1] + u0[1:-1, :-2])/dy2)
    u0 = u.copy()
    return u0, u
    


# Number of timesteps
nsteps = 10001
#timeLim = 2000
mfig = [4000,10000]
fignum = 0
time = 0
fig = plt.figure()
for m in range(nsteps):
    u0,u = do_timestep(u0,u)
    if m in mfig:
        
        time = m*dtThub*1000
        fignum += 1
        print(m, fignum, time)
        im = plt.imshow(u.copy(), cmap=plt.get_cmap('hot'), vmin=Tcool,vmax=Thot)
        plt.axis('off')
#        ax.set_axis_off()
        # ax.set_title('{:.1f} ms' .format(time))
#        if time >= timeLim:
#            break

print("time:")
print(str(time) + " ms")
fig.subplots_adjust(right=1,left=0,top=1,bottom=0,hspace=0,wspace=0)
plt.margins(0,0)
# cbar_ax = fig.add_axes([0.9, 0.15, 0.03, 0.7])
# cbar_ax.set_xlabel('$T$ / K', labelpad=20)
#fig.colorbar(im, cax=cbar_ax)
#fig.set_size_inches(10,10)
axe = plt.Axes(fig, [0., 0., 1., 1.])
axe.set_axis_off()

fig.savefig("/Users/richardzane/Documents/EAGLE/ulps/HeatMapThubber/files/heatmap.png",bbox_inches='tight',pad_inches=0)
Image.open("/Users/richardzane/Documents/EAGLE/ulps/HeatMapThubber/files/heatmap.png").save("/Users/richardzane/Documents/EAGLE/ulps/HeatMapThubber/files/heatmap.bmp")

img = Image.open("/Users/richardzane/Documents/EAGLE/ulps/HeatMapThubber/files/heatmap.bmp")
newimg = img.convert(mode='P', colors=256)
newimg.save("/Users/richardzane/Documents/EAGLE/ulps/HeatMapThubber/files/heatmap.bmp")

plt.show()
