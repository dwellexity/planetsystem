from visual import *
scene.fullscreen=True
# axisx = arrow(pos=(0,0,0), axis=(2,0,0), shaftwidth=0.001, color=color.red)
# axisy = arrow(pos=(0,0,0), axis=(0,2,0), shaftwidth=0.001, color=color.green)
# axisz = arrow(pos=(0,0,0), axis=(0,0,2), shaftwidth=0.001, color=color.blue)
# label(pos=(10,1,0),text="x-axis")
# label(pos=(0,11,1),text="y-axis")
# label(pos=(1,0,11),text="z-axis")
s=[]
r=[]
#raise SystemExit
params=open('params.dat','r')
n=int(params.readline())
rad=params.readline()
rad=rad.split()
for i in range(0,n):
	r.append(float(rad[i]))
print n,r
num_lines = sum(1 for line in open('xyz.dat'))
loop=num_lines/n
f=open('xyz.dat','r')
for i in range(0,n):
	line=f.readline()
	row = line.split()
	x = float(row[0])
	y = float(row[1])
	z = float(row[2])
	s.append(sphere(pos=[x,y,z],radius=r[i],color=color.cyan,make_trail=True))
for i in range(1,loop):
	rate(10000)
	for j in range(0,n):
		line=f.readline()
		row = line.split()
		x = float(row[0])
		y = float(row[1])
		z = float(row[2])
		s[j].pos = [x,y,z]
