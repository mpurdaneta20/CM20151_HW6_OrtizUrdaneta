all: 4body.x orbitas-1yrs.txt orbitas-1000yrs.txt orbitas-1yrs.png orbitas-1000yrs.png 
 
orbitas-1000yrs.png: Trayectoria.py orbitas-1000yrs.txt
	
orbitas-1yrs.png: Trayectoria.py orbitas-1yrs.txt
	python Trayectoria.py
orbitas-1000yrs.txt: 4body.x
	./4body.x ic.txt 0.01 1000
orbitas-1yrs.txt: 4body.x
	./4body.x ic.txt 0.001 1
4body.x: 4body.c
	cc 4body.c -o 4body.x -lm
clear:
	rm -f 4body.x orbitas-1yrs.txt orbitas-1000yrs.txt orbitas-1yrs.png orbitas-1000yrs.png 
