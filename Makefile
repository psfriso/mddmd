EXEDIR =

all: discrete
	
discrete:
	cd discrete; make

install: all
	cp discrete/discrete discrete

clean:
	cd discrete; make clean
	rm discrete$
        	
