default: ../ IFOS2D

.PHONY: aff
aff: 
	$(MAKE) -C ../contrib/aff install

.PHONY: fourier
fourier:
	$(MAKE) -C ../contrib/fourier install

.PHONY: stfinv
stfinv:
	$(MAKE) -C ../contrib/stfinv install

.PHONY: libcseife
libcseife:
	$(MAKE) -C ../contrib/libcseife install


.PHONY: IFOS2D 
IFOS2D: aff fourier stfinv libcseife
ifdef MODEL
	$(MAKE) -C ../src/ IFOS2D MODEL=$(MODEL)
else
	$(MAKE) -C ../src/ IFOS2D 
endif


.PHONY: clean
clean: 
	rm -rf ../bin/IFOS2D
	$(MAKE) -C ../contrib/aff clean
	$(MAKE) -C ../contrib/fourier clean
	$(MAKE) -C ../contrib/stfinv clean
	$(MAKE) -C ../contrib/libcseife clean
	$(MAKE) -C ../src/ clean
