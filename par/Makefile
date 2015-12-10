default: ../ denise

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


.PHONY: denise 
denise: aff fourier stfinv libcseife
ifdef MODEL
	$(MAKE) -C ../src/ denise MODEL=$(MODEL)
else
	$(MAKE) -C ../src/ denise 
endif


.PHONY: clean
clean: 
	rm -rf ../bin/denise
	$(MAKE) -C ../contrib/aff clean
	$(MAKE) -C ../contrib/fourier clean
	$(MAKE) -C ../contrib/stfinv clean
	$(MAKE) -C ../contrib/libcseife clean
	$(MAKE) -C ../src/ clean
