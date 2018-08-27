default: IFOS_INV

.phony: IFOS_FW
IFOS_FW:
	$(MAKE) -C ../../src/ IFOS2D EXEC="../bin/tmp" MODEL_AC=../genmod/toy_example_ac_true.c
	mv ../../bin/tmp/IFOS2D IFOS_FW

.phony: IFOS_INV
IFOS_INV: IFOS_FW
	$(MAKE) -C ../../src/ IFOS2D EXEC="../bin/tmp" MODEL_AC=../genmod/toy_example_ac_start.c
	mv ../../bin/tmp/IFOS2D IFOS_INV

clean: 
	rm -rf IFOS_FW
	rm -rf IFOS_INV
	rm -rf ../../bin/tmp/IFOS2D
	$(MAKE) -C ../../src clean
