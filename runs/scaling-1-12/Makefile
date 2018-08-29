

default: IFOS_INV

PHONY: IFOS_FW
IFOS_FW:
	$(MAKE) -C ../../src/ IFOS2D EXEC="../runs/scaling-1-12" MODEL_AC=../genmod/toy_example_ac_true.c
	mv IFOS2D IFOS_FW

IFOS_INV: IFOS_FW
	$(MAKE) -C ../../src/ IFOS2D EXEC="../runs/scaling-1-12" MODEL_AC=../genmod/toy_example_ac_start.c
	mv IFOS2D IFOS_INV

.PHONY: clean
clean: clean_local clean_src

.PHONY: clean_local
clean_local:
	rm -rf IFOS_FW
	rm -rf IFOS_INV

.PHONY: clean_src
clean_src:
	rm -rf ../../bin/tmp/IFOS2D
	$(MAKE) -C ../../src clean

.PHONY: clean_data
clean_data:
	rm -rf FW_model/*
	rm -rf INV_model/*
	rm -rf seis_FW/*
	rm -rf seis_INV/*
	rm -rf jacobian/*
	rm -rf taper_te_ac*
