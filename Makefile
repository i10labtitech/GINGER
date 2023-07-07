all: mapping denovo homology abinitio merge_phase0 merge_phase1 merge_phase2 summary evaluation

clean: clean_mapping clean_denovo clean_homology clean_abinitio \
clean_merge_phase0 clean_merge_phase1 clean_merge_phase2 clean_summary clean_evaluation

mapping:
	cd src/mapping && $(MAKE) all

denovo:
	cd src/denovo && $(MAKE) all

homology:
	cd src/homology && $(MAKE) all

abinitio:
	cd src/abinitio && $(MAKE) all

merge_phase0:
	cd src/merge_phase0 && $(MAKE) all

merge_phase1:
	cd src/merge_phase1 && $(MAKE) all

merge_phase2:
	cd src/merge_phase2 && $(MAKE) all

summary:
	cd src/summary && $(MAKE) all

evaluation:
	cd src/evaluation && $(MAKE) all

clean_mapping:
	cd src/mapping && $(MAKE) clean

clean_denovo:
	cd src/denovo && $(MAKE) clean

clean_homology:
	cd src/homology && $(MAKE) clean

clean_abinitio:
	cd src/abinitio && $(MAKE) clean

clean_merge_phase0:
	cd src/merge_phase0 && $(MAKE) clean

clean_merge_phase1:
	cd src/merge_phase1 && $(MAKE) clean

clean_merge_phase2:
	cd src/merge_phase2 && $(MAKE) clean

clean_summary:
	cd src/summary && $(MAKE) clean

clean_evaluation:
	cd src/evaluation && $(MAKE) clean
