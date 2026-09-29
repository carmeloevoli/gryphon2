PYTHON ?= python3
SEEDS ?= 10000
export PYTHONPATH := $(CURDIR)/../../src:$(PYTHONPATH)

.PHONY: all test publish figure2 publish-figure2 snapshot publish-snapshot figure5 publish-figure5 halos publish-halos variations publish-variations low-rate publish-low-rate
.PHONY: results population-results publish-results publish-population
.PHONY: effective-sources publish-effective-sources
all: results

results: halos variations population-results

test:
	$(PYTHON) -m unittest test_feature_statistics.py test_figure2.py test_association_snapshot.py test_realizations.py test_injection_curvature.py test_halos.py test_variations.py test_low_rate.py test_population_results.py test_effective_sources.py
	$(PYTHON) -m unittest discover -s ../../../scripts -p 'test_run_*.py'

figure2:
	$(PYTHON) plot_escape.py
	$(PYTHON) plot_spirals.py

# Only the two Figure 2 panels are copied; association results are untouched.
publish-figure2: figure2
	$(PYTHON) -c "import gryphon_plots as gp; gp.publish(pattern='gryphon_escape.pdf'); gp.publish(pattern='gryphon_spirals.pdf')"

snapshot:
	$(PYTHON) plot_association_snapshot.py

publish-snapshot: snapshot
	$(PYTHON) -c "import gryphon_plots as gp; gp.publish(pattern='gryphon_association_snapshot.pdf')"

figure5:
	$(PYTHON) plot_realizations.py

publish-figure5: figure5
	$(PYTHON) -c "from pathlib import Path; import shutil; import gryphon_plots as gp; gp.publish(pattern='gryphon_realizations_H[24].pdf'); p=gp.project(); d=p.paper.parent/'modelstats'; d.mkdir(exist_ok=True); shutil.copy2(p.figures/'figure5_summary.json', d)"

halos:
	$(PYTHON) plot_halos.py --seeds $(SEEDS)

publish-halos: halos
	$(PYTHON) -c "import shutil; import gryphon_plots as gp; gp.publish(pattern='gryphon_halo_statistics.pdf'); p=gp.project(); d=p.paper.parent/'modelstats'; d.mkdir(exist_ok=True); [shutil.copy2(p.figures/n, d/n) for n in ('halos_summary.json', 'halos_statistics.txt', 'halos_table.tex')]"

variations:
	$(PYTHON) plot_variations.py --seeds $(SEEDS)

publish-variations: variations
	$(PYTHON) -c "import shutil; import gryphon_plots as gp; gp.publish(pattern='gryphon_variations_*.pdf'); p=gp.project(); d=p.paper.parent/'modelstats'; d.mkdir(exist_ok=True); [shutil.copy2(p.figures/n, d/n) for n in ('variations_summary.json', 'variations_statistics.txt', 'variations_table.tex')]"

low-rate: population-results

population-results:
	$(PYTHON) plot_population_results.py --seeds $(SEEDS)

publish-low-rate: low-rate
	$(PYTHON) -c "import shutil; import gryphon_plots as gp; gp.publish(pattern='gryphon_low_rate_statistics.pdf'); p=gp.project(); d=p.paper.parent/'modelstats'; d.mkdir(exist_ok=True); [shutil.copy2(p.figures/n, d/n) for n in ('low_rate_summary.json', 'low_rate_statistics.txt', 'low_rate_table.tex')]"

# Publishing is explicit; also retain the machine-readable statistics in the paper.
publish: population-results
	$(PYTHON) -c "import shutil; import gryphon_plots as gp; gp.publish(pattern='gryphon_association_statistics.pdf'); p=gp.project(); d=p.paper.parent/'modelstats'; d.mkdir(exist_ok=True); [shutil.copy2(p.figures/n, d/n) for n in ('association_summary.json', 'association_statistics.txt', 'association_table.tex')]"

publish-population: publish-low-rate publish

publish-results: publish-halos publish-variations publish-population

effective-sources:
	$(PYTHON) plot_effective_sources.py

publish-effective-sources: effective-sources
	$(PYTHON) -c "import shutil; import gryphon_plots as gp; gp.publish(pattern='gryphon_effective_sources.pdf'); p=gp.project(); d=p.paper.parent/'modelstats'; d.mkdir(exist_ok=True); shutil.copy2(p.figures/'effective_sources_summary.json', d)"
