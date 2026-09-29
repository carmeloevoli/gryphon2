# Figures for the hpwne paper.
#
#   make            build every figure into ./figures
#   make publish    copy ./figures/*.pdf into the paper directory (project.toml)
#   make clean      remove ./figures

PYTHON ?= python3
export PYTHONPATH := $(CURDIR)/../../src:$(PYTHONPATH)

SCRIPTS := \
	plot_injection_spectrum.py \
	plot_spirals.py \
	plot_escape.py \
	plot_distributions.py \
	plot_proton_flux.py \
	plot_crab.py

.PHONY: all $(SCRIPTS) publish clean

all: $(SCRIPTS)

$(SCRIPTS):
	$(PYTHON) $@

publish:
	$(PYTHON) -c "import gryphon_plots as gp; gp.publish()"

clean:
	rm -rf figures
