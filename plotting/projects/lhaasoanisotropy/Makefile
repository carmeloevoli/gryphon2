PYTHON ?= python3
export PYTHONPATH := $(CURDIR)/../../src:$(PYTHONPATH)

SCRIPTS := plot_transport_vs_source.py

.PHONY: all $(SCRIPTS) publish clean

all: $(SCRIPTS)

$(SCRIPTS):
	$(PYTHON) $@

publish:
	$(PYTHON) -c "import gryphon_plots as gp; gp.publish()"

clean:
	rm -rf figures
