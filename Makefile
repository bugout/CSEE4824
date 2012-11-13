TEAM ?= none

PREFIX := w4824
DIRNAME := $(PREFIX)-team$(TEAM)
TARNAME := $(DIRNAME).tar.bz2
PROGRAM := sommelier
REPORT := report
RESULTS := results
TARBALL := tarball

ifeq "$(TEAM)" "none"
tarball:
	@echo "Error: forgot to set TEAM? [do 'make tarball TEAM=X']" && exit 1
else
tarball: tarclean
	mkdir -p $(DIRNAME)
	make -C $(PROGRAM) clean
	cp -r $(PROGRAM) $(DIRNAME)
	cp $(REPORT)/report.pdf $(DIRNAME)
	cp -r $(RESULTS) $(DIRNAME)
	cp README* $(DIRNAME)
	-cp -r .git $(DIRNAME)
	tar -cjf $(TARNAME) $(DIRNAME)
endif

tarclean:
	 rm -rf $(DIRNAME) $(TARNAME)

.PHONY: tarclean tarball
