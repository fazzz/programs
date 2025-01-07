SHELL=/bin/sh
AMBER=Amber12
#CONFIG_FILE=../config.h
CONFIG_FILE=config.h
CONFIG_COMMAND=./configure
PARALLEL_TOKEN=DMPI
CUDA_TOKEN=DCUDA
WINDOWS_TOKEN=DWINDOWS

include $(CONFIG_FILE)

install: $(INSTALLTYPE)

serial: configured_serial
	$(MAKE) -C src/ install
	@echo "Installation of pmemd complete"

parallel: configured_parallel
	$(MAKE) -C src/ parallel
	@echo "Installation of pmemd.MPI complete"

cuda: configured_cuda
	$(MAKE) -C src/ cuda
	@echo "Installation of pmemd.cuda complete"

cuda_parallel: configured_cuda configured_parallel
	$(MAKE) -C src/ cuda_parallel
	@echo "Installation of pmemd.cuda.MPI complete"

#Cannot use configured_win or configured_serial for Windows right now
#since this will not run in a dos shell.
#install_win: configured_win configured_serial
#parallel_win: configured_win configured_parallel
install_win:
	$(MAKE) -C src/ install_win

parallel_win:
	$(MAKE) -C src/ parallel_win

uninstall.serial: clean
	-rm -f $(BINDIR)/pmemd$(SFX)

uninstall.cuda: clean
	-rm -f $(BINDIR)/pmemd.cuda$(SFX)
	-rm -f $(BINDIR)/pmemd.cuda_SPFP$(SFX)
	-rm -f $(BINDIR)/pmemd.cuda_SPDP$(SFX)
	-rm -f $(BINDIR)/pmemd.cuda_DPDP$(SFX)

uninstall.parallel: clean
	-rm -f $(BINDIR)/pmemd.MPI$(SFX)

depends:
	$(MAKE) -C src/ depends

clean:
	$(MAKE) -C src/ clean

configured:
	@(if [ ! -f $(CONFIG_FILE) ] ; then \
	    echo "Error: $(CONFIG_COMMAND) must be executed before $(MAKE)!" ;\
	    exit 2 ;\ # $(CONFIG_COMMAND) ;\
	fi ;\
	)

configured_parallel: configured
	@(grep $(PARALLEL_TOKEN) $(CONFIG_FILE) > /dev/null || \
	    { echo "Error: $(CONFIG_FILE) is not of type parallel!" ;\
	      echo "  Rerun $(CONFIG_COMMAND) and specify an MPI implementation." ;\
	      exit 2 ;\
	    } ;\
	)

configured_serial: configured
	@(if grep $(PARALLEL_TOKEN) $(CONFIG_FILE) > /dev/null ; then \
	    echo "Error: $(CONFIG_FILE) is of type parallel, not serial!" ;\
	    echo "  Rerun $(CONFIG_COMMAND) and specify -nopar." ;\
	    exit 2 ;\
	fi ;\
	)
	@(if grep $(CUDA_TOKEN) $(CONFIG_FILE) > /dev/null ; then \
	    echo "Error: $(CONFIG_FILE) is of type cuda, not serial!" ;\
	    echo "  Rerun $(CONFIG_COMMAND) and specify -nopar." ;\
	    exit 2 ;\
	fi ;\
	)

configured_cuda: configured
	@(grep $(CUDA_TOKEN) $(CONFIG_FILE) > /dev/null || \
	    { echo "Error: $(CONFIG_FILE) is not of type cuda!" ;\
	      echo "  Rerun $(CONFIG_COMMAND) and specify a cuda implementation." ;\
	      exit 2 ;\
	    } ;\
	)

configured_win: configured
	@(grep $(WINDOWS_TOKEN) $(CONFIG_FILE) > /dev/null || \
	{ echo "Error: $(CONFIG_FILE) is not of type Windows!" ;\
		echo "  Rerun $(CONFIG_COMMAND) and specify -windows." ;\
		exit 2 ;\
	} ;\
	)


