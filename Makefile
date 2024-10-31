TEST 	 = src/main
# FC       = mpif90
FC       = nvfortran
EXE	 = exe

# FCFLAGS  = -fast -gpu=managed,cc75,cuda11.5.0 -acc=gpu -cuda -cudalib=cufft -Minfo=accel -Mpreprocess -D_USE_NVTX -lnvToolsExt -mcmodel=medium
FCFLAGS  = -fast -acc=gpu -cuda -cudalib=cufft -Minfo=all -Mpreprocess -mcmodel=medium#-fPIC
# FCFLAGS  = -fast -acc=gpu -cuda -stdpar=gpu -gpu=cc86 -cudalib=cufft -mcmodel=medium -Minfo=all

# -ta=tesla:managed
# FCFLAGS  = -acc=gpu -cuda -cudalib=cufft -mcmodel=medium -Minfo=all

SRCFILES = src/modules.f90 src/defk.f90 src/iniforcing.f90 src/trunc.f90 \
					 src/nlt.f90 src/produ.f90 src/diagnostics.f90 src/rkstep.f90 \
					 src/forcing.f90 src/read.f90 src/write.f90 src/inieuler.f90 \
					 src/fft_inv.f90 src/fft_dir.f90 \

all: build run verify

build: $(TEST).f90
	$(FC) $(FCFLAGS) -o $(TEST).$(EXE) $(SRCFILES) $<
	@mkdir -p gturb gturb/fields gturb/files
# 	@mkdir -p gturb/Diag/fields gturb/Diag/spectra gturb/Diag/fluxes
# 	@mkdir -p gturb/history-slurm
	@scp src/jobscript gturb/jobscript
	@scp src/main.exe gturb/main.exe
	@scp src/seed.0 gturb/files/seed.000
	@scp src/startframe.dat gturb/curframe.dat
	@scp src/params.dat gturb/params.dat
	@scp src/reset.sh gturb/reset.sh
	
	@echo 'Cleaning up...'
	@rm -rf *.o *.mod
	@rm -rf src/main.exe
	
run:
	# $(TEST).$(EXE)
# 	$(RUN) ./$(TEST).$(EXE)

verify:


clean:
	@echo 'Cleaning up...'
	@rm -rf *.$(EXE) *.o *.dwf *.pdb *.mod prof