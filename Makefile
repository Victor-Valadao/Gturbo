NSE 	 = src/nse/main
LYA 	 = src/lya/main
LAG 	 = src/lag/main
PAS		 = src/pas/main
FC       = nvfortran
EXE	     = exe

FCFLAGS  = -fast -acc=gpu -cuda -cudalib=cufft \
			-Minfo=all -Mpreprocess -mcmodel=medium

SRCFILESN=src/comm/modules.f90 src/comm/defk.f90 \
          src/comm/iniforcing.f90 src/comm/physics.f90 \
          src/comm/diagnostics.f90 src/comm/rkstep.f90 \
          src/comm/forcing.f90 src/comm/readwrite.f90 \
		  src/nse/inieuler.f90 src/comm/fft.f90 \
		  src/comm/sf_3rd.f90 \
		  
SRCFILESL=src/comm/modules.f90 src/comm/defk.f90 \
          src/comm/iniforcing.f90 src/comm/physics.f90 \
          src/comm/diagnostics.f90 src/comm/diagnosticsp.f90 \
          src/lya/diagnosticsd.f90 src/comm/rkstep.f90 \
          src/comm/forcing.f90 src/comm/readwrite.f90 \
		  src/lya/inieuler.f90 src/lya/rescale.f90 \
		  src/comm/fft.f90 \
		  
SRCFILESP=src/comm/modules.f90 src/comm/defk.f90 \
          src/comm/iniforcing.f90 src/pas/physicp.f90 \
          src/comm/diagnostics.f90 src/comm/diagnosticsp.f90 \
          src/pas/rkptep.f90 src/comm/forcing.f90 \
          src/comm/readwrite.f90 src/pas/inieuler.f90 \
          src/comm/fft.f90 \
          
SRCFILESA=src/comm/modules.f90 src/comm/defk.f90 \
          src/comm/iniforcing.f90 src/lag/physicl.f90 \
          src/comm/diagnostics.f90 src/lag/diaglag.f90 \
          src/lag/advection.f90 src/lag/rkltep.f90 \
          src/comm/forcing.f90 src/comm/readwrite.f90 \
		  src/lag/inieuler.f90 src/lag/rescale.f90 \
		  src/comm/fft.f90 src/lag/inilag.f90 \

nse: $(NSE).f90
	$(FC) $(FCFLAGS) -o $(NSE).$(EXE) $(SRCFILESN) $<
	@mkdir -p gturb gturb/fields gturb/files
	@scp src/comm/jobscript gturb/jobscript
	@scp src/nse/main.exe gturb/nse.exe
	@scp src/comm/seed.0 gturb/files/seed.000
	@scp src/comm/curframe.dat gturb/curframe.dat
	@scp src/nse/params.dat gturb/params.dat
	@scp src/reset.sh gturb/reset.sh
	
	@echo 'Cleaning up...'
	@rm -rf *.o *.mod
	@rm -rf src/nse/main.exe

lya: $(LYA).f90
	$(FC) $(FCFLAGS) -o $(LYA).$(EXE) $(SRCFILESL) $<
	@mkdir -p gturb gturb/fields gturb/files
	@scp src/comm/jobscript gturb/jobscript
	@scp src/lya/main.exe gturb/lya.exe
	@scp src/comm/seed.0 gturb/files/seed.000
# 	@scp w.000.060 gturb/fields/w.000
	@scp src/comm/curframe.dat gturb/curframe.dat
	@scp src/lya/params.dat gturb/params.dat
	@scp src/reset.sh gturb/reset.sh
	
	@echo 'Cleaning up...'
	@rm -rf *.o *.mod
	@rm -rf src/lya/main.exe

pas: $(PAS).f90
	$(FC) $(FCFLAGS) -o $(PAS).$(EXE) $(SRCFILESP) $<
	@mkdir -p gturb gturb/fields gturb/files
	@scp src/comm/jobscript gturb/jobscript
	@scp src/pas/main.exe gturb/pas.exe
	@scp src/comm/seed.0 gturb/files/seed.000
	@scp src/comm/curframe.dat gturb/curframe.dat
	@scp src/pas/params.dat gturb/params.dat
	@scp src/reset.sh gturb/reset.sh
	
	@echo 'Cleaning up...'
	@rm -rf *.o *.mod
	@rm -rf src/pas/main.exe
	
lag: $(LAG).f90
	$(FC) $(FCFLAGS) -o $(LAG).$(EXE) $(SRCFILESA) $<
	@mkdir -p gturb gturb/fields gturb/files
	@scp src/comm/jobscript gturb/jobscript
	@scp src/lag/main.exe gturb/lag.exe
	@scp src/comm/seed.0 gturb/files/seed.000
	@scp src/comm/curframe.dat gturb/curframe.dat
	@scp src/lag/params.dat gturb/params.dat
	@scp src/reset.sh gturb/reset.sh
	
	@echo 'Cleaning up...'
	@rm -rf *.o *.mod
	@rm -rf src/lya/main.exe
	

clean:
	@echo 'Cleaning up...'
	@rm -rf *.$(EXE) *.o *.dwf *.pdb *.mod prof