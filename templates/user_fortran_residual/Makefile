# GNU make PREDEFINES FC as `f77` (origin: default), so a plain `FC ?= gfortran`
# would silently do nothing and the build would call a compiler nobody has.
# Override the built-in default, but still honour an exported FC or `make FC=ifx`.
ifeq ($(origin FC),default)
  FC := gfortran
endif
FFLAGS ?= -O2

# `residual` is what resasm.yml runs (as ./residual, relative to this folder).
# On Windows gfortran writes residual.exe from this same rule and ./residual
# still resolves to it, so the shipped resasm.yml works unchanged on both.
residual: residual.f90
	$(FC) $(FFLAGS) -o residual residual.f90

.PHONY: clean
clean:
	rm -f residual residual.exe *.o *.mod
