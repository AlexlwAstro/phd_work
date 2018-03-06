SRCS := $(shell ls *.f)
  OBJS := $(SRCS:%.f=%.o)

  FC     := ifort
  FFLAGS := -O3 -f77rtl -ipo -ftz -save -warn noalign

  franec2p1p4.go: $(OBJS)
	$(FC) $(FFLAGS) $^ -o $@

  $(OBJS): $(SRCS)
	$(FC) $(FFLAGS) -c $^

  clean:
	rm -f $(OBJS)

  info:
	@echo $(OBJS)
	@echo $(SRCS)
