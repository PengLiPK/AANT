fc=gfortran


all: linklist.o strct.o fmm_adpt fmm_adpt_del fmm_forward \
	 inverseprb prepsr


fmm_adpt : fmm_adpt.o strct.o linklist.o
	$(fc) -o fmm_adpt fmm_adpt.o strct.o linklist.o
fmm_adpt.o : fmm_adpt.f90
	$(fc) -c fmm_adpt.f90
fmm_adpt_del : fmm_adpt_del.o strct.o linklist.o
	$(fc) -o fmm_adpt_del fmm_adpt_del.o strct.o linklist.o
fmm_adpt_del.o : fmm_adpt_del.f90
	$(fc) -c fmm_adpt_del.f90
fmm_forward : fmm_forward.o strct.o
	$(fc) -o fmm_forward fmm_forward.o strct.o
fmm_forward.o : fmm_forward.f90
	$(fc) -c fmm_forward.f90
inverseprb : inverseprb.o
	$(fc) -o inverseprb inverseprb.o
inverseprb.o : inverseprb.f90
	$(fc) -c inverseprb.f90
prepsr : prepsr.o
	$(fc) -o prepsr prepsr.o
prepsr.o : prepsr.f90
	$(fc) -c prepsr.f90
interpvel : interpvel.o strct.o
	$(fc) -o interpvel interpvel.o strct.o
interpvel.o : interpvel.f90
	$(fc) -c interpvel.f90
strct.o : strct.f90
	$(fc) -c strct.f90
linklist.o : linklist.f90
	$(fc) -c linklist.f90
