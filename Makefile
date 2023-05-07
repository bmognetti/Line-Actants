COMPFLAGS     = -O -Wall -fcheck=all -g -fbacktrace # -fbounds-check # -fcheck=<all|bounds|array-temps> #  -fbounds-check

CFLAGS	      = $(COMPFLAGS)
PFLAGS	      = $(COMPFLAGS)
FFLAGS	      = $(COMPFLAGS)
CCFLAGS       = $(COMPFLAGS)
CXXFLAGS      = $(COMPFLAGS)

DEST	      = .

# Uncomment this if you want system header files to be expanded
#
# SYSHDRS       =

EXTHDRS	      =

HDRS	      =

INSTALL	      = install

LD	      = gfortran

LDFLAGS	      = $(COMPFLAGS)

LIBS	      =

LINTLIBS      =

LINTFLAGS     = -u $(CFLAGS)

MAKEFILE      = Makefile

OBJS	      = line_actants.o\
                thermal.o\
		load_mask.o

PRINT	      = pr

PRINTFLAGS    =

LP	      = lp

LPFLAGS       = 

PROGRAM       = line_actants

SHELL	      = /bin/sh

SRCS	      = line_actants.f\
                thermal.f\
		load_mask.f

all:		$(PROGRAM)

$(PROGRAM):     $(OBJS) $(LIBS) $(MAKEFILE)
		@echo "Linking $(PROGRAM) ..."
		@$(LD) $(LDFLAGS) $(OBJS) $(LIBS) -o $(PROGRAM)
		@echo "done"

objects:	line_actants.f thermal.f update_tiles.f load_mask.f
		@$(LD) -c $(LDFLAGS) line_actants.f $(LIBS) -o line_actants.o
		@$(LD) -c $(LDFLAGS) thermal.f $(LIBS) -o thermal.o
		@$(LD) -c $(LDFLAGS) load_mask.f $(LIBS) -o load_mask.o
