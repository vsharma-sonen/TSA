################################################
# Makefile for TERRA standalone
################################################

####################################################
# Declaration of Variables
####################################################
#
.SILENT:
#
STDROOT      = ${PWD}
SRCDIR       = $(STDROOT)/src
OBJDIR       = $(STDROOT)/obj
WORKDIR      = $(STDROOT)/work
#
#
####################################################
#  include (or build) an appropriate compiler file
#     some compiler files are available in LOCAL
#     and can be copied to Fopts or directly used from there
####################################################

MACHINE := $(shell hostname | sed 's/[0-9]//g')
#include LOCAL/Fopts_$(MACHINE)
include Fopts

####################################################
#  Declaration of the Object Files
####################################################

#include ObjFiles

OBJS=   $(OBJDIR)/kind_parameters.o           \
        $(OBJDIR)/data_block_fields.o         \
        $(OBJDIR)/data_constants.o            \
        $(OBJDIR)/data_fields.o               \
        $(OBJDIR)/data_io.o                   \
	$(OBJDIR)/data_modelconfig.o          \
	$(OBJDIR)/data_parallel.o             \
        $(OBJDIR)/data_runcontrol.o           \
        $(OBJDIR)/meteo_utilities.o           \
        $(OBJDIR)/utilities.o                 \
        $(OBJDIR)/sfc_terra_data.o            \
        $(OBJDIR)/sfc_terra_init.o            \
        $(OBJDIR)/sfc_terra.o                 \
        $(OBJDIR)/sfc_utilities.o             \
        $(OBJDIR)/sfc_snow_data.o             \
        $(OBJDIR)/sfc_snow_utilities.o        \
        $(OBJDIR)/sfc_snow_init.o             \
        $(OBJDIR)/sfc_snow.o                  \
        $(OBJDIR)/turb_data.o                 \
	$(OBJDIR)/tsa_data.o                  \
	$(OBJDIR)/tsa_setup.o                 \
	$(OBJDIR)/tsa_input.o                 \
	$(OBJDIR)/tsa_output.o                \
	$(OBJDIR)/tsa_gribio.o                \
	$(OBJDIR)/tsa_lmparam.o               \
	$(OBJDIR)/tsa_interpol.o              \
	$(OBJDIR)/tsa_sfc_interface.o         \
	$(OBJDIR)/support_datetime.o          \
	$(OBJDIR)/fxtr_definition.o           \
	$(OBJDIR)/tsa_main.o

####################################################
#  Building the program
####################################################

terra : ${OBJS}
	$(LDSEQ) $(LDFLAGS) -o $(PROGRAM) $(OBJS) $(LIB)

clean:
	echo cleaning up
	rm -f $(PROGRAM)*
	rm -f $(OBJDIR)/*

####################################################
#  Dependencies 
####################################################

include ObjDependencies
