################################################################################
# Automatically-generated file. Do not edit!
################################################################################

-include ../makefile.init

RM := rm -rf

# Every subdirectory with source files must be described here
SUBDIRS := src 


# All of the sources participating in the build are defined here
#Library Sources
STOCHHMMLIB_SRCS := 
STOCHHMMLIB_CPP_SRCS := 
STOCHHMMLIB_OBJS :=
STOCHHMMLIB_CPP_DEPS :=

#Application Source Files
STOCHHMM_SRCS := ./src/StochHMM.o
STOCHHMM_CPP_SRCS := ./src/StochHMM.cpp
STOCHHMM_OBJS := ./src/StochHMM.o
STOCHHMM_CPP_DEPS := ./src/StochHMM.d
STOCHHMM_LIBS := ./bin/libStochHMM.a


#Library Source Files
STOCHHMMLIB_SRCS += \
./src/PDF.o \
./src/trellis.o \
./src/viterbi.o \
./src/stoch_viterbi.o \
./src/stoch_forward.o \
./src/nth_best.o \
./src/stochTable.o \
./src/backward.o \
./src/forward.o \
./src/baum_welch.o \
./src/forward_viterbi.o \
./src/posterior.o \
./src/traceback_path.o \
./src/externDefinitions.o \
./src/index.o \
./src/stochMath.o \
./src/text.o \
./src/userFunctions.o \
./src/hmm.o \
./src/state.o \
./src/lexicalTable.o \
./src/track.o \
./src/emm.o \
./src/externalFuncs.o \
./src/modelTemplate.o \
./src/transitions.o \
./src/weight.o \
./src/options.o \
./src/seqJobs.o \
./src/seqTracks.o \
./src/sequence.o \
./src/sequences.o \


STOCHHMMLIB_CPP_SRCS += \
./src/PDF.cpp \
./src/trellis.cpp \
./src/viterbi.cpp \
./src/stoch_viterbi.cpp \
./src/stoch_forward.cpp \
./src/nth_best.cpp \
./src/stochTable.cpp \
./src/backward.cpp \
./src/forward.cpp \
./src/baum_welch.cpp \
./src/forward_viterbi.cpp \
./src/posterior.cpp \
./src/traceback_path.cpp \
./src/externDefinitions.cpp \
./src/index.cpp \
./src/stochMath.cpp \
./src/text.cpp \
./src/userFunctions.cpp \
./src/hmm.cpp \
./src/state.cpp \
./src/lexicalTable.cpp \
./src/track.cpp \
./src/emm.cpp \
./src/externalFuncs.cpp \
./src/modelTemplate.cpp \
./src/transitions.cpp \
./src/weight.cpp \
./src/options.cpp \
./src/seqJobs.cpp \
./src/seqTracks.cpp \
./src/sequence.cpp \
./src/sequences.cpp 
 

STOCHHMMLIB_OBJS += \
./src/PDF.o \
./src/trellis.o \
./src/viterbi.o \
./src/stoch_viterbi.o \
./src/stoch_forward.o \
./src/nth_best.o \
./src/stochTable.o \
./src/backward.o \
./src/forward.o \
./src/baum_welch.o \
./src/forward_viterbi.o \
./src/posterior.o \
./src/traceback_path.o \
./src/externDefinitions.o \
./src/index.o \
./src/stochMath.o \
./src/text.o \
./src/userFunctions.o \
./src/hmm.o \
./src/state.o \
./src/lexicalTable.o \
./src/track.o \
./src/emm.o \
./src/externalFuncs.o \
./src/modelTemplate.o \
./src/transitions.o \
./src/weight.o \
./src/options.o \
./src/seqJobs.o \
./src/seqTracks.o \
./src/sequence.o \
./src/sequences.o 


STOCHHMMLIB_CPP_DEPS += \
./src/PDF.d \
./src/trellis.d \
./src/viterbi.d \
./src/stoch_viterbi.d \
./src/stoch_forward.d \
./src/nth_best.d \
./src/stochTable.d \
./src/backward.d \
./src/forward.d \
./src/baum_welch.d \
./src/forward_viterbi.d \
./src/posterior.d \
./src/traceback_path.d \
./src/externDefinitions.d \
./src/index.d \
./src/stochMath.d \
./src/text.d \
./src/userFunctions.d \
./src/hmm.d \
./src/state.d \
./src/lexicalTable.d \
./src/track.d \
./src/emm.d \
./src/externalFuncs.d \
./src/modelTemplate.d \
./src/transitions.d \
./src/weight.d \
./src/options.d \
./src/seqJobs.d \
./src/seqTracks.d \
./src/sequence.d \
./src/sequences.d 
 
	

# Each subdirectory must supply rules for building sources it contributes
src/%.o: ../src/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


ifneq ($(MAKECMDGOALS),clean)
ifneq ($(strip $(CPP_DEPS)),)
-include $(STOCHHMMLIB_CPP_DEPS)
endif
ifneq ($(strip $(STOCHHMM_CPP_DEPS)),)
-include $(CPP_DEPS)
endif
endif

# Add inputs and outputs from these tool invocations to the build variables 

# All Target
all: bin_dir libStochHMM.a StochHMM

#Create directory
bin_dir:
	mkdir -p bin

# Tool invocations
#Create Static Library
libStochHMM.a: $(STOCHHMMLIB_OBJS)
	@echo 'Building target: $@' 
	@echo 'Invoking: GCC Archiver'
	ar -r  "bin/libStochHMM.a" $(STOCHHMMLIB_OBJS) $(STOCHHMMLIB_LIBS)
	@echo 'Finished building target: $@'
	@echo ' '

#Compile StochHMM application using Static Library
StochHMM: $(STOCHHMM_OBJS)
	@echo 'Building target: $@'
	@echo 'Invoking: C++ Linker'
	g++  -o "bin/StochHMM" $(STOCHHMM_OBJS) $(STOCHHMM_LIBS)
	@echo 'Finished building target: $@'
	@echo ' '


# Other Targets
clean:
	-$(RM) $(STOCHHMMLIB_SRCS)$(STOCHHMMLIB_CPP_OBJS)$(STOCHHMMLIB_CPP_DEPS)$(STOCHHMM_SRCS)$(STOCHHMM_OBJS)$(STOCHHMM_CPP_DEPS) bin/
	-@echo ' '

.PHONY: all clean dependents
.SECONDARY:

