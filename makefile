################################################################################
# Automatically-generated file. Do not edit!
################################################################################

-include ../makefile.init

RM := rm -rf

# All of the sources participating in the build are defined here
O_SRCS := 
CPP_SRCS := 
C_UPPER_SRCS := 
C_SRCS := 
S_UPPER_SRCS := 
OBJ_SRCS := 
ASM_SRCS := 
CXX_SRCS := 
C++_SRCS := 
CC_SRCS := 
C++_DEPS := 
OBJS := 
C_DEPS := 
CC_DEPS := 
CPP_DEPS := 
EXECUTABLES := 
CXX_DEPS := 
C_UPPER_DEPS := 

# Every subdirectory with source files must be described here
SUBDIRS := \
src \

# Add inputs and outputs from these tool invocations to the build variables 
O_SRCS += \
./src/StochHMM.o \
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
./src/externDefinitions.o \
./src/externalFuncs.o \
./src/hmm.o \
./src/emm.o \
./src/index.o \
./src/lexicalTable.o \
./src/modelTemplate.o \
./src/options.o \
./src/seqJobs.o \
./src/seqTracks.o \
./src/sequence.o \
./src/sequences.o \
./src/state.o \
./src/stochMath.o \
./src/table.o \
./src/text.o \
./src/traceback_path.o \
./src/track.o \
./src/transitions.o \
./src/userFunctions.o \
./src/weight.o 

CPP_SRCS += \
./src/StochHMM.cpp \
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
./src/externDefinitions.cpp \
./src/externalFuncs.cpp \
./src/hmm.cpp \
./src/emm.cpp \
./src/index.cpp \
./src/lexicalTable.cpp \
./src/modelTemplate.cpp \
./src/options.cpp \
./src/seqJobs.cpp \
./src/seqTracks.cpp \
./src/sequence.cpp \
./src/sequences.cpp \
./src/state.cpp \
./src/stochMath.cpp \
./src/table.cpp \
./src/text.cpp \
./src/traceback_path.cpp \
./src/track.cpp \
./src/transitions.cpp \
./src/userFunctions.cpp \
./src/weight.cpp 

OBJS += \
./src/StochHMM.o \
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
./src/externDefinitions.o \
./src/externalFuncs.o \
./src/hmm.o \
./src/emm.o \
./src/index.o \
./src/lexicalTable.o \
./src/modelTemplate.o \
./src/options.o \
./src/seqJobs.o \
./src/seqTracks.o \
./src/sequence.o \
./src/sequences.o \
./src/state.o \
./src/stochMath.o \
./src/table.o \
./src/text.o \
./src/traceback_path.o \
./src/track.o \
./src/transitions.o \
./src/userFunctions.o \
./src/weight.o 

CPP_DEPS += \
./src/StochHMM.d \
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
./src/externDefinitions.d \
./src/externalFuncs.d \
./src/hmm.d \
./src/emm.d \
./src/index.d \
./src/lexicalTable.d \
./src/modelTemplate.d \
./src/options.d \
./src/seqJobs.d \
./src/seqTracks.d \
./src/sequence.d \
./src/sequences.d \
./src/state.d \
./src/stochMath.d \
./src/table.d \
./src/text.d \
./src/traceback_path.d \
./src/track.d \
./src/transitions.d \
./src/userFunctions.d \
./src/weight.d 


# Each subdirectory must supply rules for building sources it contributes
src/%.o: ../src/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


USER_OBJS :=

LIBS :=


ifneq ($(MAKECMDGOALS),clean)
ifneq ($(strip $(C++_DEPS)),)
-include $(C++_DEPS)
endif
ifneq ($(strip $(C_DEPS)),)
-include $(C_DEPS)
endif
ifneq ($(strip $(CC_DEPS)),)
-include $(CC_DEPS)
endif
ifneq ($(strip $(CPP_DEPS)),)
-include $(CPP_DEPS)
endif
ifneq ($(strip $(CXX_DEPS)),)
-include $(CXX_DEPS)
endif
ifneq ($(strip $(C_UPPER_DEPS)),)
-include $(C_UPPER_DEPS)
endif
endif

# Add inputs and outputs from these tool invocations to the build variables 

# All Target
all: bin_dir StochHMM libStochHMM.a

#Create directory
bin_dir:
	mkdir -p bin

# Tool invocations
StochHMM: $(OBJS) $(USER_OBJS)
	@echo 'Building target: $@'
	@echo 'Invoking: C++ Linker'
	g++  -o "bin/StochHMM" $(OBJS) $(USER_OBJS) $(LIBS)
	@echo 'Finished building target: $@'
	@echo ' '

libStochHMM.a: $(OBJS) $(USER_OBJS)
	@echo 'Building target: $@'
	@echo 'Invoking: GCC Archiver'
	ar -r  "bin/libStochHMM.a" $(OBJS) $(USER_OBJS) $(LIBS)
	@echo 'Finished building target: $@'
	@echo ' '

# Other Targets
clean:
	-$(RM) $(C++_DEPS)$(OBJS)$(C_DEPS)$(CC_DEPS)$(CPP_DEPS)$(EXECUTABLES)$(CXX_DEPS)$(C_UPPER_DEPS) bin/
	-@echo ' '

.PHONY: all clean dependents
.SECONDARY:

