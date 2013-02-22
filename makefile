################################################################################
# Automatically-generated file. Do not edit!
################################################################################

-include ../makefile.init

RM := rm -f

# Every subdirectory with source files must be described here
SUBDIRS := src 

#Application Source Files
CPP_SRCS := src/StochHMM.cpp 
CPP_DEPS := src/StochHMM.d 
OBJS := src/StochHMM.o
STOCHHMM_LIBS := src/libStochHMM.a 


# Each subdirectory must supply rules for building sources it contributes
%.o: %.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


ifneq ($(MAKECMDGOALS),clean)
ifneq ($(strip $(CPP_DEPS)),)
-include $(CPP_DEPS)
endif
ifneq ($(strip $(CPP_DEPS)),)
-include $(CPP_DEPS)
endif
endif

# Add inputs and outputs from these tool invocations to the build variables 

# All Target
all: bin_dir library StochHMM_obj StochHMM

#Create directory
bin_dir:
	mkdir -p bin

# Tool invocations
#Create Static Library
library:
	cd src; ${MAKE} libStochHMM.a


StochHMM_obj:
	@echo 'Invoking: GCC C++ Compiler'
	#g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"src/StochHMM.d" -MT"src/StochHMM.d" -o "src/StochHMM.o" "src/StochHMM.cpp"
	@echo 'Finished building'

#Compile StochHMM application using Static Library
StochHMM: $(OBJS)	
	@echo 'Building target: $@'
	@echo 'Invoking: C++ Linker'
	g++  -o "bin/StochHMM" $(OBJS) $(STOCHHMM_LIBS)
	@echo 'Finished building target: $@'
	@echo ' '


# Other Targets
clean:
	-$(RM) $(OBJS) $(CPP_DEPS) bin/StochHMM
	cd src; ${MAKE} clean;
	-@echo ' '

.PHONY: all clean dependents
.SECONDARY:

