################################################################################
# Automatically-generated file. Do not edit!
################################################################################

RM := rm -f

# All of the sources participating in the build are defined here

#Library Sources

#Library Source Files
CPP_SRCS = \
	pwm.cpp \
	PDF.cpp \
	trellis.cpp \
	viterbi.cpp \
	stoch_viterbi.cpp \
	stoch_forward.cpp \
	nth_best.cpp \
	stochTable.cpp \
	backward.cpp \
	forward.cpp \
	baum_welch.cpp \
	forward_viterbi.cpp \
	posterior.cpp \
	traceback_path.cpp \
	externDefinitions.cpp \
	index.cpp \
	stochMath.cpp \
	text.cpp \
	userFunctions.cpp \
	hmm.cpp \
	state.cpp \
	lexicalTable.cpp \
	track.cpp \
	emm.cpp \
	externalFuncs.cpp \
	modelTemplate.cpp \
	transitions.cpp \
	weight.cpp \
	options.cpp \
	seqJobs.cpp \
	seqTracks.cpp \
	sequence.cpp \
	sequences.cpp 
 

OBJS = \
	pwm.o \
	PDF.o \
	trellis.o \
	viterbi.o \
	stoch_viterbi.o \
	stoch_forward.o \
	nth_best.o \
	stochTable.o \
	backward.o \
	forward.o \
	baum_welch.o \
	forward_viterbi.o \
	posterior.o \
	traceback_path.o \
	externDefinitions.o \
	index.o \
	stochMath.o \
	text.o \
	userFunctions.o \
	hmm.o \
	state.o \
	lexicalTable.o \
	track.o \
	emm.o \
	externalFuncs.o \
	modelTemplate.o \
	transitions.o \
	weight.o \
	options.o \
	seqJobs.o \
	seqTracks.o \
	sequence.o \
	sequences.o 


CPP_DEPS = \
	pwm.d \
	PDF.d \
	trellis.d \
	viterbi.d \
	stoch_viterbi.d \
	stoch_forward.d \
	nth_best.d \
	stochTable.d \
	backward.d \
	forward.d \
	baum_welch.d \
	forward_viterbi.d \
	posterior.d \
	traceback_path.d \
	externDefinitions.d \
	index.d \
	stochMath.d \
	text.d \
	userFunctions.d \
	hmm.d \
	state.d \
	lexicalTable.d \
	track.d \
	emm.d \
	externalFuncs.d \
	modelTemplate.d \
	transitions.d \
	weight.d \
	options.d \
	seqJobs.d \
	seqTracks.d \
	sequence.d \
	sequences.d 

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


# All Target
all: libStochHMM.a

# Tool invocations
#Create Static Library
libStochHMM.a: $(OBJS)
	@echo 'Building target: $@' 
	@echo 'Invoking: GCC Archiver'
	libtool -static  $(OBJS) -o libStochHMM.a
	ranlib libStochHMM.a
	chmod 644 libStochHMM.a
	@echo 'Finished building target: $@'
	@echo ' '



# Other Targets
clean:
	-$(RM) $(OBJS) $(CPP_DEPS) libStochHMM.a
	-@echo ' '

.PHONY: all clean dependents
.SECONDARY: