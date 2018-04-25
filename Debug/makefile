-include ../makefile.init

RM := rm -rf

#CXX := # uncomment and add a specific compiler path
BCRCPPLAROOT := /path/to/#BCR_CPP_LA should be found in this directory

# All of the sources participating in the build are defined here
-include sources.mk
-include src/subdir.mk
-include subdir.mk
-include objects.mk

ifneq ($(MAKECMDGOALS),clean)
ifneq ($(strip $(CC_DEPS)),)
-include $(CC_DEPS)
endif
ifneq ($(strip $(C++_DEPS)),)
-include $(C++_DEPS)
endif
ifneq ($(strip $(C_UPPER_DEPS)),)
-include $(C_UPPER_DEPS)
endif
ifneq ($(strip $(CXX_DEPS)),)
-include $(CXX_DEPS)
endif
ifneq ($(strip $(CPP_DEPS)),)
-include $(CPP_DEPS)
endif
ifneq ($(strip $(C_DEPS)),)
-include $(C_DEPS)
endif
endif

-include ../makefile.defs

# Add inputs and outputs from these tool invocations to the build variables 

# All Target
all: DiscreteCCSOpt.d

# Tool invocations
DiscreteCCSOpt.d: $(DOBJS) $(USER_DOBJS)
	@echo 'Building target: $@'
	@echo 'Invoking: GCC C++ Linker'
	$(CXX) -L/usr/local/lib -o "DiscreteCCSOpt.d" -g $(DOBJS) $(USER_DOBJS) $(LIBS)
	@echo 'Finished building target: $@'
	@echo ' '


# Other Targets
clean:
	-$(RM) $(CC_DEPS)$(CC_DDEPS)$(C++_DEPS)$(EXECUTABLES)$(C_UPPER_DEPS)$(CXX_DEPS)$(OBJS)$(DOBJS)$(CPP_DEPS)$(C_DEPS) DiscreteCCSopt.d
	-@echo ' '

docs:
	@echo 'Building documentation'
	doxygen Doxyfile

.PHONY: all clean dependents docs
.SECONDARY:

-include ../makefile.targets
