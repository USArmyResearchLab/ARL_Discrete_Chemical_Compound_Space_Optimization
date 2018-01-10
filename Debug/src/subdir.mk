################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CC_SRCS += \
../src/DiscreteCCSOpt.cc \
../src/Library_data.cc \
../src/chem_opt.cc \
../src/chemgroup.cc \
../src/chemident.cc \
../src/entropic_aux.cc \
../src/generalbaseiterator.cc \
../src/parse.cc \
../src/simpleprune.cc \
../src/zmat.cc \
../src/zmat_opt.cc 

CC_DEPS += \
./src/DiscreteCCSOpt.d \
./src/Library_data.d \
./src/chem_opt.d \
./src/chemgroup.d \
./src/chemident.d \
./src/entropic_aux.d \
./src/generalbaseiterator.d \
./src/parse.d \
./src/simpleprune.d \
./src/zmat.d \
./src/zmat_opt.d 

OBJS += \
./src/DiscreteCCSOpt.o \
./src/Library_data.o \
./src/chem_opt.o \
./src/chemgroup.o \
./src/chemident.o \
./src/entropic_aux.o \
./src/generalbaseiterator.o \
./src/parse.o \
./src/simpleprune.o \
./src/zmat.o \
./src/zmat_opt.o 

DOBJS += \
./src/DiscreteCCSOpt.d.o \
./src/Library_data.d.o \
./src/chem_opt.d.o \
./src/chemgroup.d.o \
./src/chemident.d.o \
./src/entropic_aux.d.o \
./src/generalbaseiterator.d.o \
./src/parse.d.o \
./src/simpleprune.d.o \
./src/zmat.d.o \
./src/zmat_opt.d.o



# Each subdirectory must supply rules for building sources it contributes
src/%.o: ../src/%.cc
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	/usr/local/bin/g++ -g -DDEBUG -I/usr/local/include -I/home/berend/projects/BCR_CPP_LA -I/home/berend/projects -I"/home/berend/workspace/DiscreteOptimize/include" -I"/home/berend/workspace/DiscreteOptimize/src" -I"../" -O3 -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


