# Compiler and flags
FC = gfortran                   # Fortran compiler
FFLAGS = -w -O2 -J$(BUILD_DIR)  # Compiler flags:
                                # -w    : Disable all warnings
                                # -O2   : Optimization level 2
                                # -J    : Direct .mod files to build directory

# Directories
SRC_DIR = src
MODULES_DIR = modules
BUILD_DIR = build

# Source module files
MODULES = $(MODULES_DIR)/MOD_Select_Kind.f90 \
          $(MODULES_DIR)/MOD_Roots.f90		  

# Object files for modules
MODULE_OBJS = $(patsubst $(MODULES_DIR)/%.f90,$(BUILD_DIR)/%.o,$(MODULES))

# Executable names
EXES = main.exe test.exe

# Source files for executables
MAIN_SRC = $(SRC_DIR)/main.f90
TEST_SRC = $(SRC_DIR)/test.f90

# Object files for executables
MAIN_OBJ = $(BUILD_DIR)/main.o
TEST_OBJ = $(BUILD_DIR)/test.o

# Default target: build all executables
all: $(EXES)

# Rule to build main.exe
main.exe: $(MODULE_OBJS) $(MAIN_OBJ)
	$(FC) $(FFLAGS) -o $@ $^

# Rule to build test.exe
test.exe: $(MODULE_OBJS) $(TEST_OBJ)
	$(FC) $(FFLAGS) -o $@ $^

# Generic rule to compile module files
$(BUILD_DIR)/%.o: $(MODULES_DIR)/%.f90 | $(BUILD_DIR)
	$(FC) $(FFLAGS) -c $< -o $@

# Generic rule to compile source files (main and test)
$(BUILD_DIR)/%.o: $(SRC_DIR)/%.f90 | $(BUILD_DIR)
	$(FC) $(FFLAGS) -c $< -o $@

# Ensure build directory exists
$(BUILD_DIR):
	mkdir $(BUILD_DIR)

# Clean target to remove build artifacts and executables
clean:
	if exist $(BUILD_DIR) rmdir /s /q $(BUILD_DIR)
	if exist main.exe del /f /q main.exe
	if exist main_debug.exe del /f /q main_debug.exe
	if exist test.exe del /f /q test.exe
	del /f /q *.o *.mod > NUL 2>&1

# Phony targets to avoid conflicts with files named 'all', 'clean', etc.
.PHONY: all clean main.exe test.exe
