CXX = clang++
CXXFLAGS = -std=c++11 -fsanitize=address -Wall -Wextra -O3 -g 

# Paths
MCL_DIR = ./mcl
BUILD_DIR = ./build

# Include and Library paths
INCLUDES = -I$(MCL_DIR)/include
LDFLAGS = -L$(MCL_DIR)/lib -Wl,-rpath,$(shell pwd)/$(MCL_DIR)/lib
LIBS = -lmcl -lgmp -lgmpxx -lcrypto

# Source files
PLONK_SRC = ./src/plonk.cpp
NTT_SRC = ./src/ntt.cpp
KZG_SRC = ./src/kzg.cpp
HELPER_SRC = ./src/helper.cpp

# Test files
TEST_SRC = ./test.cpp
TEST_TARGET = $(BUILD_DIR)/test

# Default target
all: $(TEST_TARGET)

# Create build directory if it doesn't exist
$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

# Build the test executable
$(TEST_TARGET): $(TEST_SRC) $(PLONK_SRC) $(NTT_SRC) $(KZG_SRC) $(HELPER_SRC) | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -I./plonk -o $@ $^ $(LDFLAGS) $(LIBS)

# Run the test
test: $(TEST_TARGET)
	$(TEST_TARGET)

# Clean up
clean:
	rm -f $(TEST_TARGET)

.PHONY: all test clean
