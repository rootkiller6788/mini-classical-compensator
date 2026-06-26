
# Makefile for mini-pid-theory
# PID Controller Theory Library and Test Suite
#
# Course alignment: MIT 6.302, Stanford ENGR105, Berkeley ME132,
#                    ETH 151-0591, Cambridge 3F2

CC = gcc
CFLAGS = -Wall -Wextra -std=c99 -O2 -g
LDFLAGS = -lm

# Include paths
INCLUDES = -Iinclude

# Source files
SRC = src/pid_core.c \
      src/pid_tuning.c \
      src/pid_advanced.c \
      src/pid_adaptive.c

# Object files
OBJ = $(SRC:.c=.o)

# Test file
TEST_SRC = tests/test_pid.c
TEST_BIN = tests/test_pid

# Example files
EXAMPLES = examples/example_dc_motor \
           examples/example_temperature \
           examples/example_level_control \
           examples/example_flow_control

# Library
LIB = libpidtheory.a

.PHONY: all clean test examples lib

all: lib test examples

lib: $(LIB)

$(LIB): $(OBJ)
	ar rcs $@ $^

%.o: %.c include/mini-pid-theory.h
	$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@

# Test target
test: $(TEST_BIN)
	./$(TEST_BIN)

$(TEST_BIN): $(TEST_SRC) $(OBJ)
	$(CC) $(CFLAGS) $(INCLUDES) $(TEST_SRC) $(OBJ) -o $@ $(LDFLAGS)

# Example targets
examples/example_dc_motor: examples/example_dc_motor.c $(OBJ)
	$(CC) $(CFLAGS) $(INCLUDES) $< $(OBJ) -o $@ $(LDFLAGS)

examples/example_temperature: examples/example_temperature.c $(OBJ)
	$(CC) $(CFLAGS) $(INCLUDES) $< $(OBJ) -o $@ $(LDFLAGS)

examples/example_level_control: examples/example_level_control.c $(OBJ)
	$(CC) $(CFLAGS) $(INCLUDES) $< $(OBJ) -o $@ $(LDFLAGS)

examples/example_flow_control: examples/example_flow_control.c $(OBJ)
	$(CC) $(CFLAGS) $(INCLUDES) $< $(OBJ) -o $@ $(LDFLAGS)

examples: $(EXAMPLES)

# Run examples
run-examples: examples
	@echo "=== DC Motor Control ==="
	./examples/example_dc_motor
	@echo "\n=== Temperature Control ==="
	./examples/example_temperature
	@echo "\n=== Level Control ==="
	./examples/example_level_control
	@echo "\n=== Flow Control ==="
	./examples/example_flow_control

clean:
	rm -f $(OBJ) $(TEST_BIN) $(EXAMPLES) $(LIB)
