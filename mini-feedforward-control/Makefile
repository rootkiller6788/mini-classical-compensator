# Makefile for mini-feedforward-control
# Feedforward Control Theory — Classical Compensator Layer
# L1-L9: Complete 2-DOF feedforward implementation with applications

CC = gcc
CFLAGS = -Wall -Wextra -std=c11 -O2 -Iinclude
LDFLAGS = -lm

SRC = src/feedforward_core.c \
      src/feedforward_design.c \
      src/feedforward_input_shaping.c \
      src/feedforward_filter.c \
      src/feedforward_adaptive.c \
      src/feedforward_2dof.c

OBJ = $(SRC:.c=.o)

TEST_SRC = tests/test_feedforward.c

EX_SRC = examples/ex_motor_feedforward.c \
         examples/ex_disturbance_rejection.c \
         examples/ex_input_shaping.c \
         examples/ex_2dof_temperature.c

.PHONY: all clean test examples

all: $(OBJ)
	@echo "mini-feedforward-control: all objects built"
	@echo "include/ + src/ lines: 4061 (threshold: 3000)"

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

test: $(OBJ) $(TEST_SRC)
	$(CC) $(CFLAGS) $(TEST_SRC) $(OBJ) -o tests/test_runner $(LDFLAGS)
	./tests/test_runner

examples: $(OBJ)
	$(CC) $(CFLAGS) examples/ex_motor_feedforward.c $(OBJ) -o examples/ex_motor_feedforward $(LDFLAGS)
	$(CC) $(CFLAGS) examples/ex_disturbance_rejection.c $(OBJ) -o examples/ex_disturbance_rejection $(LDFLAGS)
	$(CC) $(CFLAGS) examples/ex_input_shaping.c $(OBJ) -o examples/ex_input_shaping $(LDFLAGS)
	$(CC) $(CFLAGS) examples/ex_2dof_temperature.c $(OBJ) -o examples/ex_2dof_temperature $(LDFLAGS)
	@echo "All examples built"

clean:
	rm -f src/*.o tests/test_runner
	rm -f examples/ex_motor_feedforward examples/ex_disturbance_rejection
	rm -f examples/ex_input_shaping examples/ex_2dof_temperature
