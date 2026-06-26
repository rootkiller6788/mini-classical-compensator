# Makefile for mini-lag-compensator
# Classical Control -- Lag Compensator Design and Analysis
# L1-L9: Complete engineering implementation

CC = gcc
CFLAGS = -Wall -Wextra -std=c11 -O2 -Iinclude
LDFLAGS = -lm

SRC = src/lag_compensator.c       src/lag_design.c       src/lag_frequency.c       src/lag_simulation.c       src/lag_optimization.c       src/lag_digital.c       src/lag_application.c       src/lag_identification.c

OBJ = $(SRC:.c=.o)

TEST_SRC = tests/test_lag.c
EX_SRC = examples/ex_dc_motor.c          examples/ex_temp_control.c          examples/ex_position_servo.c

.PHONY: all clean test examples

all: $(OBJ)
	@echo "mini-lag-compensator: all objects built"

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

test: $(OBJ) $(TEST_SRC)
	$(CC) $(CFLAGS) $(TEST_SRC) $(OBJ) -o tests/test_runner $(LDFLAGS)
	./tests/test_runner

examples: $(OBJ)
	$(CC) $(CFLAGS) examples/ex_dc_motor.c $(OBJ) -o examples/ex_dc_motor $(LDFLAGS)
	$(CC) $(CFLAGS) examples/ex_temp_control.c $(OBJ) -o examples/ex_temp_control $(LDFLAGS)
	$(CC) $(CFLAGS) examples/ex_position_servo.c $(OBJ) -o examples/ex_position_servo $(LDFLAGS)

clean:
	rm -f src/*.o tests/test_runner examples/ex_dc_motor examples/ex_temp_control examples/ex_position_servo
