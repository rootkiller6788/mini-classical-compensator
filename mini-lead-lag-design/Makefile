# Makefile for mini-lead-lag-design
# Classical Compensator Design ¡ª Lead, Lag, Lead-Lag
CC = gcc
CFLAGS = -Wall -Wextra -std=c11 -O2 -Iinclude
LDFLAGS = -lm

SRC = src/lead_design.c       src/lag_design.c       src/lead_lag_design.c       src/compensator_spec.c       src/compensator_freq_design.c       src/compensator_rlocus_design.c       src/compensator_digital.c       src/compensator_sensitivity.c       src/compensator_analytical.c       src/compensator_loopshaping.c

OBJ = $(SRC:.c=.o)

TEST_SRC = tests/test_lead_lag.c
EX1 = examples/ex_dc_motor_lead.c
EX2 = examples/ex_temp_control_lag.c
EX3 = examples/ex_servo_lead_lag.c
EX4 = examples/ex_antenna_tracking.c

.PHONY: all clean test examples

all: $(OBJ)
	@echo "mini-lead-lag-design: all objects built"

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

test: $(OBJ) $(TEST_SRC)
	$(CC) $(CFLAGS) $(TEST_SRC) $(OBJ) -o tests/test_runner $(LDFLAGS)
	./tests/test_runner

examples: $(OBJ) $(EX1) $(EX2) $(EX3) $(EX4)
	$(CC) $(CFLAGS) $(EX1) $(OBJ) -o examples/ex_dc_motor_lead $(LDFLAGS)
	$(CC) $(CFLAGS) $(EX2) $(OBJ) -o examples/ex_temp_control_lag $(LDFLAGS)
	$(CC) $(CFLAGS) $(EX3) $(OBJ) -o examples/ex_servo_lead_lag $(LDFLAGS)
	$(CC) $(CFLAGS) $(EX4) $(OBJ) -o examples/ex_antenna_tracking $(LDFLAGS)

clean:
	rm -f src/*.o tests/test_runner examples/ex_dc_motor_lead examples/ex_temp_control_lag examples/ex_servo_lead_lag examples/ex_antenna_tracking
