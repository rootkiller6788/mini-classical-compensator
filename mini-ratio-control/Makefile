CC = gcc
CFLAGS = -std=c99 -Wall -Wextra -pedantic -O2 -g
LDFLAGS = -lm

SRC_DIR = src
INC_DIR = include
TEST_DIR = tests
EXAMPLE_DIR = examples
BUILD_DIR = build

SRCS = $(SRC_DIR)/ratio_control_core.c $(SRC_DIR)/ratio_station.c $(SRC_DIR)/blending_control.c $(SRC_DIR)/cross_limiting.c $(SRC_DIR)/ratio_feedforward.c $(SRC_DIR)/ratio_cascade.c
OBJS = $(SRCS:.c=.o)

TEST_SRC = $(TEST_DIR)/test_ratio.c
TEST_BIN = $(BUILD_DIR)/test_ratio

EXAMPLES = $(EXAMPLE_DIR)/example_blending $(EXAMPLE_DIR)/example_combustion $(EXAMPLE_DIR)/example_cascade_ratio

.PHONY: all clean test examples run-test check lines

all: test examples

%.o: %.c
	@mkdir -p $(dir $@)
	$(CC) $(CFLAGS) -I$(INC_DIR) -c $< -o $@

test: $(OBJS) $(TEST_SRC)
	@mkdir -p $(BUILD_DIR)
	$(CC) $(CFLAGS) -I$(INC_DIR) $(OBJS) $(TEST_SRC) -o $(TEST_BIN) $(LDFLAGS)
	./$(TEST_BIN)

run-test: test

examples: $(OBJS) $(EXAMPLE_DIR)/example_blending.c $(EXAMPLE_DIR)/example_combustion.c $(EXAMPLE_DIR)/example_cascade_ratio.c
	$(CC) $(CFLAGS) -I$(INC_DIR) $(OBJS) $(EXAMPLE_DIR)/example_blending.c -o $(EXAMPLE_DIR)/example_blending $(LDFLAGS)
	$(CC) $(CFLAGS) -I$(INC_DIR) $(OBJS) $(EXAMPLE_DIR)/example_combustion.c -o $(EXAMPLE_DIR)/example_combustion $(LDFLAGS)
	$(CC) $(CFLAGS) -I$(INC_DIR) $(OBJS) $(EXAMPLE_DIR)/example_cascade_ratio.c -o $(EXAMPLE_DIR)/example_cascade_ratio $(LDFLAGS)

lines:
	@echo "Headers:"
	@wc -l $(INC_DIR)/*.h
	@echo "Sources:"
	@wc -l $(SRC_DIR)/*.c
	@echo "Total include+src:"
	@cat $(INC_DIR)/*.h $(SRC_DIR)/*.c | wc -l

check: test run-test

clean:
	rm -rf $(BUILD_DIR) $(SRC_DIR)/*.o $(EXAMPLE_DIR)/example_blending $(EXAMPLE_DIR)/example_combustion $(EXAMPLE_DIR)/example_cascade_ratio
