# Makefile — mini-pid-tuning-ziegler
# PID Tuning by Ziegler-Nichols and Extended Methods
#
# Target: make test → compile and run all tests
#         make examples → build all examples
#         make clean  → remove build artifacts

CC       = gcc
CFLAGS   = -Wall -Wextra -O2 -std=c99 -Iinclude -Isrc
LDFLAGS  = -lm

SRC_DIR  = src
INC_DIR  = include
TST_DIR  = tests
EXM_DIR  = examples
BLD_DIR  = build

SRCS     = $(wildcard $(SRC_DIR)/*.c)
TSTS     = $(wildcard $(TST_DIR)/*.c)
EXMS     = $(wildcard $(EXM_DIR)/*.c)

OBJS     = $(patsubst $(SRC_DIR)/%.c, $(BLD_DIR)/%.o, $(SRCS))
TEST_BIN = $(BLD_DIR)/test_pid_tuning
EXM_BINS = $(patsubst $(EXM_DIR)/%.c, $(BLD_DIR)/%, $(EXMS))

.PHONY: all test examples clean

all: test examples

$(BLD_DIR):
	mkdir -p $(BLD_DIR)

$(BLD_DIR)/%.o: $(SRC_DIR)/%.c | $(BLD_DIR)
	$(CC) $(CFLAGS) -c $< -o $@

$(TEST_BIN): $(OBJS) $(TSTS) | $(BLD_DIR)
	$(CC) $(CFLAGS) $(TST_DIR)/test_pid_tuning.c $(OBJS) -o $@ $(LDFLAGS)

$(BLD_DIR)/%: $(EXM_DIR)/%.c $(OBJS) | $(BLD_DIR)
	$(CC) $(CFLAGS) $< $(OBJS) -o $@ $(LDFLAGS)

test: $(TEST_BIN)
	@echo "Running PID Tuning test suite..."
	@$(TEST_BIN)

examples: $(EXM_BINS)
	@echo "All examples built."
	@for ex in $(EXM_BINS); do \
		echo "--- Running $$ex ---"; \
		$$ex; \
		echo ""; \
	done

clean:
	rm -rf $(BLD_DIR)
