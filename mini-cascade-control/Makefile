CC       = gcc
CFLAGS   = -std=c11 -Wall -Wextra -pedantic -O2 -g
LDFLAGS  = -lm
INCLUDES = -Iinclude

SRCDIR   = src
OBJDIR   = build
TESTDIR  = tests
EXAMPLEDIR = examples

SRCS     = $(wildcard $(SRCDIR)/*.c)
OBJS     = $(patsubst $(SRCDIR)/%.c, $(OBJDIR)/%.o, $(SRCS))
TESTS    = $(wildcard $(TESTDIR)/*.c)
TESTBINS = $(patsubst $(TESTDIR)/%.c, $(OBJDIR)/test_%, $(TESTS))
EXAMPLES = $(wildcard $(EXAMPLEDIR)/*.c)
EXAMPLEBINS = $(patsubst $(EXAMPLEDIR)/%.c, $(OBJDIR)/example_%, $(EXAMPLES))

.PHONY: all clean test examples

all: $(OBJS)
	@echo "All objects compiled (cascade control). Running: make test"

$(OBJDIR)/%.o: $(SRCDIR)/%.c
	@mkdir -p $(OBJDIR)
	$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@

$(OBJDIR)/test_%: $(TESTDIR)/%.c $(OBJS)
	@mkdir -p $(OBJDIR)
	$(CC) $(CFLAGS) $(INCLUDES) $< $(OBJS) -o $@ $(LDFLAGS)

$(OBJDIR)/example_%: $(EXAMPLEDIR)/%.c $(OBJS)
	@mkdir -p $(OBJDIR)
	$(CC) $(CFLAGS) $(INCLUDES) $< $(OBJS) -o $@ $(LDFLAGS)

test: $(TESTBINS)
	@for t in $(TESTBINS); do echo "Running $$t..."; ./$$t || exit 1; done
	@echo "All tests passed."

examples: $(EXAMPLEBINS)
	@for e in $(EXAMPLEBINS); do echo "Running $$e..."; ./$$e; done

clean:
	rm -rf $(OBJDIR)
