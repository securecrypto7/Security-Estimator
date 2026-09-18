# Makefile - Raccoon-based Interactive Aggregate Signatures (IAS)
# reference implementation.
#
#   make            build bench + test for all 6 parameter sets
#   make bench      build only the benchmark binaries
#   make test       build only the self-test binaries
#   make run        build and run every benchmark
#   make runtest    build and run every self-test
#   make clean
#
# Per-set binaries land in build/ as  bench_<SET>  and  test_<SET>.

CC      ?= gcc
CFLAGS  ?= -O2 -march=native -I.
LDLIBS  ?= -lm

SRC     = ias.c poly.c sample.c serialize.c fips202.c
BENCH   = bench.c $(SRC)
TEST    = test.c  $(SRC)

SETS    = IAS128_Q60 IAS128_Q30 IAS192_Q60 IAS192_Q30 IAS256_Q60 IAS256_Q30

BUILD   = build
BENCHBIN = $(addprefix $(BUILD)/bench_,$(SETS))
TESTBIN  = $(addprefix $(BUILD)/test_,$(SETS))

.PHONY: all bench test run runtest clean
all: bench test

bench: $(BENCHBIN)
test: $(TESTBIN)

$(BUILD):
	mkdir -p $(BUILD)

$(BUILD)/bench_%: $(BENCH) | $(BUILD)
	$(CC) $(CFLAGS) -DIAS_PARAM=$* $(BENCH) $(LDLIBS) -o $@

$(BUILD)/test_%: $(TEST) | $(BUILD)
	$(CC) $(CFLAGS) -DIAS_PARAM=$* $(TEST) $(LDLIBS) -o $@

run: bench
	@for b in $(BENCHBIN); do echo; $$b || exit 1; done

runtest: test
	@for t in $(TESTBIN); do echo; $$t || exit 1; done

clean:
	rm -rf $(BUILD)
