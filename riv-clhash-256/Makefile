# export ASAN_SYMBOLIZER_PATH=/usr/bin/llvm-symbolizer-3.5
# export ASAN_OPTIONS=symbolize=1

CC = clang
CFLAGS = -W -Wall -Wextra -std=c99
CNIFLAGS = -W -Wall -Wextra -std=c99 -mavx2 -maes -mpclmul -march=native

NI = ni
REF = ref
TESTS = tests
RIV = riv

NI_OBJECTS = $(NI)/riv.c $(NI)/clhash.c
# $(patsubst %.c,%.o,$(wildcard $(NI)/*.c))
REF_OBJECTS = $(patsubst %.c,%.o,$(wildcard $(REF)/*.c))

ifdef DEBUG
CFLAGS += -O0 -g3 -DDEBUG -fsanitize=undefined -fsanitize=address -fsanitize=alignment -ftrapv -fno-omit-frame-pointer -fno-optimize-sibling-calls
CNIFLAGS += -O0 -g3 -DDEBUG -fsanitize=undefined -fsanitize=address -fsanitize=alignment -ftrapv -fno-omit-frame-pointer -fno-optimize-sibling-calls
else
CFLAGS += -O3
CNIFLAGS += -O3
endif

.PHONY: clean 
.PHONY: ni 
.PHONY: ni-test
.PHONY: ni-bench
.PHONY: ref-test
.PHONY: generatetests

all: ref ref-test

$(NI)/%.o: $(NI)/%.c
	$(CC) $(CNIFLAGS) -c $< -o $@

$(REF)/%.o: $(REF)/%.c
	$(CC) $(CFLAGS) -c $< -o $@

ni-bench: $(NI_OBJECTS)
	$(CC) $(CNIFLAGS) -I$(NI) $(TESTS)/benchmark.c $^ -o $(RIV)-$@

siv-bench: $(NI)/clhash.c $(NI)/siv.c
	$(CC) $(CNIFLAGS) -I$(NI) $(TESTS)/benchmark.c $^ -o $@

ni-test: $(NI_OBJECTS)
	$(CC) $(CNIFLAGS) -I$(NI) $(TESTS)/tests.c $^ -o $(RIV)-$@

ref-test: $(REF_OBJECTS)
	$(CC) $(CFLAGS) -I$(REF) $(TESTS)/tests.c $^  -o $(RIV)-$@

generatetests: $(REF_OBJECTS)
	$(CC) $(CFLAGS) -I$(REF) $(TESTS)/generatetests.c $(TESTS)/encrypt.c $^  -o $@

clean:
	@echo cleaning $(RIV)
	cd $(NI) && rm -f *.o *.gch *.out
	cd $(REF) && rm -f *.o *.gch *.out
	cd $(TESTS) && rm -f *.o *.gch *.out
	rm -f vgcore.* core
	rm -f *~ $(RIV)-ref-test $(RIV)-ni-test $(RIV)-ni-bench siv-bench generatetests

