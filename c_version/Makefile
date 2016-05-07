RUN := run_TTVFast
COMPARE := compare_output
SOURCES := $(wildcard *.c)
HEADERS := $(wildcard *.h)

all: $(RUN)

$(RUN): $(SOURCES) $(HEADERS)
	gcc -o $@ run_TTVFast.c TTVFast.c -lm -O3  -Wall -Wextra

test: $(COMPARE) $(RUN)
	./run_TTVFast setup_file Times RV_file RV_out
	./$(COMPARE)

$(COMPARE): compare_output.c
	gcc -o $@ $< -lm -Wall -Wextra


clean:
	@-rm $(RUN)

.PHONY: test
