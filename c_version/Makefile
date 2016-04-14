RUN := run_TTVFast
SOURCES := $(wildcard *.c)
HEADERS := $(wildcard *.h)

all: test

$(RUN): $(SOURCES) $(HEADERS)
	gcc -o $@ run_TTVFast.c TTVFast.c -lm -O3  -Wall -Wextra

test: $(RUN)
	./run_TTVFast setup_file Times RV_file RV_out
	diff Times Times.expected
	diff RV_out RV_out.expected


clean:
	@-rm $(RUN)

.PHONY: test
