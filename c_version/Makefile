RUN := run_TTVFast
SOURCES := $(wildcard *.c)
HEADERS := $(wildcard *.h)

all: $(RUN)

$(RUN): $(SOURCES) $(HEADERS)
	gcc -o $@ run_TTVFast.c TTVFast.c -lm -O3  -Wall -Wextra

test:
	./run_TTVFast setup_file Times RV_file RV_out

clean:
	@-rm $(RUN)

.PHONY: test
