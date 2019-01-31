TARGETS = main.x

CPPFLAGS=-Wall -Wextra -pedantic -Wno-sign-compare --std=c++11 -O3 -pthread

all: $(TARGETS)

$(TARGETS): %.x: %.cpp
	g++ $(CPPFLAGS) $< -o $@
	
clean:
	rm -f $(TARGETS)

