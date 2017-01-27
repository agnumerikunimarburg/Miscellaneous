CXX = g++
CXXFLAGS = -O3 -Wall -pipe 
LDFLAGS = 

OBJ = subdivision.o 


EXE = $(OBJ:.o=)


all:: compile

compile:: $(EXE)

$(OBJ): %.o: %.cpp
	$(CXX) $(CXXFLAGS) -c -DNDEBUG -o $@ $<

$(EXE): %: %.o
	$(CXX) $(LDFLAGS) $< -o $@

.PHONY: clean

clean::
	rm -f $(OBJ) $(EXE)
	rm -f *~

