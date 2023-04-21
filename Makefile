
CC=mpicxx
CFLAGS=-pg -O3


TARGET=main
SRCS = main.cpp \
	./lu.cpp \
	./mesh.cpp\
	./element.cpp\
	./point.cpp\
	./function.cpp\
	./bc.cpp\
	./solver.cpp\
	./assembling.cpp



INC = -I./include
OBJS = $(SRCS:.cpp=.o)
$(TARGET):$(OBJS)
#	@echo TARGET:$@
#	@echo OBJECTS:$^
	$(CC) -o $@ $^ $(CFLAGS)

clean:
	rm -rf $(TARGET) $(OBJS)

%.o:%.cpp
	$(CC) $(CFLAGS) $(INC) -o $@ -c $<
