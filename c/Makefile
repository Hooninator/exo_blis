OBJS	= blocking.o
SOURCE	= blocking.c
HEADER	= 
OUT	= mat
CC	 = gcc
FLAGS	 = -g -c -Wall
LFLAGS	 = 

all: $(OBJS)
	$(CC) -g $(OBJS) -o $(OUT) $(LFLAGS)

blocking.o: blocking.c
	$(CC) $(FLAGS) blocking.c 


clean:
	rm -f $(OBJS) $(OUT)