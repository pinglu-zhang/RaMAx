
BIN = ./bin
SRC = ./src
INC = ./include
OBJ = ./obj

COMMON = $(OBJ)/align.o \
	$(OBJ)/cluster.o \
	$(OBJ)/distancemat.o \
	$(OBJ)/options.o \
	$(OBJ)/util.o \
	$(OBJ)/sequence.o

QUICKTREEOBJ = $(COMMON) $(OBJ)/quicktree.o $(OBJ)/buildtree.o $(OBJ)/tree.o
SONLIB_QUICKTREEOBJ = $(OBJ)/buildtree.o \
	$(OBJ)/cluster.o \
	$(OBJ)/distancemat.o \
	$(OBJ)/options.o \
	$(OBJ)/sequence.o \
	$(OBJ)/tree.o \
	$(OBJ)/util.o

CFLAGS  += -c -I$(INC) -O3 -Wall

LDFLAGS += -g 

all: $(SONLIB_QUICKTREEOBJ)

quicktree : $(BIN)/quicktree

$(BIN)/quicktree : $(QUICKTREEOBJ)
	@mkdir -p $(dir $@)
	$(CC) $(DLFLAGS) -o $@ $(QUICKTREEOBJ) -lm

$(OBJ)/quicktree.o : $(SRC)/quicktree.c 
	@mkdir -p $(dir $@)
	$(CC) $(CFLAGS) $(CPPFLAGS) -o $(OBJ)/quicktree.o $(SRC)/quicktree.c

$(OBJ)/align.o : $(SRC)/align.c $(INC)/align.h
	@mkdir -p $(dir $@)
	$(CC) $(CFLAGS) $(CPPFLAGS) -o $(OBJ)/align.o $(SRC)/align.c

$(OBJ)/buildtree.o : $(SRC)/buildtree.c $(INC)/buildtree.h
	@mkdir -p $(dir $@)
	$(CC) $(CFLAGS) $(CPPFLAGS) -o $(OBJ)/buildtree.o $(SRC)/buildtree.c

$(OBJ)/cluster.o : $(SRC)/cluster.c $(INC)/cluster.h
	@mkdir -p $(dir $@)
	$(CC) $(CFLAGS) $(CPPFLAGS) -o $(OBJ)/cluster.o $(SRC)/cluster.c

$(OBJ)/distancemat.o : $(SRC)/distancemat.c $(INC)/distancemat.h
	@mkdir -p $(dir $@)
	$(CC) $(CFLAGS) $(CPPFLAGS) -o $(OBJ)/distancemat.o $(SRC)/distancemat.c

$(OBJ)/sequence.o : $(SRC)/sequence.c $(INC)/sequence.h
	@mkdir -p $(dir $@)
	$(CC) $(CFLAGS) $(CPPFLAGS) -o $(OBJ)/sequence.o $(SRC)/sequence.c

$(OBJ)/tree.o : $(SRC)/tree.c $(INC)/tree.h
	@mkdir -p $(dir $@)
	$(CC) $(CFLAGS) $(CPPFLAGS) -o $(OBJ)/tree.o $(SRC)/tree.c

$(OBJ)/util.o : $(SRC)/util.c $(INC)/util.h
	@mkdir -p $(dir $@)
	$(CC) $(CFLAGS) $(CPPFLAGS) -o $(OBJ)/util.o $(SRC)/util.c

$(OBJ)/options.o : $(SRC)/options.c $(INC)/options.h
	@mkdir -p $(dir $@)
	$(CC) $(CFLAGS) $(CPPFLAGS) -o $(OBJ)/options.o $(SRC)/options.c

clean:
	rm -rf $(OBJ)/*.o
