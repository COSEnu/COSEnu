include Makefile.inc

TARGET= main
OBJS  = main.o

${TARGET}: ${OBJS}
	$(CXX) $(OPT) $< -o $@ $(LIBS)
	


run:${TARGET}
	./${TARGET}

clean:
	rm ${TARGET} ${OBJS}
