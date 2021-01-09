CC=g++
ODIR=build
IDIR=build/include
LDIR=build/lib
CFLAGS = -std=gnu++17 -O2 -I$(IDIR)


all:
	# Create build directories
	mkdir -p $(ODIR)
	mkdir -p $(IDIR)
	mkdir -p $(LDIR)

	# Copy headers
	cp core/*.hpp $(ODIR)/include

	# Build core libraries
	$(CC) -fPIC -shared -rdynamic core/ruben.cpp -o $(ODIR)/lib/libruben.so -lquadmath -O2 
	$(CC) -fPIC -shared -rdynamic core/davies.cpp -o $(ODIR)/lib/libdavies.so -lquadmath -O2 
	
	$(CC) -fPIC core/davies.cpp core/ruben.cpp core/wchissum.cpp -shared -rdynamic -o $(ODIR)/lib/libwchissum.so -lquadmath $(CFLAGS)
#	$(CC) -fPIC -shared core/wchissum.cpp -o $(ODIR)/lib/libwchissum.so $(CFLAGS) -L$(LDIR) -lruben -ldavies -lquadmath

#	make test


test:
	# Build tests
	$(CC) -o $(ODIR)/test_ruben.o tests/test_ruben.cpp $(CFLAGS) -L$(LDIR) -lruben -lquadmath
	./build/test_ruben.o
	$(CC) -o $(ODIR)/test_davies.o tests/test_davies.cpp $(CFLAGS) -L$(LDIR) -ldavies -lquadmath	
	./build/test_davies.o

	#$(CC) -o $(ODIR)/test_wchissum.o tests/test_wchissum.cpp $(CFLAGS) -L$(LDIR) -lwchissum -lquadmath
	#./build/test_wchissum.o

clean:
	rm -f $(ODIR)/*
	rm -f $(IDIR)/*
	rm -f $(LDIR)/*
