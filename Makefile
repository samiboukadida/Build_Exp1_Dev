# Variables de compilation
CXX = g++
CXXFLAGS = -m64 -fpermissive -g -Wall -Wextra -pedantic -Wno-unused-parameter
LDFLAGS = -lopengl32 -lglu32 -lglut
# Fichiers d'en-tête (headers)
#GLINS = -I/cygwin64/usr/include

# Fichiers sources et objets
SRCS = sfem.cpp exp1.tab.c
HEADERS = sfem.hpp RNM.hpp GC.hpp LexicoEdge.hpp HeapSort.hpp assertion.hpp gnuplotiso.hpp
OBJS = sfem.o exp1.o

# Nom de l'exécutable
TARGET = exp1

# Règle par défaut : compilation complète
all: $(TARGET)

# Cible pour l'exécutable
$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDFLAGS)

# Compilation de sfem.o
sfem.o: sfem.cpp sfem.hpp RNM.hpp LexicoEdge.hpp HeapSort.hpp assertion.hpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Compilation de exp1.o
exp1.o: exp1.tab.c sfem.hpp RNM.hpp GC.hpp LexicoEdge.hpp HeapSort.hpp gnuplotiso.hpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Génération des fichiers Bison
exp1.tab.c exp1.tab.h: exp1.y
	bison -d $<

# Nettoyage des fichiers temporaires
clean:
	rm -f $(OBJS) exp1.tab.c exp1.tab.h $(TARGET) *.o

# Nettoyage complet
realclean: clean
	rm -f *~ *.bak
