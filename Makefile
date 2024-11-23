# ============================
# Makefile pour projet exp1
# ============================

# Pour le test:
OUTPUT_DIR = outputs
DATA_FILE = data.txt
GNUPLOT_SCRIPT = ../scripts/uh3_3D.gnu

# Variables de compilation
CXX = g++
CXXFLAGS = -m64 -fpermissive -g -Wall -Wextra -pedantic -Wno-unused-parameter  -DYYDEBUG
LDFLAGS = -lopengl32 -lglu32 -lglut

# Fichiers d'en-tête (headers)
HEADERS = sfem.hpp RNM.hpp GC.hpp LexicoEdge.hpp HeapSort.hpp assertion.hpp gnuplotiso.hpp

# Fichiers sources
SRCS = sfem.cpp exp1.tab.c

# Fichiers objets
OBJS = $(SRCS:.cpp=.o)

# Nom de l'exécutable
TARGET = exp1

# Règle par défaut : compilation complète
all: $(TARGET)

# Règle pour créer l'exécutable
$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDFLAGS)

# Règle générique pour compiler les fichiers objets
%.o: %.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Règle pour générer les fichiers Bison
exp1.tab.c exp1.tab.h: exp1.y
	bison -d --debug $<

# Nettoyage des fichiers temporaires
clean:
	rm -f $(OBJS) exp1.tab.c exp1.tab.h $(TARGET)

# Nettoyage complet (inclut les fichiers sauvegardés)
realclean: clean
	rm -f *~ *.bak

# Test (exécute l'application après compilation)
test: $(TARGET)
	mkdir -p $(OUTPUT_DIR)
	cd $(OUTPUT_DIR) && \
	../$(TARGET) ../$(DATA_FILE) > q11.out 2>&1 && \
	([ -f uh3 ] && gnuplot $(GNUPLOT_SCRIPT) || echo "Erreur : uh3 n'existe pas.") && \
	([ -f uh3 ] && gnuplot -e "splot 'uh3' using 1:2:3 with lines; pause -1" || echo "Erreur : uh3 n'existe pas pour affichage.")
