all: serialPageRank ompPageRank

serialPageRank:serialPageRank.c
	gcc serialPageRank.c -lm -O3 -o serialPageRank
	
ompPageRank:ompPageRank.c
	gcc ompPageRank.c -lm -O3 -o ompPageRank -fopenmp
	
clear:
	rm serialPageRank ompPageRank
	


