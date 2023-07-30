all:
		gcc -shared -fpic -O3 -ffast-math -o verhulst2012/tridiag.so verhulst2012/cochlea_utils.c
		gcc -shared -fpic -O3 -ffast-math -o verhulst2015/tridiag.so verhulst2015/cochlea_utils.c
		gcc -shared -fpic -O3 -ffast-math -o verhulst2018/tridiag.so verhulst2018/cochlea_utils.c
	
.PHONY: clean	
clean:
		rm verhulst2012/tridiag.so
		rm verhulst2015/tridiag.so
		rm verhulst2018/tridiag.so
	   