build:
	g++ -Wall -Wextra -Werror -pedantic -Wmissing-declarations -Wunreachable-code \
	-ftree-vectorize -ftree-vectorizer-verbose=2 \
	-std=c++20  \
	-o test \
	./test.cpp
	# -ggdb \
	# -O3 \
