CC = g++
# CFLAGS = -Wall -Wextra -Werror -std=c11 -g -fsanitize=leak -fsanitize=address
CFLAGS = -Wall -Wextra -Werror -std=c++17 -g #-fsanitize=address

OS:=$(shell uname -s)
ifeq ($(OS), Darwin)
LDFLAGS = -lcheck -lm
endif
ifeq ($(OS),Linux)
LDFLAGS=-lcheck_pic -lrt -lpthread -lsubunit -lm
endif

all: clean
	$(CC) $(CFLAGS) s21_matrix_oop.cc
	./a.out

s21_matrix_oop.a: clean
	$(CC) $(CFLAGS) -c s21_matrix_oop.cc 
	ar rcs s21_matrix_oop.a *.o
	ranlib s21_matrix_oop.a

test: clean s21_matrix_oop.a
	# g++ tests.cc -lgtest && ./a.out
	$(CC) tests.cc -lgtest  -L. s21_matrix_oop.a
	./a.out

gcov_report: clean s21_matrix_oop.a
	# $(CC) $(CFLAGS) --coverage tests.cc s21_matrix_oop.cc s21_matrix_oop.a $(LDFLAGS) -o gcov_tests
	# ./gcov_tests
	# lcov -t "s21_matrix_oop" -o s21_matrix_oop.info -c -d . --rc lcov_branch_coverage=1
	# genhtml -o report s21_matrix_oop.info --rc lcov_branch_coverage=1
	# open report/index.html

	$(CC) $(CFLAGS)  -fprofile-arcs -ftest-coverage s21_matrix_oop.cc tests.cc -lgtest -o test
	./test
	gcov --no-output tests.cc
	lcov --no-external --capture --directory . --output-file coverage.info --rc lcov_branch_coverage=1
	genhtml coverage.info -o report --rc lcov_branch_coverage=1
	open ./report/index.html

style_check:
	cppcheck  *.h *.cc  --language=c++
	cp ../materials/linters/.clang-format .clang-format
	clang-format -n tests/*.cc *.cc *h 
	rm -rf .clang-format

style_make:
	cppcheck  *.h *.cc  --language=c++
	cp ../materials/linters/.clang-format .clang-format
	clang-format -i tests/*.cc *cc *.h исправить
	rm -rf .clang-format

valgrind: test #LINUX only
	valgrind -s --leak-check=full --show-leak-kinds=all  ./a.out

clean:
	rm -rf vgcheck ./*.o ./*.a ./*.out ./*.dSYM gcov_test *.html *.css ./GcovReport ./test ./*.gcda ./*.gch ./*.gcno ./report ./calc_lib_test
	rm -rf ./*.info ./gcov_tests ./gcov_tests ./res_leak.txt ./s21_matrix_test ./*.LcdoKh
	rm -rf ./.clang-format

p: 
	git add .
	git commit -m "final"
	git push origin develop 