NAME := s21_matrix_oop.a
CC = gcc
CPPFLAGS = -Wall -Werror -Wextra -std=c++17 -g
TEST_LIBS = -lgtest -lstdc++ -lgtest_main -lpthread
SRC := S21_Matrix_oop.cc
SRC_TEST := test.cc

all: s21_matrix_oop.a

build: S21_Matrix_oop.cc
	$(CC) -c $(CPPFLAGS) -lstdc++ -std=c++17 S21_Matrix_oop.cc -o S21_Matrix_oop.o

s21_matrix_oop.a: build
	ar rcs s21_matrix_oop.a S21_Matrix_oop.o
	ranlib s21_matrix_oop.a

clean:
	rm -rf *.o *.a *.out *.so *.gcda *.gcno *.gch rep.info *.html *.css calc *.dSYM test *.info report/ .DS_Store ../.DS_Store
	rm -rf *gz

test: s21_matrix_oop.a test.cc
	$(CC) $(TEST_LIBS) $(CPPFLAGS) s21_matrix_oop.a test.cc -o test
	./test

gcov_report: $(NAME)
	$(CC) $(CPPFLAGS) --coverage $(SRC) $(SRC_TEST) $(TEST_LIBS) -o gcov_test
	chmod +x *
	./gcov_test
	lcov -t "gcov_test" -o gcov_test.info --no-external -c -d .
	genhtml -o report/ gcov_test.info
	open ./report/index.html
	rm -rf *.gcno *gcda *.gcov
	
check: s21_matrix_oop.a
	cppcheck *.h *.cc
	mv ../materials/linters/.clang-format .clang-format
	clang-format -i *.cc *.h
	clang-format -n *.cc *.h
	mv .clang-format ../materials/linters/.clang-format
	CK_FORK=no leaks --atExit -- ./test